"""Manage the gene-set collections used by iGEA (command line).

    python code/libraries.py list
    python code/libraries.py activate "C5: Human Phenotype Ontology" "C7: VAX: Vaccine response"
    python code/libraries.py deactivate "C2: BioCarta"
    python code/libraries.py add my_sets.gmt --name "Custom: My sets" [--inactive] [--color "#1F77B4"]
    python code/libraries.py remove "Custom: My sets"
    python code/libraries.py update-gene-data --gene-info Homo_sapiens.gene_info.gz --gene-history gene_history.gz
    python code/libraries.py update-msigdb --zip msigdb_v2026.1.Hs_files_to_download_locally.zip [--dry-run]

Collections are listed in data/libraries/alias.json (name, file, active, colour). The command-line tool and the
app use the active ones by default; the CLI can also take any GMT file with -l.

Gene symbols. Every collection is harmonised with the NCBI gene data in data/recent_release: each symbol is
mapped to its current official symbol (synonyms and gene history); symbols that cannot be mapped are kept
unchanged. This is how the bundled MSigDB 2025.1 files were prepared, so input lists, backgrounds and gene sets
use the same symbols. Update the NCBI gene data first (update-gene-data), then the collections (update-msigdb).

Nothing is deleted: replaced or removed GMT files are moved to data/libraries/previous/.
Changing collections or gene data changes analysis results; record the versions you use.
"""
import argparse
import gzip
import io
import json
import re
import shutil
import sys
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MSIGDB_FILE = re.compile(r"^(?P<prefix>.+)\.v(?P<version>\d{4}\.\d+)\.Hs\.symbols\.gmt$")


# ----------------------------------------------------------------------------- alias.json
def lib_dir(data: Path) -> Path:
    return data / "libraries"


def load_aliases(data: Path):
    p = lib_dir(data) / "alias.json"
    return json.load(open(p)) if p.exists() else []


def save_aliases(data: Path, aliases) -> None:
    p = lib_dir(data) / "alias.json"
    tmp = p.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(aliases, f, indent=4)
        f.write("\n")
    tmp.replace(p)


def find_entries(aliases, keys):
    """Entries matching each key (exact name or file name; otherwise a unique case-insensitive name match)."""
    out = []
    for k in keys:
        hit = [a for a in aliases if a["name"] == k or a["file"] == k]
        if not hit:
            hit = [a for a in aliases if k.lower() in a["name"].lower()]
        if len(hit) != 1:
            names = ", ".join(a["name"] for a in hit) or "none"
            sys.exit(f'"{k}" matches {len(hit)} collections ({names}); use the exact name from `list`.')
        out.append(hit[0])
    return out


def read_gmt(path: Path):
    sets = []
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\n\r")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError(f"{path.name}, line {n}: expected name<TAB>description<TAB>gene1<TAB>gene2...")
            sets.append((parts[0], parts[1], [g for g in parts[2:] if g.strip()]))
    return sets


def write_gmt(path: Path, sets) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for name, desc, genes in sets:
            f.write("\t".join([name, desc] + genes) + "\n")


def move_to_previous(data: Path, filename: str, label: str) -> None:
    src = lib_dir(data) / filename
    if not src.exists():
        return
    dst_dir = lib_dir(data) / "previous" / label
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / filename
    if dst.exists():
        dst = dst_dir / f"{src.stem}.{datetime.now():%Y%m%d_%H%M%S}{src.suffix}"
    shutil.move(str(src), str(dst))


# ----------------------------------------------------------------------------- gene symbols
def converter(data: Path):
    sys.path.insert(0, str(ROOT / "code"))
    from gene_converter import GeneConverter
    rr = data / "recent_release"
    info = next((rr / n for n in ("Homo_sapiens.gene_info", "Homo_sapiens.gene_info.gz") if (rr / n).exists()), None)
    hist = next((rr / n for n in ("gene_history", "gene_history_9606.gz", "gene_history.gz") if (rr / n).exists()), None)
    if info is None:
        sys.exit(f"NCBI gene info not found in {rr}; run update-gene-data first.")
    return GeneConverter(str(info), str(hist) if hist else None)


def harmonise(sets, conv):
    """Map each symbol to its current official symbol; keep symbols that cannot be mapped."""
    out, renamed, unknown = [], 0, 0
    for name, desc, genes in sets:
        new = []
        for g in genes:
            m = conv.validate_and_map_symbol(g)
            if m is None:
                unknown += 1
                new.append(g)
            else:
                renamed += m != g
                new.append(m)
        out.append((name, desc, new))
    return out, renamed, unknown


# ----------------------------------------------------------------------------- commands
def cmd_list(a):
    aliases = load_aliases(a.data)
    if not aliases:
        print("No collections (data/libraries/alias.json is empty or missing).")
        return
    w = max(len(x["name"]) for x in aliases)
    print(f"{'':2}{'Collection':{w}}  {'Sets':>6}  File")
    for x in aliases:
        p = lib_dir(a.data) / x["file"]
        n = sum(1 for line in open(p, encoding="utf-8") if line.strip()) if p.exists() else "missing"
        print(f"{'*' if x.get('active') else ' ':2}{x['name']:{w}}  {n:>6}  {x['file']}")
    print(f"\n* = active ({sum(bool(x.get('active')) for x in aliases)} of {len(aliases)}); "
          "active collections are used by default in the CLI and the app.")


def set_active(a, value):
    aliases = load_aliases(a.data)
    for e in find_entries(aliases, a.names):
        e["active"] = value
        print(f"{'Activated' if value else 'Deactivated'}: {e['name']}")
    save_aliases(a.data, aliases)


def cmd_add(a):
    src = Path(a.gmt)
    sets = read_gmt(src)
    aliases = load_aliases(a.data)
    if any(x["name"] == a.name for x in aliases):
        sys.exit(f'A collection named "{a.name}" already exists.')
    dest = lib_dir(a.data) / src.name
    if dest.exists():
        sys.exit(f"{dest.name} already exists in {lib_dir(a.data)}; rename the file first.")
    note = "symbols kept as given"
    if not a.no_harmonise:
        sets, renamed, unknown = harmonise(sets, converter(a.data))
        note = f"{renamed} symbols updated to current NCBI symbols, {unknown} not recognised (kept)"
    write_gmt(dest, sets)
    entry = {"name": a.name, "file": dest.name, "active": not a.inactive}
    if a.color:
        entry["color"] = a.color
    aliases.append(entry)
    save_aliases(a.data, aliases)
    sizes = [len(set(g)) for _, _, g in sets]
    in_range = sum(10 <= s <= 600 for s in sizes)
    print(f'Added "{a.name}" ({dest.name}): {len(sets)} gene sets, {in_range} with 10-600 genes '
          f"(default size filter); {note}; {'active' if entry['active'] else 'inactive'}.")


def cmd_remove(a):
    aliases = load_aliases(a.data)
    for e in find_entries(aliases, a.names):
        aliases.remove(e)
        move_to_previous(a.data, e["file"], "removed")
        print(f"Removed: {e['name']} (file moved to data/libraries/previous/removed/)")
    save_aliases(a.data, aliases)


def cmd_update_gene_data(a):
    rr = a.data / "recent_release"
    rr.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    old = rr / "previous" / stamp
    for name in ("Homo_sapiens.gene_info", "Homo_sapiens.gene_info.gz", "gene_history", "gene_history_9606.gz"):
        if (rr / name).exists():
            old.mkdir(parents=True, exist_ok=True)
            shutil.move(str(rr / name), str(old / name))
    opener = gzip.open if str(a.gene_info).endswith(".gz") else open
    with opener(a.gene_info, "rb") as fi, gzip.open(rr / "Homo_sapiens.gene_info.gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    kept = 0
    opener = gzip.open if str(a.gene_history).endswith(".gz") else open
    with opener(a.gene_history, "rt", encoding="utf-8") as fi, gzip.open(rr / "gene_history_9606.gz", "wt", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith("#") or line.startswith("9606\t"):
                fo.write(line)
                kept += not line.startswith("#")
    print(f"NCBI gene data updated in {rr} (previous files in {old if old.exists() else 'n/a'}); "
          f"gene history: {kept:,} human records.")
    print("Next: harmonise the collections with the new symbols, e.g. update-msigdb with your MSigDB zip.")


def cmd_update_msigdb(a):
    z = zipfile.ZipFile(a.zip)
    in_zip = {}
    for n in z.namelist():
        m = MSIGDB_FILE.match(Path(n).name)
        if m:
            in_zip[m["prefix"]] = (n, m["version"])
    if not in_zip:
        sys.exit("No *.vYYYY.N.Hs.symbols.gmt files found in the zip.")
    new_version = sorted({v for _, v in in_zip.values()})[-1]
    aliases = load_aliases(a.data)
    targets = []
    for e in aliases:
        m = MSIGDB_FILE.match(e["file"])
        if not m:
            continue                                   # not an MSigDB file (e.g. custom collections)
        if a.only and m["prefix"] not in a.only and e["name"] not in a.only:
            continue
        if a.active_only and not e.get("active"):
            continue
        targets.append((e, m["prefix"], m["version"]))
    conv = None if a.no_harmonise else converter(a.data)
    print(f"MSigDB {new_version} from {Path(a.zip).name}: {len(targets)} collections to update"
          + (" (dry run, nothing written)" if a.dry_run else ""))
    for e, prefix, old_version in targets:
        if prefix not in in_zip:
            print(f"  - {e['name']}: {prefix} not in this release; kept at v{old_version}")
            continue
        member, version = in_zip[prefix]
        new_sets = read_gmt_bytes(z.read(member), Path(member).name)
        renamed = unknown = 0
        if conv:
            new_sets, renamed, unknown = harmonise(new_sets, conv)
        old_path = lib_dir(a.data) / e["file"]
        old_names = {s[0] for s in read_gmt(old_path)} if old_path.exists() else set()
        new_names = {s[0] for s in new_sets}
        print(f"  - {e['name']}: v{old_version} -> v{version}, {len(old_names)} -> {len(new_names)} sets "
              f"(+{len(new_names - old_names)} / -{len(old_names - new_names)})"
              + (f"; {renamed} symbols updated, {unknown} not recognised (kept)" if conv else ""))
        if a.dry_run:
            continue
        new_file = f"{prefix}.v{version}.Hs.symbols.gmt"
        write_gmt(lib_dir(a.data) / new_file, new_sets)
        if new_file != e["file"]:
            move_to_previous(a.data, e["file"], f"v{old_version}")
        e["file"] = new_file
    if not a.dry_run:
        save_aliases(a.data, aliases)
        print(f"Done. Previous files are in data/libraries/previous/. Record the MSigDB version ({new_version}) "
              "and the NCBI gene data you used with your results.")


def read_gmt_bytes(raw: bytes, label: str):
    sets = []
    for n, line in enumerate(io.StringIO(raw.decode("utf-8")), 1):
        line = line.rstrip("\n\r")
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            raise ValueError(f"{label}, line {n}: not a GMT line")
        sets.append((parts[0], parts[1], [g for g in parts[2:] if g.strip()]))
    return sets


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", type=Path, default=ROOT / "data", help="data folder (default: the repository's data/)")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list collections; * marks active ones")
    for name, helptext in (("activate", "use these collections by default"), ("deactivate", "stop using them by default")):
        s = sub.add_parser(name, help=helptext)
        s.add_argument("names", nargs="+", help="collection names (as in `list`) or file names")
    s = sub.add_parser("add", help="add a custom GMT collection")
    s.add_argument("gmt")
    s.add_argument("--name", required=True, help='display name, e.g. "Custom: My sets"')
    s.add_argument("--inactive", action="store_true", help="add without activating")
    s.add_argument("--color", help="colour for the network figures, e.g. #1F77B4")
    s.add_argument("--no-harmonise", action="store_true", help="keep gene symbols exactly as in the file")
    s = sub.add_parser("remove", help="remove collections (files are moved to data/libraries/previous/removed/)")
    s.add_argument("names", nargs="+")
    s = sub.add_parser("update-gene-data", help="install new NCBI gene info and gene history files")
    s.add_argument("--gene-info", required=True, type=Path, help="Homo_sapiens.gene_info(.gz) from NCBI")
    s.add_argument("--gene-history", required=True, type=Path, help="gene_history(.gz) from NCBI (filtered to human)")
    s = sub.add_parser("update-msigdb", help="update MSigDB collections from the downloaded release zip")
    s.add_argument("--zip", required=True, type=Path, help="msigdb_vYYYY.N.Hs_files_to_download_locally.zip")
    s.add_argument("--only", nargs="+", help="update only these (file prefixes such as c2.cp.reactome, or names)")
    s.add_argument("--active-only", action="store_true", help="update only active collections")
    s.add_argument("--no-harmonise", action="store_true", help="keep MSigDB symbols as released")
    s.add_argument("--dry-run", action="store_true", help="report what would change without writing")
    a = p.parse_args(argv)
    {"list": cmd_list, "activate": lambda x: set_active(x, True), "deactivate": lambda x: set_active(x, False),
     "add": cmd_add, "remove": cmd_remove, "update-gene-data": cmd_update_gene_data,
     "update-msigdb": cmd_update_msigdb}[a.cmd](a)


if __name__ == "__main__":
    main()

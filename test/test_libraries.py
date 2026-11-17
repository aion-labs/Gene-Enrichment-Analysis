import json, shutil, sys, zipfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
import libraries as L


def _data(tmp):
    d = tmp / "data"
    (d / "libraries").mkdir(parents=True)
    (d / "recent_release").mkdir()
    for f in ("Homo_sapiens.gene_info.gz", "gene_history_9606.gz"):
        shutil.copy(ROOT / "data" / "recent_release" / f, d / "recent_release" / f)
    (d / "libraries" / "h.all.v2025.1.Hs.symbols.gmt").write_text("SET_A\td\tTP53\tMDM2\tSLC22A18\nSET_B\td\tEGFR\tAKT1\n")
    json.dump([{"name": "H: Hallmark", "file": "h.all.v2025.1.Hs.symbols.gmt", "active": True}],
              open(d / "libraries" / "alias.json", "w"))
    return d


def test_add_activate_remove(tmp_path):
    d = _data(tmp_path)
    gmt = tmp_path / "mine.gmt"
    gmt.write_text("MINE\td\tTP53\tSLC22A18\tNOTAGENE\n")
    L.main(["--data", str(d), "add", str(gmt), "--name", "Custom: Mine", "--inactive"])
    a = L.load_aliases(d)
    assert a[-1] == {"name": "Custom: Mine", "file": "mine.gmt", "active": False}
    assert (d / "libraries" / "mine.gmt").read_text().split("\t")[3:] == ["SLC67A1", "NOTAGENE\n"]  # renamed, unknown kept
    L.main(["--data", str(d), "activate", "Custom: Mine"])
    assert L.load_aliases(d)[-1]["active"] is True
    L.main(["--data", str(d), "remove", "Custom: Mine"])
    assert [x["name"] for x in L.load_aliases(d)] == ["H: Hallmark"]
    assert (d / "libraries" / "previous" / "removed" / "mine.gmt").exists()


def test_update_msigdb(tmp_path):
    d = _data(tmp_path)
    z = tmp_path / "m.zip"
    with zipfile.ZipFile(z, "w") as f:
        f.writestr("x/h.all.v2026.1.Hs.symbols.gmt", "SET_A\td\tTP53\tSLC22A18\nSET_C\td\tKRAS\tNRAS\n")
        f.writestr("x/c9.all.v2026.1.Hs.symbols.gmt", "X\td\tA1BG\n")
    L.main(["--data", str(d), "update-msigdb", "--zip", str(z), "--dry-run"])
    assert L.load_aliases(d)[0]["file"] == "h.all.v2025.1.Hs.symbols.gmt"
    L.main(["--data", str(d), "update-msigdb", "--zip", str(z)])
    e = L.load_aliases(d)[0]
    assert e["file"] == "h.all.v2026.1.Hs.symbols.gmt" and e["name"] == "H: Hallmark" and e["active"]
    assert "SLC67A1" in (d / "libraries" / e["file"]).read_text()
    assert (d / "libraries" / "previous" / "v2025.1" / "h.all.v2025.1.Hs.symbols.gmt").exists()

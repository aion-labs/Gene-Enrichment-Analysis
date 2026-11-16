"""Send the iGEA LLM prompt to a language model (OpenAI API) and save the response.

Standard library only, so it also runs outside the iGEA environment:

    export OPENAI_API_KEY="sk-..."        # never put the key in a file or on the command line
    python3 code/llm_client.py <run folder or llm_prompt.txt> --model <model-id> [--repeats 3]
    python3 code/llm_client.py --list-models
    python3 code/llm_client.py <run folder> --model <model-id> --context "siRNA screen for ..."   # second turn

The CLI calls the same function when ``--llm-model`` is given. Each call is independent (no
conversation, no stored memory): the request contains only the prompt, and responses are not stored
by the provider (store=false). With --context, a second turn is sent for each saved response: the
first prompt, that response and a follow-up revealing how the gene list was generated (llm_export.FOLLOWUP),
sent together in one stateless request; results are saved as followup_r<N>. For every call the folder
``llm/`` receives the full request and response (JSON), the response text (Markdown) and a line in
``manifest.tsv`` (requested and returned model, response id, time, tokens).
"""
import argparse
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

API = "https://api.openai.com/v1"


def _request(method, path, key, body=None, timeout=900):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(API + path, data=data, method=method,
                                 headers={"Authorization": "Bearer " + key, "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def text_of(resp):
    parts = []
    for item in resp.get("output", []):
        if item.get("type") == "message":
            for c in item.get("content", []):
                if c.get("type") in ("output_text", "text"):
                    parts.append(c.get("text", ""))
    if not parts and "choices" in resp:  # chat completions
        parts = [resp["choices"][0]["message"]["content"]]
    return "\n".join(parts).strip()


class ProviderError(RuntimeError):
    """The provider refused or failed the request (HTTP error after retries)."""
    def __init__(self, code, message):
        super().__init__(f"HTTP {code}: {message}")
        self.code, self.message = code, message


POLICY_CODES = ("policy", "flagged", "safety")


def _is_refusal(err) -> bool:
    """A 400 caused by a provider content policy (not a malformed request)."""
    return err.code == 400 and any(k in err.message.lower() for k in POLICY_CODES)


def _record_refusal(base, body, err, log):
    """Save a refused request (e.g. a provider content policy) so the run continues and the refusal is reported."""
    with open(base + ".refused.json", "w") as f:
        json.dump({"request": body, "http_status": err.code, "error": err.message,
                   "utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}, f, indent=2)
    log(f"  refused by the provider (HTTP {err.code}); saved {os.path.basename(base)}.refused.json and continuing")


def _body(model, turns, api):
    msgs = [{"role": r, "content": c} for r, c in turns]
    if api == "responses":
        return {"model": model, "input": msgs if len(msgs) > 1 else msgs[0]["content"], "store": False}
    return {"model": model, "messages": msgs}


def _post(body, key, api, log):
    path = "/responses" if api == "responses" else "/chat/completions"
    for attempt in range(4):
        try:
            return _request("POST", path, key, body)
        except urllib.error.HTTPError as e:
            msg = e.read().decode()[:500]
            if e.code in (429, 500, 502, 503) and attempt < 3:
                log(f"  retry after HTTP {e.code} ...")
                time.sleep(20 * (attempt + 1))
                continue
            raise ProviderError(e.code, msg)


def run_followup(prompt, context, model, out_dir, key=None, api="responses", log=print):
    """Second turn for every saved response_r<N>: reveal the experiment and rebuild the interpretation."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from llm_export import followup_prompt
    key = key or os.environ.get("OPENAI_API_KEY", "")
    if not key:
        raise RuntimeError('No API key: run  export OPENAI_API_KEY="sk-..."  in this terminal first')
    follow = followup_prompt(context)
    manifest = os.path.join(out_dir, "followup_manifest.tsv")
    new = not os.path.exists(manifest)
    saved = []
    with open(manifest, "a", newline="") as mf:
        w = csv.writer(mf, delimiter="\t")
        if new:
            w.writerow(["repeat", "requested_model", "returned_model", "response_id", "started_utc", "seconds",
                        "input_tokens", "output_tokens", "context", "response_file"])
        r = 1
        while os.path.exists(os.path.join(out_dir, f"response_r{r}.md")):
            base = os.path.join(out_dir, f"followup_r{r}")
            if os.path.exists(base + ".json"):
                log(f"skip (done) follow-up r{r}")
            elif os.path.exists(base + ".refused.json"):
                log(f"skip (refused earlier) follow-up r{r}")
            else:
                first = open(os.path.join(out_dir, f"response_r{r}.md")).read().strip()
                body = _body(model, [("user", prompt), ("assistant", first), ("user", follow)], api)
                started = datetime.now(timezone.utc).isoformat(timespec="seconds")
                t0 = time.time()
                try:
                    resp = _post(body, key, api, log)
                except ProviderError as err:
                    if not _is_refusal(err):
                        raise
                    _record_refusal(base, body, err, log)
                    r += 1
                    continue
                secs = round(time.time() - t0, 1)
                with open(base + ".json", "w") as f:
                    json.dump({"request": body, "response": resp, "started_utc": started, "seconds": secs}, f, indent=2)
                with open(base + ".md", "w") as f:
                    f.write(text_of(resp) + "\n")
                u = resp.get("usage", {})
                w.writerow([r, model, resp.get("model", ""), resp.get("id", ""), started, secs,
                            u.get("input_tokens", u.get("prompt_tokens", "")),
                            u.get("output_tokens", u.get("completion_tokens", "")), context, os.path.basename(base) + ".json"])
                mf.flush()
                log(f"done follow-up r{r} ({secs} s)")
            saved.append(base + ".md")
            r += 1
    if r == 1:
        raise RuntimeError(f"No response_r<N>.md in {out_dir}: run the first turn first")
    return saved


def list_models(key):
    return sorted(m["id"] for m in _request("GET", "/models", key)["data"])


def run_llm(prompt, model, out_dir, repeats=1, key=None, api="responses", log=print):
    """Send ``prompt`` ``repeats`` times (independent calls); returns the list of saved .md files.
    Finished repeats are skipped, so an interrupted run can be resumed."""
    key = key or os.environ.get("OPENAI_API_KEY", "")
    if not key:
        raise RuntimeError('No API key: run  export OPENAI_API_KEY="sk-..."  in this terminal first')
    os.makedirs(out_dir, exist_ok=True)
    manifest = os.path.join(out_dir, "manifest.tsv")
    new = not os.path.exists(manifest)
    saved = []
    with open(manifest, "a", newline="") as mf:
        w = csv.writer(mf, delimiter="\t")
        if new:
            w.writerow(["repeat", "requested_model", "returned_model", "response_id", "api", "started_utc",
                        "seconds", "input_tokens", "output_tokens", "response_file"])
        for r in range(1, repeats + 1):
            base = os.path.join(out_dir, f"response_r{r}")
            if os.path.exists(base + ".json"):
                log(f"skip (done) r{r}")
                saved.append(base + ".md")
                continue
            if os.path.exists(base + ".refused.json"):
                log(f"skip (refused earlier) r{r}")
                continue
            body = _body(model, [("user", prompt)], api)
            started = datetime.now(timezone.utc).isoformat(timespec="seconds")
            t0 = time.time()
            try:
                resp = _post(body, key, api, log)
            except ProviderError as err:
                if not _is_refusal(err):
                    raise
                _record_refusal(base, body, err, log)
                continue
            secs = round(time.time() - t0, 1)
            with open(base + ".json", "w") as f:
                json.dump({"request": body, "response": resp, "started_utc": started, "seconds": secs}, f, indent=2)
            with open(base + ".md", "w") as f:
                f.write(text_of(resp) + "\n")
            u = resp.get("usage", {})
            w.writerow([r, model, resp.get("model", ""), resp.get("id", ""), api, started, secs,
                        u.get("input_tokens", u.get("prompt_tokens", "")),
                        u.get("output_tokens", u.get("completion_tokens", "")), os.path.basename(base) + ".json"])
            mf.flush()
            log(f"done r{r} ({secs} s, model {resp.get('model', '?')})")
            saved.append(base + ".md")
    return saved


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target", nargs="?", help="iGEA run folder (containing llm_prompt.txt) or a prompt file")
    ap.add_argument("--model")
    ap.add_argument("--repeats", type=int, default=1)
    ap.add_argument("--api", choices=["responses", "chat"], default="responses")
    ap.add_argument("--list-models", action="store_true")
    ap.add_argument("--context", help="How the gene list was generated; sends the second turn for each saved response")
    a = ap.parse_args()
    key = os.environ.get("OPENAI_API_KEY", "")
    if not key:
        sys.exit('Set your key first:  export OPENAI_API_KEY="sk-..."')
    try:
        if a.list_models:
            print("\n".join(list_models(key)))
            return
        if not a.target or not a.model:
            sys.exit("Usage: python3 llm_client.py <run folder or prompt file> --model <model-id>")
        prompt_file = os.path.join(a.target, "llm_prompt.txt") if os.path.isdir(a.target) else a.target
        out_dir = os.path.join(os.path.dirname(os.path.abspath(prompt_file)), "llm")
        prompt = open(prompt_file).read()
        if a.context:
            run_followup(prompt, a.context, a.model, out_dir, key, a.api)
        else:
            run_llm(prompt, a.model, out_dir, a.repeats, key, a.api)
        print("Responses are in", out_dir)
    except urllib.error.URLError as e:
        hint = ""
        if "CERTIFICATE_VERIFY_FAILED" in str(e):
            hint = ("\nPython cannot verify SSL certificates: run /usr/bin/python3 instead, or the "
                    "'Install Certificates.command' in /Applications/Python 3.x.")
        sys.exit(f"Network error: {e}{hint}")
    except RuntimeError as e:
        sys.exit(str(e))


if __name__ == "__main__":
    main()

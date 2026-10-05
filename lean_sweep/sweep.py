"""Tactic sweep: try every handwritten script plus a generic battery on each target.
One Lean file per (target, group); each variant is its own theorem so one failure
does not block the others. Errors are mapped back to variants by line range.
Output: results/sweep.json and results/SWEEP.md. A closed variant is a kernel-checked
proof at the pinned Mathlib (v4.33.1, same as GPPVerify), modulo `import Mathlib`,
which GPPVerify forbids: trim imports when porting."""
import json, os, re, subprocess, sys, time
from targets import TARGETS, BATTERY, SEARCH

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.join(HERE, "gen"); OUT = os.path.join(HERE, "results")
os.makedirs(GEN, exist_ok=True); os.makedirs(OUT, exist_ok=True)
only = set(sys.argv[1:])


def indent(script):
    return "\n".join("  " + ln for ln in script.split("\n"))


def build(t, group, scripts):
    lines = ["import Mathlib", "set_option maxHeartbeats 400000", "set_option linter.all false", ""]
    spans = []
    for i, sc in enumerate(scripts):
        start = len(lines) + 1
        lines.append(f"theorem {t['name']}_{group}{i} {t['stmt']} := by")
        lines.extend(indent(sc).split("\n"))
        lines.append("")
        spans.append((start, len(lines), sc))
    path = os.path.join(GEN, f"{t['name']}_{group}.lean")
    open(path, "w").write("\n".join(lines) + "\n")
    return path, spans


def run(path, timeout):
    t0 = time.time()
    try:
        p = subprocess.run(["lake", "env", "lean", path], cwd=HERE, capture_output=True, text=True, timeout=timeout)
        return p.stdout + p.stderr, time.time() - t0, False
    except subprocess.TimeoutExpired as e:
        return (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), time.time() - t0, True


def parse(out, spans, timed_out):
    res = []
    msgs = [(int(m.group(1)), m.group(2), m.group(3)) for m in
            re.finditer(r":(\d+):\d+: (error|info|warning): ?(.*)", out)]
    for (a, b, sc) in spans:
        errs = [m for (ln, kind, m) in msgs if a <= ln <= b and kind == "error"]
        infos = [m for (ln, kind, m) in msgs if a <= ln <= b and kind == "info"]
        ok = (not errs) and not timed_out
        res.append(dict(script=sc, ok=ok, errors=errs[:3], info=infos[:3]))
    return res


def main():
    report = {"generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "targets": []}
    for t in TARGETS:
        if only and t["name"] not in only:
            continue
        entry = {"name": t["name"], "map": t["map"], "stmt": t["stmt"], "groups": {}}
        for group, scripts, to in [("h", t.get("tries", []), 900), ("b", BATTERY, 900), ("s", SEARCH, 1200)]:
            if not scripts:
                continue
            path, spans = build(t, group, scripts)
            out, dt, to_hit = run(path, to)
            entry["groups"][group] = dict(seconds=round(dt, 1), timed_out=to_hit, variants=parse(out, spans, to_hit),
                                          raw_tail=out[-1500:])
            print(t["name"], group, round(dt), "s", [v["ok"] for v in entry["groups"][group]["variants"]], flush=True)
        entry["closed_by"] = [v["script"] for g in entry["groups"].values() for v in g["variants"] if v["ok"]]
        report["targets"].append(entry)
    json.dump(report, open(os.path.join(OUT, "sweep.json"), "w"), indent=1, ensure_ascii=False)
    md = ["# Tactic sweep results", "", f"Generated {report['generated']} (Mathlib v4.33.1, `import Mathlib`).", "",
          "| target | map item | closed? | first closing script |", "|---|---|---|---|"]
    for e in report["targets"]:
        c = e["closed_by"]
        md.append(f"| {e['name']} | {e['map']} | {'yes' if c else 'no'} | `{c[0].splitlines()[0] if c else ''}`{' ...' if c and len(c[0].splitlines()) > 1 else ''} |")
    open(os.path.join(OUT, "SWEEP.md"), "w").write("\n".join(md) + "\n")


if __name__ == "__main__":
    main()

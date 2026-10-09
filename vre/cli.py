import argparse, csv, json, sys, os
from .scoring import score_row

ORDER = {"Very High": 0, "High": 1, "Medium": 2, "Low": 3, "Very Low": 4}


def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def infer_source(path):
    n = os.path.basename(path).lower()
    if "unmatched" in n: return "unmatched_hpo"
    if "secondary" in n: return "secondary_findings"
    return "matched_hpo"


def run(paths, patient_hpo=None):
    out = []
    for p in paths:
        src = infer_source(p)
        out += [score_row(r, src, patient_hpo) for r in load(p)]
    return sorted(out, key=lambda x: (ORDER[x.tier], -x.score))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="vre", description="Rank germline variants Very High -> Very Low with reasons.")
    ap.add_argument("tsv", nargs="+", help="Annotated TSVs (filename containing 'unmatched' or 'secondary' selects the source type)")
    ap.add_argument("--hpo", help="Comma-separated patient HPO IDs (optional)")
    ap.add_argument("--json", action="store_true", help="JSON output")
    a = ap.parse_args(argv)
    hpo = [h.strip() for h in a.hpo.split(",")] if a.hpo else None
    res = run(a.tsv, hpo)
    if a.json:
        json.dump([r.__dict__ for r in res], sys.stdout, indent=2); print()
    else:
        for r in res:
            print(f"[{r.tier:<9}] {r.score:>4}  {r.gene:<8} {r.variant}  ({r.source})")
            for why in r.reasons: print(f"             - {why}")
    return 0

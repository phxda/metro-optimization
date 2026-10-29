"""Degeneracy / constraint audit over all finished runs. Prints a summary; exits nonzero on any violation."""
import sys, os, json, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import cities, model

def main():
    bad, n, warn = [], 0, []
    for p in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "runs", "*.json"))):
        r = json.load(open(p)); n += 1
        c = cities.get_city(r["city"]); L = [[tuple(x) for x in l] for l in r["lines"]]
        if model.net_cost(c, L) > r["budget_units"] + 1e-6: bad.append((r["id"], "over budget"))
        for l in L:
            if any(c.water[a] for a in l): bad.append((r["id"], "station in water"))
            if len(set(l)) != len(l): bad.append((r["id"], "duplicate station"))
            if not all(c.min_hop - 1e-9 <= c.seg(a, b)[0] <= c.max_hop + 1e-9 for a, b in zip(l[:-1], l[1:])): bad.append((r["id"], "hop out of range"))
        if len(L) * model.V_MIN > r["fleet"]: bad.append((r["id"], "fleet infeasible"))
        m = r["metrics"]
        if not m: warn.append((r["id"], "empty network")); continue
        if m["longest_share"] > 0.8 and m["n_lines"] > 1: warn.append((r["id"], "one line dominates"))
        if m["n_lines"] == 1: warn.append((r["id"], "single line"))
        if m["fitness"] <= 0: warn.append((r["id"], "non-positive fitness"))
    print(f"{n} runs audited; {len(bad)} violations; {len(warn)} warnings")
    for b in bad[:20]: print("VIOLATION", b)
    from collections import Counter
    print(Counter(w[1] for w in warn))
    sys.exit(1 if bad else 0)
main()

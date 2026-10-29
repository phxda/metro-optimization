"""Evaluate the stylized real/proposed comparison networks with the same fitness function."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
import cities, model
out = {}
for name, sets in {
    "baltimore": {"Metro SubwayLink": ["Metro SubwayLink"], "Light RailLink": ["Light RailLink"],
                  "Red Line (planned)": ["Red Line (planned)"],
                  "Metro + Light Rail": ["Metro SubwayLink", "Light RailLink"],
                  "Metro + Light Rail + Red Line": ["Metro SubwayLink", "Light RailLink", "Red Line (planned)"]},
    "triangle": {"Proposed commuter rail": ["Proposed commuter rail (Durham-Garner-Clayton)"],
                 "S-Line": ["S-Line (Raleigh-Wake Forest)"],
                 "Proposed + S-Line": ["Proposed commuter rail (Durham-Garner-Clayton)", "S-Line (Raleigh-Wake Forest)"]}}.items():
    c = cities.get_city(name)
    rn = c.real_networks()
    out[name] = {}
    for label, ks in sets.items():
        lines = [l for k in ks for l in rn[k]]
        for pk, pen in (("lowpen", 5.0), ("highpen", 20.0)):
            f, m = model.evaluate(c, lines, pen, model.fleet_for(c.mid_budget), detail=True)
            out[name][f"{label}|{pk}"] = dict(label=label, pen=pk, lines=lines, metrics=m)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "real_eval.json"), "w"))
for n, d in out.items():
    for k, v in d.items():
        print(n, k, round(v["metrics"]["fitness"]), round(v["metrics"]["cost"]), round(v["metrics"]["mode_share"], 3))

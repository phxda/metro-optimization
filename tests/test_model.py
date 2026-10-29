"""Hand-computed checks on tiny networks. Run: python3 tests/test_model.py"""
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import cities, model
from cities import City

def tiny(R=1, C=11, water_cols=(), cell=1.0, **kw):
    water = np.zeros((R, C), bool)
    for c in water_cols:
        water[:, c] = True
    H = np.ones((R, C)); J = np.ones((R, C))
    return City("tiny", "tiny", R, C, cell, water, H, J, block=1, transit_kmh=60., dwell=1.0, wait_w=2.0, access_w=2.0,
                total_workers=1000, **kw)

def approx(a, b, tol=1e-4):
    assert abs(a - b) < tol, (a, b)

def test_cost():
    c = tiny(water_cols=(5,))
    # line of 3 stations at cols 0,3,8 (r=0): hops 3 km and 5 km; second hop crosses water col 5
    ln = [(0, 0), (0, 3), (0, 8)]
    k1, w1 = c.seg((0, 0), (0, 3)); k2, w2 = c.seg((0, 3), (0, 8))
    approx(k1, 3.0); approx(k2, 5.0); approx(w1, 0.0)
    assert 0.1 < w2 < 2.0     # water cell sampled along hop: about 1 of 5 km
    exp = 2.5 * 8 + 1.5 * 3 + 4.0 * w2
    approx(model.net_cost(c, [ln]), exp)
    # second line sharing station (0,3): stations charged once, track charged per hop
    ln2 = [(0, 3), (0, 1)]
    approx(model.net_cost(c, [ln, ln2]), exp + 2.5 * 2 + 1.5 * 1)

def test_travel_time():
    c = tiny()
    ln = [(0, 0), (0, 3), (0, 8)]
    pen = 10.0
    # fleet gives headway: cycle = 2*(3/60*60+1 + 5/60*60+1) + 2*4 = 2*(4+6)+8 = 28 min; fleet=14 -> veh 14 -> hw 2 (floor 2)
    hw = model.headways(c, [ln], 14)
    approx(hw[0], 2.0)
    idx, hw, T, pred, nid = model.station_times(c, [ln], pen, 14, want_pred=True)
    i0, i3, i8 = idx[(0, 0)], idx[(0, 3)], idx[(0, 8)]
    # 0->3 : wait*2 (=2*2/2=2) + ride 4 = 6 ; 0->8: 2 + 4 + 6 = 12
    approx(T[i0, i3], 6.0); approx(T[i0, i8], 12.0); approx(T[i3, i8], 2 + 6)
    approx(T[i8, i0], 12.0)
    # headway with fleet 4: cycle 28/4 = 7 min -> wait weight 2 * 3.5 = 7
    hw4 = model.headways(c, [ln], 4)
    approx(hw4[0], 7.0)
    _, _, T4, _, _ = model.station_times(c, [ln], pen, 4)
    approx(T4[i0, i3], 7 + 4)

def test_transfers():
    c = tiny(R=5, C=5)
    # A: horizontal row 2, cols 0,2,4 ; B: vertical col 2, rows 0,2,4 ; share (2,2)
    A = [(2, 0), (2, 2), (2, 4)]; B = [(0, 2), (2, 2), (4, 2)]
    pen = 10.0
    idx, hw, T, pred, nid = model.station_times(c, [A, B], pen, 20, want_pred=True)
    P = len(idx); assert P == 5
    X = model.transfer_counts(P, pred, nid)
    a0, b0, hub, a4, b4 = idx[(2, 0)], idx[(0, 2)], idx[(2, 2)], idx[(2, 4)], idx[(4, 2)]
    assert X[a0, a4] == 0 and X[a0, hub] == 0 and X[b0, b4] == 0
    assert X[a0, b0] == 1 and X[a4, b4] == 1 and X[b0, a4] == 1
    # times: speeds 60km/h dwell 1: hop of 2 km = 2+1 = 3 min. cycle = 2*(3+3)+8 = 20; fleet 20 split into 2 lines
    # -> 10 veh each -> hw 2 -> wait term 2. a0->hub = 2 + 3 = 5 ; a0->b0 = 2+3 + pen + (2+3) = 5+10+5 = 20
    approx(T[a0, hub], 5.0); approx(T[a0, b0], 20.0); approx(T[a0, a4], 2 + 6)
    # same hub: single-line vs transfer counts via detail path
    lines = [A, B]
    # a line count above fleet/V_MIN is infeasible
    assert model.headways(c, [A, B], 3) is None

def test_fitness_hand():
    # 2 zones far apart (cols 0 and 9), single line 0->9 with 1 hop. Check fitness against closed form.
    c = tiny(R=1, C=10, drive_kmh=30., circuity=1.0, park_base=0., park_max=0., cong=0., access_radius_km=0.5)
    ln = [(0, 0), (0, 9)]
    fleet = 20
    cyc = 2 * (9 + 1) + 8   # 28
    hw = max(2.0, cyc / fleet)   # 1.4 -> 2.0
    ride = 9 + 1
    tr = 2 * hw / 2 + ride + 0   # wait weighted 2 * (hw/2)
    # access: zone centroid at the cell center; offset sqrt(d^2+0.3^2) with d=0 -> 0.3 km
    acc = 2.0 * 0.3 / 4.8 * 60
    ct = tr + 2 * acc
    f = model.evaluate(c, [ln], 10., fleet)
    # drive time between cell 0 and cell 9 = 9 km /30 *60 = 18 min
    dt = 18.0
    d = dt - ct
    # only zones 0 and 9 can be served. Flow comes from the gravity model
    z0 = int(np.argmin(np.abs(c.zx - 0.5))); z9 = int(np.argmin(np.abs(c.zx - 9.5)))
    exp = 0.0
    for a, b in ((z0, z9), (z9, z0)):
        p = 1 / (1 + math.exp(-(cities.BETA * d + c.asc_z[a])))
        exp += c.flow[a, b] * p * d
    approx(f, exp / 60, 1e-4)

def test_budget_and_validity():
    c = cities.get_city("baltimore")
    nets = c.real_networks()
    for name, net in nets.items():
        for ln in net:
            assert all(not c.water[r, cc] for r, cc in ln), name
            assert len(set(ln)) == len(ln)

def test_symmetry_and_monotone():
    c = cities.get_city("mono")
    ln = [(14, c) for c in range(2, 27, 2)]
    f1 = model.evaluate(c, [ln], 5., 40)
    f2 = model.evaluate(c, [ln], 20., 40)
    assert f1 > 0 and f2 > 0
    assert f1 >= f2 - 1e-9          # penalty cannot help a single line (no transfers)... equal
    approx(f1, f2, 1e-4)

if __name__ == "__main__":
    n = 0
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); n += 1; print("ok", k)
    print(n, "tests passed")

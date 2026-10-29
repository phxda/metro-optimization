"""Synthetic and stylized-real city definitions. All real-city geography is a rough
hand-built approximation from memory (no data downloaded)."""
import math
import numpy as np

# ---- global cost / operating parameters (model units: 1 unit ~ $100M) ----
COST_KM = 2.5      # per track km
COST_STATION = 1.5  # per station
COST_WATER_KM = 4.0  # extra per km of water crossed
FLEET_PER_UNIT = 0.2  # vehicles per budget unit (fleet is fixed given the budget)
LAYOVER_MIN = 4.0
V_MIN = 2.0          # min vehicles per line
HEADWAY_MIN, HEADWAY_MAX = 2.0, 30.0
BETA = 0.15          # logit slope per generalized minute
ASC = -0.5           # transit alternative-specific constant
ACCESS_K = 4         # nearest stations considered per zone


def blob(X, Y, x0, y0, s, w=1.0):
    return w * np.exp(-((X - x0) ** 2 + (Y - y0) ** 2) / (2 * s * s))


class City:
    def __init__(self, name, label, R, C, cell_km, water, H, J, carless=None, side=None,
                 block=2, drive_kmh=28., circuity=1.3, park_base=6., park_max=16., cong=1.2, side_pen=0.,
                 transit_kmh=40., dwell=0.5, walk_kmh=4.8, access_radius_km=1.2,
                 access_w=1.2, wait_w=1.5, gamma_km=8., min_hop_km=0.9, max_hop_km=6.,
                 max_lines=10, asc=ASC, beta=BETA, fleet_per_unit=FLEET_PER_UNIT, mid_budget=220., real=None, total_workers=200000, desc=""):
        self.name, self.label, self.R, self.C, self.cell = name, label, R, C, cell_km
        self.water = water.astype(bool)
        self.H = np.where(self.water, 0, H)
        self.J = np.where(self.water, 0, J)
        self.carless = np.zeros((R, C)) if carless is None else carless
        self.side = np.zeros((R, C), int) if side is None else side
        self.block = block
        self.drive_kmh, self.circ, self.side_pen = drive_kmh, circuity, side_pen
        self.park_base, self.park_max, self.cong = park_base, park_max, cong
        self.transit_kmh, self.dwell, self.walk_kmh = transit_kmh, dwell, walk_kmh
        self.access_radius, self.access_w, self.wait_w = access_radius_km, access_w, wait_w
        self.gamma, self.min_hop, self.max_hop, self.max_lines = gamma_km, min_hop_km, max_hop_km, max_lines
        self.mid_budget = mid_budget
        self.asc = asc
        self.beta, self.fleet_per_unit = beta, fleet_per_unit
        self.real_xy = real or {}
        self.desc = desc
        s = total_workers / self.H.sum()
        self.H = self.H * s
        self.J = self.J * (total_workers / self.J.sum())
        self._seg = {}
        self._prep_zones()
        # demand weight for mutation sampling
        w = self.H + self.J
        k = np.pad(w, 1)
        self.demand_w = sum(k[i:i + R, j:j + C] for i in range(3) for j in range(3))
        self.demand_w = np.where(self.water, 0, self.demand_w)
        self.land = ~self.water

    def xy(self, r, c):
        return (c + 0.5) * self.cell, (r + 0.5) * self.cell

    def cell_of(self, x, y):
        return int(min(self.R - 1, max(0, y / self.cell))), int(min(self.C - 1, max(0, x / self.cell)))

    def _prep_zones(self):
        b, R, C = self.block, self.R, self.C
        Rz, Cz = -(-R // b), -(-C // b)
        pr, pc = Rz * b - R, Cz * b - C

        def pad(a):
            return np.pad(a, ((0, pr), (0, pc)))
        H, J, Q = pad(self.H), pad(self.J), pad(self.carless * self.H)
        rr, cc = np.meshgrid(np.arange(Rz * b), np.arange(Cz * b), indexing="ij")
        X, Y = (cc + .5) * self.cell, (rr + .5) * self.cell

        def blk(a):
            return a.reshape(Rz, b, Cz, b).sum((1, 3))
        Hz, Jz, Qz = blk(H), blk(J), blk(Q)
        W = H + J + 1e-9
        wz = blk(W)
        zx, zy = blk(W * X) / wz, blk(W * Y) / wz
        side = pad(self.side)
        sz = side.reshape(Rz, b, Cz, b)[:, 0, :, 0]
        keep = ((Hz + Jz) > 1e-6 * (Hz.sum())).ravel()
        self.zx, self.zy = zx.ravel()[keep], zy.ravel()[keep]
        self.zH, self.zJ = Hz.ravel()[keep], Jz.ravel()[keep]
        self.zq = (Qz.ravel() / np.maximum(Hz.ravel(), 1e-9))[keep]
        self.zside = sz.ravel()[keep]
        self.Z = len(self.zx)
        d = np.hypot(self.zx[:, None] - self.zx[None], self.zy[:, None] - self.zy[None])
        bk = b * self.cell
        np.fill_diagonal(d, 0.5 * bk)
        f = np.exp(-d / self.gamma)
        f *= self.zJ[None, :]
        np.fill_diagonal(f, 0)  # intrazone trips are walked, excluded
        flow = self.zH[:, None] * f / f.sum(1, keepdims=True).clip(1e-12)
        self.flow = 2.0 * flow  # one-way trips per day (2 per worker)
        area = bk * bk
        g = np.sqrt((Hz.ravel()[keep] + Jz.ravel()[keep]) / area)
        g = np.clip(g / g.max(), 0, 1)           # congestion / parking-pressure index in [0,1]
        self.zg = g
        self.drive_t = (self.park_base + self.park_max * g[None, :]
                        + 60 * self.circ * d / self.drive_kmh * (1 + self.cong * (g[:, None] + g[None, :]) / 2)
                        + self.side_pen * (self.zside[:, None] != self.zside[None, :]))
        self.asc_z = self.asc + 3.0 * self.zq
        self.total_trips = self.flow.sum()

    # --- geometry helpers ---
    def seg(self, a, b):
        """(km, water_km) for a straight hop between cells a, b."""
        k = (a, b) if a <= b else (b, a)
        v = self._seg.get(k)
        if v is None:
            (r1, c1), (r2, c2) = k
            n = max(2, int(4 * max(abs(r1 - r2), abs(c1 - c2))))
            t = (np.arange(n) + .5) / n
            rr = np.clip(np.round(r1 + (r2 - r1) * t).astype(int), 0, self.R - 1)
            cc = np.clip(np.round(c1 + (c2 - c1) * t).astype(int), 0, self.C - 1)
            km = math.hypot(r1 - r2, c1 - c2) * self.cell
            v = (km, km * float(self.water[rr, cc].mean()))
            self._seg[k] = v
        return v

    def snap(self, x, y):
        r0, c0 = self.cell_of(x, y)
        best, bd = None, 1e9
        for dr in range(-4, 5):
            for dc in range(-4, 5):
                r, c = r0 + dr, c0 + dc
                if 0 <= r < self.R and 0 <= c < self.C and not self.water[r, c]:
                    d = dr * dr + dc * dc
                    if d < bd:
                        best, bd = (r, c), d
        return best

    def real_networks(self):
        out = {}
        for name, lines in self.real_xy.items():
            net = []
            for ln in lines:
                st = []
                for x, y in ln:
                    p = self.snap(x, y)
                    if p is not None and (not st or st[-1] != p) and p not in st:
                        st.append(p)
                if len(st) >= 2:
                    net.append(st)
            out[name] = net
        return out


def _grid(R, C, cell):
    rr, cc = np.meshgrid(np.arange(R), np.arange(C), indexing="ij")
    return (cc + .5) * cell, (rr + .5) * cell


def make_mono():
    R = C = 28
    X, Y = _grid(R, C, 1.)
    d = np.hypot(X - 14, Y - 14)
    mask = d < 13.5
    H = (blob(X, Y, 14, 14, 6.5) + .05) * mask
    J = (blob(X, Y, 14, 14, 2.5) + .12 * blob(X, Y, 14, 14, 7)) * mask
    return City("mono", "Monocentric (round)", R, C, 1., ~mask, H, J, mid_budget=220.,
                desc="Round city, single CBD.")


def make_coastal():
    R, C = 18, 40
    X, Y = _grid(R, C, 1.)
    coast = 12.5 + 1.2 * np.sin(X / 6.)
    water = Y > coast
    inland = coast - Y
    H = np.exp(-np.clip(inland, 0, None) / 3.5) * (1 + .25 * np.sin(X / 4.)) * (X > 1) * (X < 39)
    J = (blob(X, Y, 20, 9.5, 3, 1.0) + blob(X, Y, 8, 10, 2.5, .45)
         + blob(X, Y, 32, 10, 2.5, .45)) * 1.0
    return City("coastal", "Linear coastal", R, C, 1., water, H, J, mid_budget=220.,
                desc="Elongated strip along a coast; main center mid-strip plus two secondary centers.")


def make_river():
    R = C = 28
    X, Y = _grid(R, C, 1.)
    rc = 13.5 + 1.5 * np.sin(Y / 5.)
    water = np.abs(X - rc) < 0.8
    side = (X > rc).astype(int)
    H = blob(X, Y, 14, 14, 7.5) + .05
    J = blob(X, Y, 18, 14, 2.5, 1.0) + blob(X, Y, 7, 12, 2.2, .45) + .1 * blob(X, Y, 14, 14, 8)
    return City("river", "River-split", R, C, 1., water, H, J, side=side, side_pen=8.,
                mid_budget=220., desc="Vertical river with costly crossings; big CBD east bank, minor center west bank.")


def make_poly():
    R = C = 28
    X, Y = _grid(R, C, 1.)
    H = blob(X, Y, 14, 14, 9) + .15
    J = (blob(X, Y, 7, 7, 2.2, 1.0) + blob(X, Y, 21, 7, 2.2, .8) + blob(X, Y, 8, 21, 2.2, .7)
         + blob(X, Y, 20, 20, 2.2, 1.0) + .05)
    return City("poly", "Polycentric (4 centers)", R, C, 1., np.zeros((R, C), bool), H, J, mid_budget=220.,
                desc="Four job centers on a broad homes field.")


def _patapsco_water(X, Y):
    pts = [(10, 29.2), (15, 28.4), (19, 28.0), (22, 28.6), (26, 29.6), (30, 31.6), (34, 33.5), (37, 35)]
    xs, ys = zip(*pts)
    yc = np.interp(X, xs, ys)
    hw = np.interp(X, [10, 22, 34], [1.0, 1.1, 2.4])
    w = (np.abs(Y - yc) < hw) & (X > 12.5)
    return w | (X > 33.5)


def make_baltimore():
    R, C = 38, 36
    X, Y = _grid(R, C, 1.)
    water = _patapsco_water(X, Y)
    Hb = [(10, 23, 3.5, 1.0), (7, 25, 2.5, .6), (22, 22.5, 3, .9), (17.5, 17.5, 3.5, 1.0), (17.5, 25, 2.0, .6),
          (27, 24.5, 3, .7), (11, 15, 4, .8), (20, 9.5, 3.5, .7), (25, 14, 4, .6), (15, 33.5, 3.5, .45),
          (8, 29.5, 3, .45), (6, 13, 3, .35), (20, 3.5, 2.5, .15), (12, 20, 3, .5)]
    Jb = [(17.3, 24.3, 1.8, 1.0), (20.6, 23.6, 1.2, .38), (17.5, 18.5, 1.2, .22), (20, 10, 1.8, .25),
          (10, 35, 2.0, .12), (5.5, 21.5, 1.5, .12), (28, 23, 1.3, .12), (20.5, 3.5, 1.8, .1),
          (13, 19.5, 1.2, .12), (14, 15, 3, .1)]
    H = sum(blob(X, Y, *b) for b in Hb) + .02
    J = sum(blob(X, Y, *b) for b in Jb) + .02 * blob(X, Y, 17, 18, 11)
    q = (blob(X, Y, 10, 23, 3.5, .30) + blob(X, Y, 22, 22.5, 3.0, .30) + blob(X, Y, 17.5, 25, 2.5, .25)
         + blob(X, Y, 17.5, 17.5, 3.5, .14) + blob(X, Y, 12, 20, 3, .2))
    q = np.clip(q, 0, .40) + .03
    real = {
        "Metro SubwayLink": [[(5.5, 12.8), (6.6, 14.3), (8, 16), (9.6, 17.7), (10.8, 18.9), (12, 19.5), (12.8, 20.8),
                              (13.8, 21.8), (14.8, 22.9), (15.8, 23.5), (16.5, 24.3), (17.3, 24.5), (18.4, 24.6),
                              (20.6, 24.1)]],
        "Light RailLink": [[(20.5, 2.5), (20.3, 6), (20, 8), (19, 10), (17, 12.5), (14.8, 14.8), (15.4, 16.4),
                            (16.2, 17.8), (16.4, 21.2), (16.7, 22.4), (16.5, 24.3), (16.8, 25.5), (16.0, 30.2),
                            (15.2, 31.8), (13.0, 33.0), (11.5, 34.2), (9.5, 36)]],
        "Red Line (planned)": [[(4.5, 22.5), (6.6, 22.6), (9, 23), (11, 23.2), (12.8, 23.5), (14.6, 24), (16.8, 24.7),
                                (18.8, 25.3), (20.8, 25.9), (22.6, 25.6), (23.8, 25.2), (25.2, 24.7), (27.8, 23.4)]],
    }
    return City("baltimore", "Baltimore (stylized)", R, C, 1., water, H, J, carless=q, drive_kmh=26.,
                side_pen=0., max_hop_km=7., mid_budget=None, real=real, total_workers=420000,
                desc="Stylized approximation of Baltimore: downtown/Inner Harbor, JHH, Homewood, water to S and E.")


def make_triangle():
    R = C = 36
    cell = 1.7
    X, Y = _grid(R, C, cell)
    Hb = [(38, 29, 9, 1.0), (14, 12, 6, .6), (6, 22, 4, .3), (28, 34, 6, .6), (22, 43, 5, .25), (43, 6, 4, .25),
          (44, 40, 6, .35), (24, 24, 4, .15), (52, 26, 5, .3), (33, 17, 6, .5)]
    Jb = [(37.5, 29, 3.2, 1.0), (24.5, 22, 3.5, .9), (14, 12, 2.8, .8), (6, 22, 2.5, .5), (29, 33, 3, .35),
          (43, 6, 2, .1), (36, 20, 4, .5), (44, 40, 2.5, .08)]
    H = sum(blob(X, Y, *b) for b in Hb) + .03
    J = sum(blob(X, Y, *b) for b in Jb) + .03 * blob(X, Y, 30, 25, 20)
    q = np.full((R, C), .04) + blob(X, Y, 37.5, 29, 3., .06)
    real = {
        "Proposed commuter rail (Durham-Garner-Clayton)": [[(14, 12), (19, 17), (24, 22), (26.5, 28), (29, 33), (34, 31),
                                                              (38, 29), (42, 37), (48, 45)]],
        "S-Line (Raleigh-Wake Forest)": [[(38, 29), (39.5, 22), (41, 15), (42.5, 9), (43.5, 5.5)]],
    }
    return City("triangle", "Triangle (stylized)", R, C, cell, np.zeros((R, C), bool), H, J, carless=q, block=2,
                drive_kmh=50., circuity=1.25, park_base=5., park_max=6., cong=0.5, asc=-1.5, transit_kmh=60., dwell=1.0, access_radius_km=6.,
                walk_kmh=25., access_w=1.2, gamma_km=18., min_hop_km=1.6, max_hop_km=14., mid_budget=None,
                real=real, total_workers=700000,
                desc="Stylized Raleigh-Durham-Chapel Hill region; park-and-ride/feeder access (15 km/h, 5 km).")


CITY_MAKERS = {"mono": make_mono, "coastal": make_coastal, "river": make_river, "poly": make_poly,
               "baltimore": make_baltimore, "triangle": make_triangle}
CITY_ORDER = ["mono", "coastal", "river", "poly", "baltimore", "triangle"]
_cache = {}


BUDGET_NETS = {"baltimore": ["Metro SubwayLink", "Light RailLink", "Red Line (planned)"],
               "triangle": ["Proposed commuter rail (Durham-Garner-Clayton)"]}


def get_city(name):
    if name not in _cache:
        c = CITY_MAKERS[name]()
        if c.mid_budget is None:       # real cities: mid budget = model cost of the real/proposed network
            import model
            rn = c.real_networks()
            c.mid_budget = float(math.ceil(model.net_cost(c, [l for k in BUDGET_NETS[name] for l in rn[k]])))
        _cache[name] = c
    return _cache[name]


VARIANTS = {  # one-at-a-time perturbations used by the sensitivity sweep (v2)
    "base": {},
    "drive_slow": {"drive_kmh_mult": 0.9}, "drive_fast": {"drive_kmh_mult": 1.1},
    "asc_lo": {"asc_add": -0.5}, "asc_hi": {"asc_add": 0.5},
    "beta_lo": {"beta": 0.10}, "beta_hi": {"beta": 0.22},
    "fleet_lo": {"fleet_per_unit": 0.1}, "fleet_hi": {"fleet_per_unit": 0.4},
}
_vcache = {}


def get_variant_city(name, variant="base"):
    """City with perturbed behavioral parameters; budgets (mid_budget) stay those of the base city."""
    if variant == "base":
        return get_city(name)
    k = (name, variant)
    if k not in _vcache:
        import copy
        c = copy.deepcopy(get_city(name))
        v = VARIANTS[variant]
        if "drive_kmh_mult" in v:
            c.drive_kmh *= v["drive_kmh_mult"]
        if "asc_add" in v:
            c.asc += v["asc_add"]
        if "beta" in v:
            c.beta = v["beta"]
        if "fleet_per_unit" in v:
            c.fleet_per_unit = v["fleet_per_unit"]
        c._seg = get_city(name)._seg
        c._prep_zones()
        _vcache[k] = c
    return _vcache[k]

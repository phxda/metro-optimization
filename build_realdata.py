"""Builds data/{baltimore,triangle}_lodes.npz from downloaded LEHD LODES8 (2022) OD + block crosswalk and TIGER area-water.
Run with: python3 -I build_realdata.py   (raw inputs are untrusted data in data/raw/, read as data only)"""
import os, sys, io, zipfile, glob, math
ROOT = os.path.dirname(os.path.abspath(__file__))
import numpy as np, pandas as pd, shapefile
from matplotlib.path import Path

FRAMES = {  # name: (state, grid R, C, cell km, anchor lat, anchor lon, anchor x km, anchor y km, water county FIPS)
    "baltimore": ("md", 38, 36, 1.0, 39.2904, -76.6122, 17.3, 24.3, ["24510", "24005", "24003", "24025", "24027"]),
    "triangle": ("nc", 36, 36, 1.7, 35.7796, -78.6382, 40.0, 32.0, ["37183", "37063", "37135", "37037", "37101"]),
}


def project(lat, lon, f):
    _, R, C, cell, la0, lo0, ax, ay, _ = f
    x = ax + (lon - lo0) * 111.32 * math.cos(math.radians(la0))
    y = ay - (lat - la0) * 110.57
    return x, y


def build(name):
    f = FRAMES[name]
    st, R, C, cell = f[0], f[1], f[2], f[3]
    xw = pd.read_csv(os.path.join(ROOT, "data/raw", st, "xwalk.csv.gz"), usecols=["tabblk2020", "blklatdd", "blklondd"], dtype={"tabblk2020": str})
    la, lo = xw.blklatdd.values, xw.blklondd.values
    x = f[6] + (lo - f[5]) * 111.32 * math.cos(math.radians(f[4]))
    y = f[7] - (la - f[4]) * 110.57
    r = np.floor(y / cell).astype(int); c = np.floor(x / cell).astype(int)
    inside = (r >= 0) & (r < R) & (c >= 0) & (c < C)
    cellid = pd.Series(np.where(inside, r * C + c, -1), index=xw.tabblk2020.values)
    od = pd.read_csv(os.path.join(ROOT, "data/raw", st, "od.csv.gz"), usecols=["w_geocode", "h_geocode", "S000"], dtype={"w_geocode": str, "h_geocode": str})
    od["h"] = cellid.reindex(od.h_geocode.values).values
    od["w"] = cellid.reindex(od.w_geocode.values).values
    od = od.dropna(subset=["h", "w"])
    od = od[(od.h >= 0) & (od.w >= 0)]
    g = od.groupby(["h", "w"]).S000.sum().reset_index()
    src, dst, cnt = g.h.values.astype(np.int64), g.w.values.astype(np.int64), g.S000.values.astype(float)
    H = np.bincount(src, cnt, R * C).reshape(R, C)
    J = np.bincount(dst, cnt, R * C).reshape(R, C)
    # water: TIGER area water polygons >= 0.3 km2, rasterized at cell centres; cells with real activity stay land
    wat = np.zeros((R, C), bool)
    cx = (np.arange(C) + .5) * cell; cy = (np.arange(R) + .5) * cell
    X, Y = np.meshgrid(cx, cy)
    pts = np.c_[X.ravel(), Y.ravel()]
    for fips in f[8]:
        z = zipfile.ZipFile(os.path.join(ROOT, "data/raw/water", fips + ".zip"))
        tmp = os.path.join(ROOT, "data/raw/water", "x" + fips)
        z.extractall(tmp)
        sf = shapefile.Reader(glob.glob(tmp + "/*.shp")[0])
        for shp in sf.shapes():
            if len(shp.points) < 4:
                continue
            P = np.array(shp.points)
            px = f[6] + (P[:, 0] - f[5]) * 111.32 * math.cos(math.radians(f[4])); py = f[7] - (P[:, 1] - f[4]) * 110.57
            if px.max() < 0 or px.min() > C * cell or py.max() < 0 or py.min() > R * cell:
                continue
            parts = list(shp.parts) + [len(P)]
            for a, b in zip(parts[:-1], parts[1:]):
                ring = np.c_[px[a:b], py[a:b]]
                area = 0.5 * abs(np.dot(ring[:, 0], np.roll(ring[:, 1], 1)) - np.dot(ring[:, 1], np.roll(ring[:, 0], 1)))
                if area < 0.3:
                    continue
                m = Path(ring).contains_points(pts)
                wat |= m.reshape(R, C)
    active = (H + J) >= 60
    water = wat & ~active
    np.savez_compressed(os.path.join(ROOT, "data", name + "_lodes.npz"), R=R, C=C, cell=cell, H=H, J=J, water=water, src=src, dst=dst, cnt=cnt,
                        lat0=f[4], lon0=f[5], ax=f[6], ay=f[7])
    print(name, "workers in grid pairs", int(cnt.sum()), "water cells", int(water.sum()), "pairs", len(cnt),
          "H max", int(H.max()), "J max", int(J.max()))


if __name__ == "__main__":
    for n in FRAMES:
        build(n)

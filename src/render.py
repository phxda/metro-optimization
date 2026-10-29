import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAL = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]


def draw(ax, city, lines, title=None, legend=False, lw=2.2, ms=3.5):
    ext = [0, city.C * city.cell, city.R * city.cell, 0]
    dens = np.log1p((city.H + city.J) / (city.cell ** 2))
    img = np.where(city.water, np.nan, dens)
    ax.imshow(img, extent=ext, cmap="Greys", vmin=0, vmax=np.nanmax(dens) * 1.3, alpha=.75)
    wm = np.where(city.water, 1.0, np.nan)
    ax.imshow(wm, extent=ext, cmap="Blues", vmin=0, vmax=1.6, alpha=.9)
    for i, ln in enumerate(lines):
        pts = np.array([city.xy(*c) for c in ln])
        ax.plot(pts[:, 0], pts[:, 1], "-", color=PAL[i % len(PAL)], lw=lw, solid_capstyle="round", label=f"L{i+1}")
    cnt = {}
    for ln in lines:
        for c in ln:
            cnt[c] = cnt.get(c, 0) + 1
    for c, k in cnt.items():
        x, y = city.xy(*c)
        if k >= 2:
            ax.plot(x, y, "o", mfc="white", mec="k", ms=ms + 2, mew=1.2, zorder=5)
        else:
            ax.plot(x, y, "o", color="k", ms=ms * .7, zorder=5)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(0, city.C * city.cell); ax.set_ylim(city.R * city.cell, 0)
    if title:
        ax.set_title(title, fontsize=8)

#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Vergleich der Kelvin-Korrektur in der Evaporation (L500):
#   kelvin_cap : klassischer Kelvin-Faktor, unter d_min = 1 nm gedeckelt (bisheriger Stand)
#                K = exp(4*sigma*v_m / (max(d_p, d_min)*k_B*T))
#   tolman     : Fe mit Tolman-korrigierter Oberflächenspannung, Khrabry et al. (2025),
#                Aerosol Sci. Technol. 59(6), 743-764, Gl. 7, 15-17
#                sigma_i = sigma / (1 + 2*(delta_T/r_W)*i^(-1/3))
#                K = exp((E_s(i) - E_s(i-1)) / (k_B*T)),  E_s(i) = sigma_i*4*pi*r_W^2*(i^(2/3)-1)
#                delta_T/r_W aus dem MD-Fit (T auf 1600-2400 K begrenzt); FeO weiter mit Deckel
# Beide Varianten: gleicher Code, gleiche Flamme (cantera/nuclest aus RESULTS/L500_results),
# Umschaltung über YAML evaporation.kelvin_mode. Läufe in RUNS_Tolman-Vergleich_2026-11-26/.
#
# Usage: python SRC/scripts/compare_tolman.py [--window-mm 0.1]
# -----------------------------------------------------------------------------
import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--case", default="L500")
p.add_argument("--runs-dir", default="RUNS_Tolman-Vergleich_2026-11-26")
p.add_argument("--window-mm", type=float, default=0.1, help="Glättungsfenster in mm (0 = aus)")
args = p.parse_args()
CASE, WIN = args.case, args.window_mm

BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RUNS = os.path.join(BASE, args.runs_dir)
SETS = {
    "Kelvin (1-nm-Deckel)": dict(dir=os.path.join(RUNS, "kelvin_cap"), color="0.3"),
    "Kelvin + Tolman":      dict(dir=os.path.join(RUNS, "tolman"), color="tab:blue"),
}
OUT_DIR = os.path.join(BASE, "Plots", "Comparison")
RAW = dict(lw=0.5, alpha=0.2)


def res(d, f):
    return os.path.join(d, "RESULTS", f"{CASE}_results", f"{CASE}-{f}")


def smooth(z_mm, y):
    """Centred moving average over WIN mm; edge-normalised and NaN-safe."""
    if WIN <= 0 or len(y) < 3:
        return y
    n = max(1, int(round(WIN / np.median(np.diff(z_mm)))))
    if n < 2:
        return y
    k = np.ones(n)
    valid = np.isfinite(y)
    num = np.convolve(np.where(valid, y, 0.0), k, mode="same")
    den = np.convolve(valid.astype(float), k, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def read_runs(d):
    """Per run: end height, final n_Fe, max n_Fe, z of every evaporation event."""
    out = []
    for f in sorted(glob.glob(res(d, "SPMC-*.csv")), key=lambda s: int(s.rsplit("-", 1)[1][:-4])):
        z_end, nfe_end, nfe_max, z_ev = np.nan, 0, 0, []
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = [x.strip() for x in line.split(",")]
                if len(c) < 14:
                    continue
                z_end, nfe_end = float(c[1]), int(c[10])
                nfe_max = max(nfe_max, nfe_end)
                if c[13] == "ev":
                    z_ev.append(float(c[1]))
        out.append(dict(z_end=z_end * 1e3, nfe_end=nfe_end, nfe_max=nfe_max, z_ev=np.array(z_ev) * 1e3))
    return out


data = {}
for name, cfg in SETS.items():
    d = cfg["dir"]
    data[name] = dict(avg=np.genfromtxt(res(d, "average.csv"), delimiter=","),
                      cnt=np.genfromtxt(res(d, "count.csv"), delimiter=","),
                      runs=read_runs(d))
    cfg["n"] = len(data[name]["runs"])
cant = np.genfromtxt(res(SETS["Kelvin + Tolman"]["dir"], "cantera.csv"), delimiter=",", names=True)

fig, axs = plt.subplots(2, 3, figsize=(16, 9))
(ax_T, ax_D, ax_ev), (ax_c, ax_r, ax_run) = axs
ax_n = ax_run.twinx()
T_s_all = []
zmax = 0.0

for j, (name, cfg) in enumerate(SETS.items()):
    avg, cnt, runs = data[name]["avg"], data[name]["cnt"], data[name]["runs"]
    n, c = cfg["n"], cfg["color"]
    lab = f"{name} ({n} Läufe)"
    z = avg[:, 1] * 1e3
    zmax = max(zmax, z[-1])
    Tp, Dp = avg[:, 2], avg[:, 3] * 1e9
    # runningAver writes summed squared deviations -> divide by (active runs - 1);
    # band only where >= 3 runs are still alive (runs end at different heights)
    n_act = np.array([sum(r["z_end"] >= zz for r in runs) for zz in z])
    ok = n_act >= 3
    with np.errstate(invalid="ignore", divide="ignore"):
        sTp = np.where(ok, np.sqrt(avg[:, 4] / np.maximum(n_act - 1, 1)), np.nan)
        sDp = np.where(ok, np.sqrt(avg[:, 5] / np.maximum(n_act - 1, 1)) * 1e9, np.nan)
    if (~ok).any():
        for ax in (ax_T, ax_D):
            ax.axvline(z[np.argmax(~ok)], color=c, ls="--", lw=0.8, alpha=0.6)

    Tp_s, Dp_s = smooth(z, Tp), smooth(z, Dp)
    T_s_all.append(Tp_s)
    for ax, y, ys, s in ((ax_T, Tp, Tp_s, sTp), (ax_D, Dp, Dp_s, sDp)):
        ax.plot(z, y, color=c, **RAW)
        ax.plot(z, ys, color=c, lw=1.8, label=lab)
        ss = smooth(z, s)
        ax.fill_between(z, ys - ss, ys + ss, color=c, alpha=0.15, lw=0)

    # collisions: total reactive rate
    zc = cnt[:, 0] * 1e3
    ax_c.plot(zc, cnt[:, 9], color=c, **RAW)
    ax_c.plot(zc, smooth(zc, cnt[:, 9]), color=c, lw=1.8, label=lab)
    ax_r.plot(zc, cnt[:, 10], color=c, **RAW)
    ax_r.plot(zc, smooth(zc, cnt[:, 10]), color=c, lw=1.8, label=lab)

    # evaporation events per mm and run
    z_ev = np.concatenate([r["z_ev"] for r in runs]) if runs else np.array([])
    bins = np.arange(1.2, 7.0 + 1e-9, 0.1)
    h, _ = np.histogram(z_ev, bins=bins)
    ax_ev.step(bins[:-1], h / (n * 0.1), where="post", color=c, lw=1.6, label=lab)

    # per-run summary: end height and max n_Fe
    x = np.arange(1, n + 1) + (j - 0.5) * 0.3
    ax_run.scatter(x, [r["z_end"] for r in runs], color=c, s=30, label=f"{name}: Ende z")
    ax_n.scatter(x, [r["nfe_max"] for r in runs], color=c, marker="x", s=30)

ax_T.plot(cant["grid"] * 1e3, cant["T"], "k:", lw=1.2, label="$T_g$ (Cantera)")
lo = min(np.nanmin(t) for t in T_s_all + [cant["T"]])
hi = max(np.nanmax(t) for t in T_s_all)
ax_T.set_ylim(0.9 * lo, 1.15 * hi)

for ax in (ax_T, ax_D, ax_ev, ax_c, ax_r):
    ax.set_xlim(1.2, zmax)
    ax.set_xlabel("Höhe über Brenner z [mm]")
for ax in axs.flat:
    ax.grid(alpha=0.3)

ax_T.set(ylabel="Temperatur [K]", title="Partikeltemperatur $T_p$ (Band: ±1σ)")
ax_D.set(ylabel="$D_p$ [nm]", title="Partikeldurchmesser $D_p$ (Band: ±1σ)")
ax_ev.set(ylabel="Verdampfungen pro mm und Lauf", title="Verdampfungsereignisse")
ax_c.set(ylabel="Kollisionsrate [1/s]", title="Reaktive Kollisionen gesamt", yscale="log")
ax_r.set(ylabel="ratio [-]", title="ratio (count.csv, Spalte 11)")
ax_run.set(xlabel="Lauf-Nr.", ylabel="Ende z [mm] (●)", title="Pro Lauf: Ende z (●) und max. $n_{Fe}$ (×)")
ax_n.set_ylabel("max. $n_{Fe}$ (×)")
for ax in (ax_T, ax_D, ax_ev, ax_c, ax_r):
    ax.legend(fontsize=8)
ax_run.legend(fontsize=8, loc="lower right")

win_txt = f"gleitender Mittelwert {WIN:g} mm, Rohdaten blass" if WIN > 0 else "ungeglättet"
fig.suptitle(f"{CASE}: Evaporation mit Kelvin-Deckel vs. Kelvin + Tolman-Korrektur (Khrabry et al. 2025) – {win_txt}",
             fontsize=12)
fig.tight_layout()

os.makedirs(OUT_DIR, exist_ok=True)
out = os.path.join(OUT_DIR, f"{CASE}_kelvin-cap_vs_tolman")
fig.savefig(out + ".png", dpi=150)
fig.savefig(out + ".pdf")
print(f"Gespeichert: {out}.png / .pdf\n")

print(f"{'':22s}{'Läufe':>6s}{'Ende z [mm]':>16s}{'max n_Fe':>12s}{'Verd./Lauf':>12s}"
      f"{'Tp mean':>9s}{'Dp mean':>9s}{'Koll. ges.':>12s}")
for name in SETS:
    a, cn, runs = data[name]["avg"], data[name]["cnt"], data[name]["runs"]
    ze = np.array([r["z_end"] for r in runs])
    nm = np.array([r["nfe_max"] for r in runs])
    ev = np.array([len(r["z_ev"]) for r in runs])
    print(f"{name:22s}{len(runs):6d}{ze.mean():9.2f} ± {ze.std():4.2f}{nm.mean():8.1f} ± {nm.std():3.1f}"
          f"{ev.mean():12.1f}{a[:, 2].mean():9.0f}{a[:, 3].mean()*1e9:9.3f}{np.nanmean(cn[:, 9]):12.3e}")
print("(Tp in K, Dp in nm über die gemittelte Kurve; Koll. ges. = mittlere Gesamtrate in 1/s)")
print("Gestrichelte Linien in T_p/D_p: ab hier < 3 Läufe aktiv, kein σ-Band mehr.\n")

# Significance of the per-run differences (Welch t-test)
try:
    from scipy.stats import ttest_ind
    a_, b_ = (data[k]["runs"] for k in SETS)
    for lab_, key in (("Ende z", lambda r: r["z_end"]), ("Verdampfungen/Lauf", lambda r: len(r["z_ev"])),
                      ("Verd. pro mm", lambda r: len(r["z_ev"]) / max(r["z_end"] - 1.247, 1e-3)),
                      ("max n_Fe", lambda r: r["nfe_max"])):
        x_, y_ = [key(r) for r in a_], [key(r) for r in b_]
        t, pv = ttest_ind(x_, y_, equal_var=False)
        print(f"{lab_:20s} Deckel {np.mean(x_):7.2f} ± {np.std(x_, ddof=1):5.2f} | "
              f"Tolman {np.mean(y_):7.2f} ± {np.std(y_, ddof=1):5.2f} | p = {pv:.2f}")
except ImportError:
    pass

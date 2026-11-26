#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Vergleich RESULTS_original (23.11.2026, 30 MC-Läufe, Referenzstand)
#        vs. RESULTS          (26.11.2026, 1 MC-Lauf, mit erweiterter Evaporation)
#
# Die Evaporation in SRC/scripts/SPMC-0D-main-200326.py (evap_rates() und der
# Evaporationsblock in der MC-Schleife) wurde um folgende Modelle erweitert:
#
#  1. Hertz-Knudsen-Verdampfungsfluss
#       r = a * alpha * p_sat * K / sqrt(2*pi*m*k_B*T) * 4*pi*R^2
#     alpha = Verdampfungskoeffizient (YAML: evaporation.alpha = 1.0)
#
#  2. Sättigungsdampfdruck aus der NASA-Thermodynamik des Mechanismus
#       p_sat = 101325 * exp(-dG(T) / (R*T)),  dG = G_gas - G_flüssig
#     validiert: p_sat(Fe, 1811 K) ~ 3,5 Pa, p_sat(Fe, 3134 K) ~ 1 bar
#
#  3. Kelvin-Faktor (Krümmungseinfluss kleiner Cluster)
#       K = exp(4*sigma(T)*v_m / (max(d_p, d_min)*k_B*T))
#     gedeckelt bei d_min = 1 nm (YAML: evaporation.d_min_kelvin),
#     weil die Kapillarnäherung unter ~1 nm nicht mehr gilt
#
#  4. Verdampfende Spezies: Fe(l) -> Fe(g) mit Gewicht (1 - x),
#     FeO(l) -> FeO(g) mit Gewicht x, x = min(O/Fe, 1)
#
#  5. Stochastischer MC-Schritt: P = 1 - exp(-sum(r) * dt),
#     nur wenn n_Fe > 8
#
#  6. Energiebilanz: Das verdampfte Molekül nimmt latente Wärme (dH_vap)
#     UND seine fühlbare Wärme (m_ev * c_p * T_p) mit
#
# Details und Validierung: README_Evaporation_2026-11-24.md
# (Tolman-Korrektur nach Khrabry et al. 2025 ist dort nur vorgeschlagen,
#  aber noch nicht eingebaut)
# -----------------------------------------------------------------------------
"""Compare SPMC results: RESULTS_original vs. RESULTS (current run, Evaporation got updated (Kelvin factor, Hertz-Knudsen and )).

Usage: python SRC/scripts/compare_results.py [case] [n_runs_original] [n_runs_new] [--window-mm 0.1]
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("case", nargs="?", default="L500")
p.add_argument("n_orig", nargs="?", type=int, default=30)
p.add_argument("n_new", nargs="?", type=int, default=1)
p.add_argument("--window-mm", type=float, default=0.1, help="Glättungsfenster in mm (0 = aus)")
args = p.parse_args()
CASE, N_ORIG, N_NEW, WIN = args.case, args.n_orig, args.n_new, args.window_mm

BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SETS = {
    "original": dict(dir=os.path.join(BASE, "RESULTS_original", f"{CASE}_results"), n=N_ORIG, color="0.35"),
    "neu":      dict(dir=os.path.join(BASE, "RESULTS", f"{CASE}_results"), n=N_NEW, color="tab:red"),
}
OUT_DIR = os.path.join(BASE, "Plots", "Comparison")
SPECIES = {"O2": 2, "H2O": 4, "OH": 7, "O": 3, "H": 6, "H2": 5}
RAW = dict(lw=0.5, alpha=0.2)


def load(d):
    avg = np.genfromtxt(os.path.join(d, f"{CASE}-average.csv"), delimiter=",")
    cnt = np.genfromtxt(os.path.join(d, f"{CASE}-count.csv"), delimiter=",")
    cant = np.genfromtxt(os.path.join(d, f"{CASE}-cantera.csv"), delimiter=",", names=True)
    return avg, cnt, cant


def std_from_sum(s, n):
    # runningAver.py writes the summed squared deviations, not the std
    return np.sqrt(s / (n - 1)) if n > 1 else np.zeros_like(s)


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


data = {k: load(v["dir"]) for k, v in SETS.items()}

fig, axs = plt.subplots(2, 2, figsize=(12, 8.5), sharex=True)
(ax_T, ax_D), (ax_c, ax_r) = axs
T_smooth_all = []

for name, cfg in SETS.items():
    avg, cnt, cant = data[name]
    z = avg[:, 1] * 1e3
    Tp, Dp = avg[:, 2], avg[:, 3] * 1e9
    sTp, sDp = std_from_sum(avg[:, 4], cfg["n"]), std_from_sum(avg[:, 5], cfg["n"]) * 1e9
    lab = f"{name} ({cfg['n']} {'Läufe' if cfg['n'] > 1 else 'Lauf'})"
    c = cfg["color"]

    Tp_s, Dp_s = smooth(z, Tp), smooth(z, Dp)
    T_smooth_all.append(Tp_s)
    ax_T.plot(z, Tp, color=c, **RAW)
    ax_D.plot(z, Dp, color=c, **RAW)
    ax_T.plot(z, Tp_s, color=c, lw=1.8, label=f"$T_p$ {lab}")
    ax_D.plot(z, Dp_s, color=c, lw=1.8, label=lab)
    if cfg["n"] > 1:
        sT, sD = smooth(z, sTp), smooth(z, sDp)
        ax_T.fill_between(z, Tp_s - sT, Tp_s + sT, color=c, alpha=0.15, lw=0)
        ax_D.fill_between(z, Dp_s - sD, Dp_s + sD, color=c, alpha=0.15, lw=0)

    zc = cnt[:, 0] * 1e3
    ls = "-" if name == "neu" else "--"
    for i, (sp, col) in enumerate(SPECIES.items()):
        r = smooth(zc, cnt[:, col])  # smooth incl. zeros, then hide zeros for log axis
        r = np.where(r > 0, r, np.nan)
        ax_c.plot(zc, r, ls=ls, color=f"C{i}", lw=1.4, label=f"{sp} ({name})")
    ax_r.plot(zc, cnt[:, 10], color=c, **RAW)
    ax_r.plot(zc, smooth(zc, cnt[:, 10]), color=c, lw=1.8, label=lab)

# Gas temperature (cantera of the new run)
cant = data["neu"][2]
ax_T.plot(cant["grid"] * 1e3, cant["T"], "k:", lw=1.2, label="$T_g$ (Cantera)")

zmin = min(d[0][0, 1] for d in data.values()) * 1e3
zmax = max(d[0][-1, 1] for d in data.values()) * 1e3
for ax in axs.flat:
    ax.set_xlim(zmin, zmax)
    ax.grid(alpha=0.3)

# Limit T axis to the smoothed range so single spikes don't dominate
lo = min(np.nanmin(t) for t in T_smooth_all + [cant["T"]])
hi = max(np.nanmax(t) for t in T_smooth_all)
ax_T.set_ylim(0.9 * lo, 1.15 * hi)

ax_T.set_ylabel("Temperatur [K]")
ax_T.set_title("Partikeltemperatur (Band: ±1σ)")
ax_D.set_ylabel("$D_p$ [nm]")
ax_D.set_title("Partikeldurchmesser (Band: ±1σ)")
ax_c.set_yscale("log")
ax_c.set_ylabel("Kollisionsrate [1/s]")
ax_c.set_title("Reaktive Kollisionen (— neu, -- original)")
ax_r.set_ylabel("ratio [-]")
ax_r.set_title("ratio (count.csv, Spalte 11)")
for ax in (ax_c, ax_r):
    ax.set_xlabel("Höhe über Brenner z [mm]")
ax_T.legend(fontsize=8)
ax_D.legend(fontsize=8)
ax_c.legend(fontsize=7, ncol=2)
ax_r.legend(fontsize=8)
win_txt = f"gleitender Mittelwert {WIN:g} mm, Rohdaten blass" if WIN > 0 else "ungeglättet"
fig.suptitle(f"{CASE}: RESULTS_original vs. RESULTS – {win_txt}\n"
             f"(original {N_ORIG} MC-Läufe, neu {N_NEW} Lauf – statistisch nicht gleichwertig)", fontsize=11)
fig.tight_layout()

os.makedirs(OUT_DIR, exist_ok=True)
out = os.path.join(OUT_DIR, f"{CASE}_original_vs_new")
fig.savefig(out + ".png", dpi=150)
fig.savefig(out + ".pdf")
print(f"Gespeichert: {out}.png / .pdf\n")

# Summary table (raw data)
print(f"{'':10s}{'z-Bereich [mm]':>18s}{'Tp mean':>10s}{'Tp max':>10s}"
      f"{'Dp mean':>10s}{'Dp end':>10s}{'Koll. ges.':>12s}")
for name in SETS:
    avg, cnt, _ = data[name]
    z = avg[:, 1] * 1e3
    print(f"{name:10s}{z[0]:8.3f}-{z[-1]:<9.3f}{avg[:, 2].mean():10.0f}{avg[:, 2].max():10.0f}"
          f"{avg[:, 3].mean()*1e9:10.3f}{avg[-1, 3]*1e9:10.3f}{np.nanmean(cnt[:, 9]):12.3e}")
print("(Tp in K, Dp in nm, Koll. ges. = mittlere Gesamtrate in 1/s)")

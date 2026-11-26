# Evaporation im SPMC-Modell – Analyse vom 24.11.2026

Zusammenfassung einer Arbeitssitzung mit Claude Code: Testrechnung mit dem neuen
Evaporationsmodell, Prüfung der Werte, Ursachenanalyse, Code-Korrektur und Diskussion
der Monte-Carlo-Abfrage.

Präsentation (privat): https://claude.ai/artifact/7tBQpqpdcTkexsMDUW29wG

---

## 1. Ausgangsfrage

> „Führe eine Testrechnung durch und schaue nach, ob die Evaporationswerte nun Sinn machen.“

**Kurzantwort:** Die Evaporationswerte sind physikalisch korrekt. Der Cluster wächst
trotzdem nicht, weil kleine Cluster durch die Reaktionswärme zu heiß werden und deshalb
schneller verdampfen, als Fe aus dem Gas anlagert.

## 2. Evaporationsmodell (Stand im Code)

Datei: `SRC/scripts/SPMC-0D-main-200326.py`, Funktion `evap_rates()` und Evaporationsblock
in der MC-Schleife.

```
p_sat = 101325 · exp(−ΔG(T) / (R·T))          ΔG = G_gas − G_flüssig  (NASA-Thermo aus dem Mechanismus)
K     = exp(4·σ(T)·v_m / (max(d_p, d_min)·k_B·T))   d_min = 1 nm (part_params.yaml)
r     = a · α · p_sat · K / sqrt(2π·m·k_B·T) · 4π·R²
```

- Spezies: Fe(l) → Fe(g) mit Gewicht `1 − x`, FeO(l) → FeO(g) mit Gewicht `x`, `x = min(O/Fe, 1)`
- MC-Schritt: `P = 1 − exp(−Σr · dt)`, nur wenn `ni_FE_pMC > 8`
- Parameter in `YAML/part_params.yaml` → `evaporation: alpha, d_min_kelvin`

## 3. Testumgebung

Gerechnet wurde in einer Kopie im Scratchpad, die Projekt-`RESULTS/` blieben unberührt.

Das Skript erwartet relativ zum Arbeitsverzeichnis:
- `../part_params.yaml`
- `./YAML/<mechanismus>.yaml` (aus `L500.case`)
- `RESULTS/L500_results/L500-cantera.csv` und `L500-nuclest.csv`

Aufruf: `python SRC/scripts/SPMC-0D-main-200326.py L500.case <name>`

> **Hinweis:** `exe.sh` zeigt auf `SRC/SPMC-0D-main-200326.py`, die Datei liegt aber in
> `SRC/scripts/`. Außerdem liegt `part_params.yaml` in `YAML/`, gelesen wird `../part_params.yaml`.

Bei `awk`-Auswertungen `LC_ALL=C` setzen (sonst Dezimalkomma durch deutsche Locale).

## 4. Validierung der Thermodynamik

| Größe | Modell | Literatur |
|---|---|---|
| p_sat(Fe) bei 1811 K (Schmelzpunkt) | 3,5 Pa | ≈ 3,6 Pa |
| p_sat(Fe) bei 3134 K (Siedepunkt) | 0,99 bar | 1 bar |
| ΔH_vap(Fe) bei T_b | 349 kJ/mol | ≈ 340–350 kJ/mol |
| p_sat(FeO) / p_sat(Fe) | ≈ 1/100 – 1/200 | plausibel |

Kelvin-Faktor bei d = 1 nm: ≈ 58 bei 1811 K, ≈ 34 bei 2000 K, ≈ 6 bei 3000 K.

## 5. Erste Testrechnung (vor der Korrektur)

| Kennwert | Mit Evaporation | Ohne (Lauf 11.11.) |
|---|---|---|
| Ende | t = 1,95 ms, z = 3,3 mm | z = 6,3 mm |
| Größter Cluster | 8–10 Fe | ≈ 90 Fe bei 1,95 ms |
| Mittlere T_p | 1646 K | 1412 K |
| Anteil Zeit T_p > 2000 K | 17 % | 0,8 % |
| Fe angelagert / verdampft / geätzt (bis 1,95 ms) | 60 / 41 (32 Fe, 9 FeO) / 22 | 134 / – / 50 |

- Abbruch durch OH-Ätzen (< 6 Fe), nicht durch Verdampfung.
- Jede Verdampfung kühlt den Cluster um ≈ 1300–1600 K. Bei 8–10 Atomen ist ΔH_vap ≈ 380 kJ/mol
  etwa so groß wie die gesamte innere Energie, das ist also energetisch korrekt.
- T_p fällt nach Verdampfung teils unter T_g (z. B. 322 K bei T_g = 737 K).

## 6. Ursache: Übersättigung bezogen auf T_p ≪ 1

| | Wert |
|---|---|
| Fe-Partialdruck im Gas (z = 2–5 mm) | ≈ 0,6–1,7 Pa |
| p_sat · K bei T_p = 1650 K | ≈ 30 Pa |
| p_sat · K bei T_p = 2000 K | ≈ 1260 Pa |

Treiber sind die hohen Clustertemperaturen: Mit Bulk-c_p = 400 J/(kg·K) heizt eine einzelne
Reaktion einen Cluster aus 8–12 Atomen um mehrere hundert bis ≈ 3000 K (O-Oxidation) auf.
Das ist eine Folge des Kollisions-/Energiemodells, nicht der Evaporationsformel.

**Korrigierte Einschätzung:** Die Fe-Aktivität auf oxidierter Oberfläche ist mit der
Flächengewichtung `(1 − x)` in Ordnung, da Fe im Gleichgewicht mit FeO eine Aktivität von
≈ 1 hat. Kein Änderungsbedarf.

## 7. Code-Korrektur (umgesetzt)

Nur im Evaporationsblock von `SRC/scripts/SPMC-0D-main-200326.py`:

1. **Energiebilanz:** Das verdampfte Molekül nimmt zusätzlich zur latenten Wärme seine
   fühlbare Wärme mit. Vorher blieb sie im Rest-Cluster → T nach Verdampfung ≈ 10–15 % zu hoch.
   ```python
   dH  = (g.thermo.h(T_pMC) - c.thermo.h(T_pMC)) / 1e3 / N_A
   dH += m_ev[k] * C_pMC * T_pMC      # + sensible energy carried away by the leaving molecule
   ```
2. **Ausgabe:** Letzte Spalte der `ev`-Zeilen zeigt jetzt `T_ev0` (T direkt vor der
   Verdampfung) statt `Tbefore_pMC` (T vor der vorhergehenden Kollision).

Nicht geändert: Schwelle `ni_FE_pMC > 8`, Kollisions-/Energiemodell.

**Handrechnung zur Kontrolle:** 9 Fe + 6 O verlieren ein Fe bei 2833 K:
ΔH ≈ 5,9·10⁻¹⁹ J + m·c_p·T ≈ 1,05·10⁻¹⁹ J → T_nach = 1196 K. Simulation: 1196,7 K.

## 8. Läufe nach der Korrektur

| Lauf | Ende bei z | max. Fe | angelagert | verdampft Fe/FeO | geätzt |
|---|---|---|---|---|---|
| test2 | 1,5 mm | 8 | 1 | 0 / 0 | 4 |
| test3 | 2,2 mm | 12 | 26 | 15 / 5 | 9 |
| test4 | 3,9 mm | 10 | 75 | 40 / 12 | 26 |
| test5 | 3,5 mm | 11 | 66 | 39 / 8 | 22 |
| test6 | 1,5 mm | 9 | 4 | 1 / 0 | 6 |

- Alle Läufe enden mit 5 Fe durch OH-Ätzen.
- Keine negativen Energien; Abkühlung ≈ 1600–1660 K pro Verdampfung.
- Einzelne MC-Läufe streuen stark → immer mehrere Realisierungen vergleichen.

## 9. Liegt es an der Poisson-Abfrage?

Frage: Verhindert `if random() < 1.0 - mt.exp(-sum(r_ev)*dt)` das Wachstum? Was wäre mit
einer Summenbildung bis zu einem Schwellenwert?

**Antwort: Nein.** Jede korrekte Methode liefert im Mittel gleich viele Verdampfungen,
nämlich ∫ r_ev(T_p(t)) dt.

| Variante | Bewertung |
|---|---|
| `1 − exp(−r·dt)` (aktuell) | Exakt bei konstanter Rate pro Schritt. Max. 1 Ereignis pro Schritt → eher zu wenig, nie zu viel Verdampfung. |
| `Σ r·dt ≥ −ln(u)`, Zufallsschwelle, Reset nach Ereignis | Standardverfahren für zeitabhängige Raten, statistisch gleichwertig. |
| `Σ r·dt ≥ 1`, feste Schwelle (wie `coll_comp`) | Löst oft erst nach dem Abkühlen aus → `E_pMC > dH` scheitert → weniger Verdampfung. Wachstum dann nur durch künstliches Gedächtnis, physikalisch nicht korrekt. |

Vorgeschlagener (noch nicht durchgeführter) numerischer Test: Variante mit Zufallsschwelle in
einer Kopie des Skripts, zusätzlich Σ r_ev·dt als erwartete Ereigniszahl mitzählen,
je 4 Läufe vergleichen.

## 10. Fazit und offene Punkte

- **Gesichert:** p_sat, ΔH_vap, Kelvin-Faktor, Energiebilanz und Zufallsabfrage sind konsistent.
- **Ursache des Nicht-Wachstums:** heiße kleine Cluster im Kollisions-/Energiemodell.
- **Nächste Schritte:**
  - Kollisions-/Energiemodell prüfen (c_p kleiner Cluster, Reaktionswärme pro Ereignis)
  - Sensitivität mit `evaporation.alpha: 0.1` (reine Konfigurationsänderung)
  - Konsistenz Schwelle `ni_FE_pMC > 8` (Evaporation) vs. Abbruch `ni_FE_pMC < 6` (Ätzen)
  - Gültigkeit: Hertz-Knudsen + Kelvin unter ≈ 1 nm bzw. 20 Atomen nur noch näherungsweise
  - Optional: Tolman-Korrektur per Schalter einbauen (siehe Abschnitt 11). Erwarteter Effekt
    auf das Wachstum gering.
  - Ursache der heißen Cluster und Lösungsansätze: siehe Abschnitt 12 (zuerst Wärmekapazität
    pro Atom statt pro Masse).

## 11. Kelvin-Korrektur bei kleinen Clustern und Tolman-Korrektur

### Was der Code heute macht

```python
K = mt.exp(4.0*σ*(m/rho) / (max(2.0*R, d_min_K)*k_B*T))     # d_min_kelvin = 1 nm
```

Die Kelvin-Korrektur wird bei kleinen Clustern **nicht abgeschaltet**, sondern **gedeckelt**:
Unter 1 nm wird so gerechnet, als wäre der Cluster genau 1 nm groß. Alle kleinen Cluster
bekommen also denselben Kelvin-Faktor (bei 2000 K ≈ 35). Ohne Deckel würde er bei 9 Atomen
auf über 300 steigen, weil die klassische Kelvin-Gleichung eine konstante Oberflächenspannung
annimmt, was bei wenigen Atomen nicht mehr stimmt. Die 1-nm-Grenze ist eine willkürliche Notlösung.

### Was das Paper sagt

Khrabry, Hoffenberg, Kaganovich, Barsukov, Graves (2025): *The importance of accounting for the
Tolman correction to surface tension for nucleation and growth modeling of Fe clusters.*
Aerosol Sci. Technol. 59(6), 743–764, DOI 10.1080/02786826.2025.2452895

- Gibbs-Energien von Fe-Clustern mit 2–100 Atomen aus Molekulardynamik (MD) bestimmt.
- Kleine Cluster haben eine deutlich **kleinere Oberflächenspannung** als eine flache Fe-Oberfläche.
- Das beschreibt die Tolman-Korrektur gut (Gl. 15), außer bei Di- und Trimeren:
  ```
  σ_i = σ / (1 + 2·(δ_T/r_W)·i^(-1/3))         i = Atomzahl, r_W = Wigner-Seitz-Radius
  ```
- Fits (Gl. 16, 17), gültig etwa 1600–2400 K:
  ```
  σ_flach,MD   = −0,00152·T + 4,54   N/m
  σ·δ_T/r_W    = −0,0021 ·T + 5,46   N/m
  ```
- Mit r_W(Fe, ρ = 7000 kg/m³) = 0,147 nm: **δ_T ≈ 0,15 nm (1600 K) … 0,07 nm (2400 K)**, positiv.
- Mit Tolman stimmen Kondensationsverlauf, Clustergröße und Verteilungsbreite gut mit dem vollen
  MD-Modell überein; mit konstantem σ dagegen nicht.

### Nachgerechnet: Kelvin-Faktor für 8–12 Fe-Atome

Diskret aus der Änderung der Oberflächenenergie beim Verdampfen eines Atoms (Paper Gl. 7):
`K = exp((E_s(i) − E_s(i−1)) / (k_B·T))`, `E_s(i) = σ_i · 4π r_W² · (i^(2/3) − 1)`.

| T | heute (1-nm-Deckel) | Tolman: exp. σ + δ aus Paper | Tolman: σ und δ aus Paper |
|---|---|---|---|
| 1650 K | 101 | 83–116 | 93–130 |
| 2000 K | 34 | 36–49 | 18–23 |
| 2300 K | 18 | 25–34 | 7–9 |

Zum Vergleich ohne Deckel und ohne Tolman: 1060–3080 (1650 K), 210–475 (2000 K).

### Bewertung

- Tolman ist **physikalisch besser begründet** als der 1-nm-Deckel (MD-Daten für Fe).
- Für unsere Clustergröße liegt der 1-nm-Deckel **zufällig nahe** am Tolman-Ergebnis
  (Unterschied höchstens etwa Faktor 2).
- Das Grundproblem (zu heiße Cluster, Abschnitt 6) löst Tolman nicht.

### Haken beim Einbau

1. **Temperaturbereich:** Die Fits gelten nur für ≈ 1600–2400 K. Oberhalb ≈ 2600 K wird δ_T ≤ 0,
   oberhalb ≈ 2990 K sogar σ < 0. Im SPMC erreichen Cluster bis 3600 K → T für den Fit auf
   1600–2400 K begrenzen.
2. **Konsistenz mit dem Dampfdruck:** σ_flach des Papers gehört zum MD-Eisen mit eigenem Dampfdruck
   (L = 423 kJ/mol, T_boil = 3126 K). Wir nutzen den Cantera-/Experiment-Dampfdruck. Deshalb
   nur δ_T/r_W aus dem Paper übernehmen und das bisherige σ(T) behalten.
3. **Nur reines Fe:** Für FeO gibt es keine Tolman-Daten → dort 1-nm-Deckel beibehalten.

### Implementierungsvorschlag (noch nicht umgesetzt)

`YAML/part_params.yaml`, Abschnitt `evaporation`:
```yaml
kelvin_mode: tolman      # 'cap' = bisher (max(d, d_min_kelvin)), 'tolman' = Khrabry et al. 2025
tolman_a: -0.0021        # sigma*delta/rW = a*T + b  (N/m/K, N/m)
tolman_b: 5.46
tolman_Tmin: 1600.0      # Gültigkeitsbereich des Fits; T wird darauf begrenzt
tolman_Tmax: 2400.0
```

In `evap_rates()` nur für Fe (k = 0):
```python
Tt  = min(max(T, Tmin), Tmax)
sig = s0 + ds*(T - Tr)                                   # bisheriges sigma(T)
dr  = (a*Tt + b) / (-0.00152*Tt + 4.54)                  # delta_T / r_W aus dem Paper
i   = (R / r_W)**3                                       # äquivalente Atomzahl
Es  = lambda n: sig/(1 + 2*dr*n**(-1/3)) * 4*PI*r_W**2 * (n**(2/3) - 1)
K   = mt.exp((Es(i) - Es(i-1)) / (k_B*T))
```
`kelvin_mode: cap` reproduziert exakt das bisherige Verhalten. Danach 4 Testläufe und Vergleich
mit test2–test6.

### Lohnt sich der Einbau?

- **Für das aktuelle Problem: kaum.** Für 8–12 Atome ändert sich die Verdampfungsrate höchstens
  um etwa Faktor 2. Damit das Partikel bei 2000 K wächst, müsste die Verdampfung aber etwa
  20- bis 800-mal schwächer sein.
- **Dafür spricht:** physikalisch begründet statt willkürlicher 1-nm-Grenze; bei größeren Clustern
  (50–100 Atome) stärkerer Effekt (2000 K: K = 26 → 13 bzw. 13 → 9); wenig Aufwand (Schalter).
- **Dagegen spricht:** neue Unsicherheiten (Temperaturbereich, nur Fe, Dampfdruck-Konsistenz).
- **Empfehlung:** Zuerst die Ursache der heißen Cluster angehen (Abschnitt 12), Tolman danach
  als Option einbauen, spätestens wenn das Partikel wieder wächst oder für eine Veröffentlichung.

## 12. Warum werden die kleinen Cluster so heiß, und wie lässt sich das lösen?

### Einfach erklärt

Ein brennendes Streichholz an einem großen Stein erwärmt ihn kaum; dasselbe Streichholz an einem
Sandkorn macht es glühend heiß. Unser Cluster aus ≈ 10 Atomen ist noch viel kleiner.

Im Modell passiert bei jeder Reaktion:
1. Ein Gasmolekül trifft den Cluster und reagiert (z. B. O₂ oder ein O-Atom mit Fe).
2. Die **gesamte** Reaktionswärme landet im Cluster.
3. Weil der Cluster so wenig Masse hat, steigt seine Temperatur stark:

| Reaktion | Reaktionswärme (Hr) | Erwärmung, Cluster aus 10 Fe |
|---|---|---|
| Oxidation mit O₂ | 160 kJ/mol | ≈ +700 K |
| Oxidation mit O-Atom | 658 kJ/mol | ≈ +2900 K |

**Abkühlen dauert:** Wärme geht nur über Stöße mit Gasmolekülen weg. Jeder Stoß kühlt nur um
≈ 50–100 K, es braucht also 10–30 Stöße. In dieser heißen Phase dampfen Atome ab, schneller als
neues Eisen ankommt.

Dass winzige Cluster durch einzelne Reaktionen kurz sehr heiß werden, ist grundsätzlich echte
Physik. Zwei Stellen im Modell übertreiben das aber vermutlich.

### Was im Modell vermutlich zu stark ist

**1. Wärmekapazität wird über die Masse gerechnet.**
- Im Code: `E_pMC = C_pMC · M_pMC · T_pMC` mit `C_p = 400 J/(kg·K)` (passt für reines Fe).
- Physikalisch nimmt jedes Atom etwa gleich viel Wärme auf, egal ob Fe oder O
  (Dulong-Petit: ≈ 3·k_B pro Atom).
- Ein O-Atom ist ein Drittel so schwer wie Fe und bekommt im Modell deshalb nur ein Drittel der
  Wärmekapazität, die es tatsächlich hat.
- Unsere Cluster sind meist stark oxidiert (z. B. 10 Fe + 8 O). Für sie ist die Wärmekapazität
  im Modell ≈ 40 % zu klein, die Temperaturspitzen entsprechend zu hoch:

| Cluster 10 Fe + 8 O | Modell heute (pro Masse) | pro Atom (Dulong-Petit) |
|---|---|---|
| Erwärmung bei O₂-Oxidation | +580 K | +360 K |
| Erwärmung bei O-Atom-Oxidation | +2400 K | +1470 K |

**2. Die gesamte Reaktionswärme bleibt im Cluster.**
Bei vielen Reaktionen fliegt ein Produkt weg und nimmt Energie mit (bei O₂ ein O-Atom, bei der
Reduktion mit H ein OH). Das Modell gibt trotzdem 100 % der Wärme an den Cluster.

### Lösungsansätze (vom einfachsten zum aufwendigsten)

1. **Wärmekapazität pro Atom statt pro Masse.** Kleine Änderung an wenigen Stellen: Energie aus
   der Atomzahl (`3·k_B` pro Fe- und O-Atom) statt aus der Masse. Bei reinem Fe fast keine
   Änderung, bei oxidierten Clustern ≈ 40 % kleinere Spitzen. **Empfohlener erster Schritt.**
2. **Nur einen Teil der Reaktionswärme in den Cluster geben.** Neuer YAML-Parameter (Anteil der
   Reaktionswärme, der im Partikel bleibt); Wert muss aus der Literatur begründet werden.
3. **Reaktionswärmen prüfen.** Besonders die 658 kJ/mol für O-Atome sind sehr hoch; Quelle prüfen.
4. **Cluster-Cluster-Stöße (Koagulation) ergänzen.** Khrabry et al. (2025, Abb. 11) zeigen, dass
   kleine Fe-Cluster vor allem durch Zusammenstöße mit anderen Clustern wachsen, nicht nur Atom für
   Atom. Im SPMC wächst der Cluster bisher nur atomweise. Größte Änderung, aber möglicherweise der
   eigentliche Grund, warum echte Partikel trotzdem wachsen.

**Vorschlag:** Punkt 1 einbauen und mit 4 Testläufen prüfen, ob der Cluster dann wächst.
Reicht das nicht, sind Punkt 2 oder 4 die nächsten Kandidaten.

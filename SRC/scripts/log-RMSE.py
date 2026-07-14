# storing log-RMSE for evaluation

import pandas as pd
import os

#### Read the last columns of all CSV files and store them in a new CSV file ####
# alle Testnummern
tests = range(1, 12)

# Speicher für die letzten Spalten
result = {}

# alle Dateien durchgehen
for j in tests:
    files = [
        f"./ArrheniusRates/rates_calc_test{j}.csv",
        f"./ArrheniusRates/rates_calc_test{j}.1.csv",
        f"./ArrheniusRates/rates_calc_test{j}.2.csv",
        f"./ArrheniusRates/rates_calc_test{j}.3.csv"
    ]

    for file in files:
        # check if file exists
        if os.path.exists(file):

            # CSV einlesen
            data = pd.read_csv(file, header=None)

            # letzte Spalte auswählen
            letzte_spalte = data.iloc[:, -1]

            # Spaltenname aus Dateiname erzeugen
            name = os.path.basename(file).replace(".csv", "")

            # speichern
            result[name] = letzte_spalte


# alles zu einer großen Tabelle zusammenfügen
output = pd.DataFrame(result)

# speichern
output.to_csv("./ArrheniusRates/LastColumns_all_tests.csv", index=False)

#print(output)

#### evaluate the minimum values ####

# Datei einlesen
data = pd.read_csv("./ArrheniusRates/LastColumns_all_tests.csv")

# kleinster Wert pro Zeile
min_values = data.min(axis=1)

# Name der Datei/Testspalte mit dem kleinsten Wert
min_test = data.idxmin(axis=1)

# Ergebnis speichern
result = pd.DataFrame({
    "Minimum": min_values,
    "Test": min_test
})

result.to_csv("./ArrheniusRates/min_values_with_test.csv", index=False)

print(result)
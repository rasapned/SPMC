
rates_calc.csv = [A, b, Ea, alpha_k, m_k, E_k, logRMSE]

### Tests 1–11

| Test      | 1    | 2     | 3     | 4     | 5    | 6     | 7     | 8     | 9    | 10    | 11    |
| --------- | ---- | ----- | ----- | ----- | ---- | ----- | ----- | ----- | ---- | ----- | ----- |
| lin       | x    | x     | x     | x     |      |       |       |       |      |       |       |
| log       |      |       |       |       | x    | x     | x     | x     |      |       |       |
| method    |      |       |       |       |      |       |       |       | x    | x     | x     |
| b_switch  | True | True  | False | False | True | True  | False | False | True | True  | False |
| blockE_a  | True | False | True  | False | True | False | True  | False | True | False | True  |

---

### Test 12

| Test      | 12    |
| --------- | ----- |
| lin       |       |
| log       |       |
| method    | x     |
| b_switch  | False |
| blockE_a  | False |

### T_ref
---       
| T_ref  |_______________________________Test________________________________________
| ------ | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---  | ---  | ---  |
| 1000 K | 1.1 | 2.1 | 3.1 | 4.1 | 5.1 | 6.1 | 7.1 | 8.1 | 9.1 | 10.1 | 11.1 | 12.1 |
| 2000 K | 1.2 | 2.2 | 3.2 | 4.2 | 5.2 | 6.2 | 7.2 | 8.2 | 9.2 | 10.2 | 11.2 | 12.2 |
| 3000 K | 1.3 | 2.3 | 3.3 | 4.3 | 5.3 | 6.3 | 7.3 | 8.3 | 9.3 | 10.3 | 11.3 | 12.3 |

The Größe logRMSE in rates_calc.csv gibt auskunft über die Güte des Fittings, je kleiner der Wert desto besser. Für jeden Test wird mittels log-RMSE.py diese Größe ausgewertet und mit anderen Bedingungen ausgetestet, schlißlich werden die besten Bedingungen für das Fitting ermittelt und in min_values_with_test.csv gespeichert.

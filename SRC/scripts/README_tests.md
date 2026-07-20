
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

The logRMSE value in rates_calc.csv indicates the quality of the fitting; the smaller the value, the better the fit. For each test, this metric is evaluated using log-RMSE.py under different conditions. The optimal fitting conditions are then identified and stored in min_values_with_test.csv.

Conclusion: The results in min_values_with_test.csv were reviewed, and it was found that the corresponding fits were not visually satisfactory. One of the main issues is the activation energy, which is controlled by blockE_a. In several fits, the activation energy falls below 1, which is not a physically meaningful value. Through manual evaluation, the best fit was identified as Test 12.
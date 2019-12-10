# Ratio calculations using perebor.py.

**div** -- number of divisions on each side, i.e. genus = div^dim

**cdim** -- minimum allowed dimension for neighbour cubes intersection

## 2D-curves

|div|entr-exit      |cdim|#paths|l1   |l2   |linf |
|---|---------------|----|------|-----|-----|-----|
|2  |side           |0   |1     |9.000|6.000|6.000|
|3  |side           |0   |5     |9.800|5.666|**4.500**|
|3  |diag           |0   |2     |10.00|5.666|5.333|
|4  |side           |1   |5     |9.000|5.909|5.4545|
|4  |side           |0   |162   |9.000|5.909|5.294|
|5  |side           |1   |43    |10.00|**5.589**|5.333|
|5  |side           |0   |24850 |     |5.589|     |
|5  |diag           |1   |11    |10.66|6.0  |5.77 |
|5  |diag           |0   |659   |     |6.0  |     |
|5  |(0,0)->(1/2,1) |0   |557   |11.85|6.535|6.25 |
|6  |side           |1   |897   |9.000|5.62 |5.294|
|6  |side           |0   |      |     |     |     |

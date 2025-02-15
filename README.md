# iED.jl
repo for exact diagonalization in cmt

## Spin Model
### Heisenberg Chain
- Symmtery: $U(1)$ spin, translation symmetry.
- Results: calculating ground state with quantum number

$$m=0,\quad k=0\ \mathrm{or}\ \pi$$

gives:
| N    | $E_A$               | CPU-time  |
|------|--------------------|-----------|
| 6   |  -0.7171292729553314 | 0.000602s    |
| 12   |  -0.6989492431204336 | 0.003069s    |
| 16  |  -0.6963935225385482 | 0.031542s    |
| 18 |  -0.6957082826129863 | 0.366234s    |
| 20 |  -0.6952193264938217 | 3.278323s   |

Accuracy within 1e-15 compared to the Bethe-ansatz.





# Kopriva-benchmark

Toy project for DG method.


### Timeline
- [x] Quadrature nodes and weights
- [x] Nodal representation
  - [x] Lagrange interpolation
  - [x] Differentiation matrix
- [x] DG approximation
- [x] Time evolution
- [ ] Splitting domain
- [ ] Flux
- [ ] ...
- [ ] Completed!


#### Benchmark for $N=6$ Gauss points (Kopriva p.67 Table 3.1)

* Gauss
  
| $j$ | $x_j$ | $w_j$ |
|:---|---:|---:|
| 0 | -0.949107912342758  | 0.129484966168870 |
| 1 | -0.741531185599394  | 0.279705391489277 |
| 2 | -0.405845151377397  | 0.381830050505119 |
| 3 |  0.000000000000000  | 0.417959183673469 |

* Gauss-Lobatto

| $j$ | $x_j$ | $w_j$ |
|:---|---:|---:|
| 0  | -1.000000000000000 | 0.047619047619048 |
| 1  | -0.830223896278567 | 0.276826047361566 |
| 2  | -0.468848793470714 | 0.431745381209863 |
| 3  |  0.000000000000000 | 0.487619047619048 |




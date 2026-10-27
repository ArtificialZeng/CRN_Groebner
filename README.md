# CRN Steady-State Ideals (Gröbner Basis Computation)

A collection of Macaulay2 scripts that build steady-state ideals for biochemical reaction-network models and compute their reduced Gröbner bases.  Each script follows the same structure:

1. declare a rational-function coefficient field in the rate constants and conserved totals,
2. create a polynomial ring ordered lexicographically,
3. define the mass-action steady-state equations,
4. add linear conservation laws,
5. generate the ideal and compute a reduced Gröbner basis.

In some examples, certain rate constants are assigned values a priori.

---

## Repository layout

| file            | model                              | key reference                                      |
| --------------- | ---------------------------------- | -------------------------------------------------- |
| `ex_ERK.m2`              | ERK core         | Loman *et al.* 2023  |
| `ex_TGFb.m2`              | TGF-β receptor network  | Samal *et al.* 2016                                |
| `ex_GRN.m2`              | Gene Regulatory Network  | Conradi *et al.* 2017  |
| `ex_minimal_mapk.m2`              | Minimal MAPK Cascade  | Hadač *et al.* 2017  |
| `ex_Wnt.m2`              | WNT signalling pathway  | Gross *et al.* 2016  |
| `ex_toy.m2`              | Toy Network  | Berra *et al.* 2024  |
| `ex_toy2.m2`              | Another Toy Network  | Berra *et al.* 2024 |


## Requirements

- Macaulay2 ≥ 1.21


## Proposed positive regular diagnostic

[Original proof and exact Python checks](counterexamples/positive_regular/) for an indecomposable bimolecular family with a positive globally stable state and regular complex zeros, yet no species-permutation shape basis. This source-specific note credits the earlier counterclaim and the standard mixed-coordinate theorem; its full-text attribution targets arXiv v2. The existing Macaulay2 examples are unchanged.

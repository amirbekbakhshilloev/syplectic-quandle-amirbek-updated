# Symplectic Quandles Distinguish the Allen-Swenberg Links

Code and data for the paper *Symplectic Quandles Distinguish the Allen-Swenberg Links*
by Amirbek Baxshilloyev (arXiv:2508.18323). The paper source and PDF are in `paper/`.

Version 1 of the paper computed the enhanced quandle counting polynomial with the number of
distinct arc colors in place of the size of the image subquandle. This repository computes the
correct invariant, reproduces the version 1 numbers for comparison, and contains all computations
of the revised paper.

The link presentations of H, L1 and L2 were transcribed from the Mathematica notebooks of
version 1: https://github.com/harryfrences-svg/Horizon-Research

## Requirements

Python 3.10 or newer, the packages in `requirements.txt`, and a C compiler (gcc, cc or clang).
The file `symq/colorsolve.c` is compiled automatically the first time it is needed.
`numba` and `snappy` are only used in the `links` stage.

```
pip install -r requirements.txt
python3 reproduce.py
```

## Running

```
python3 reproduce.py                  # all stages, about one hour
python3 reproduce.py links v1         # selected stages
python3 reproduce.py --full           # larger prime ranges and more quandles, several hours
python3 big_primes.py 29,31           # p = 29, 31 for L1 in Table 2, about 45 minutes
```

Each stage writes `results/<stage>.json`. The files in `results/` are the ones used in the paper.

| stage    | paper                                         | contents |
|----------|-----------------------------------------------|----------|
| `links`  | Section 5.2                                   | components, linking numbers, Alexander polynomials, planar reconstruction, Kauffman bracket spans of H, L1, L2, L3 |
| `v1`     | Table 1, Example 5.1, Proposition 5.2         | version 1 arc color counts against the corrected polynomial, Reidemeister II example, formula for H |
| `primes` | Table 2, Theorem 5.3                          | counts over (Z_p)^2 |
| `family` | Table 3, Lemma 5.4, Corollary 5.5             | transfer matrix counts for L_n, n = 1, ..., 8 |
| `closed` | Theorems 5.6, 5.7, 5.8, equation (2)          | minimal polynomials, closed formulas, enhanced polynomial of every L_n |
| `other`  | Table 4                                       | other symplectic quandles |

## Files

- `symq/links.py`: Wirtinger presentations. A relation `(k, a, b, c)` means `x_a |> x_b = x_c` at crossing `c_k`. `allen_swenberg(n)` builds L_n.
- `symq/quandles.py`: rings Z/n and GF(p^k), symplectic quandles `x |> y = x + <x,y> y`.
- `symq/coloring.py`, `symq/colorsolve.c`: enumeration of Hom(Q(L), T), in Python and in C.
- `symq/invariants.py`: generated subquandle, corrected enhanced polynomial, version 1 arc color count.
- `symq/reidemeister.py`: the Reidemeister II example.
- `symq/topology.py`, `symq/planar.py`, `symq/jones.py`: topological checks of the link data.
- `symq/transfer.py`, `symq/closedform.py`: transfer matrices along the family L_n and exact closed forms.
- `run_enhanced_allN.py`: enhanced polynomial of L_n for all n via the lattice of subquandles.
- `scan_other.py`: list of further symplectic quandles.
- `paper/`: LaTeX source, figures and PDF of the revised paper.

## License

MIT, see `LICENSE`.

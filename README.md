# chol_incomplete
Matlab MAX code for incomplete Cholesky decomposition of dense matrices


## Installation

- Compile MEX files, by running make. This produces DLL files.
- Copy all DLL and *.m files somewhere into your Matlab path.

## How to use it

The Matlab help for CHOL_INCOMPLETE should be clear enough, and the code is well-documented. The incomplete Cholesky decomposition works like a normal Cholesky decomposition, but you stop after less than n iterations. Furthermore, you reorder the rows/columns as you go, s.t. you always choose the next column to have the largest possible pivot element (diagonal of L).

If you want to know more, you might want to check

```bibtex
@article{Fine:01,
  author      = {Fine, S. and Scheinberg, K.},
  title       = {Efficient {SVM} Training using Low-Rank Kernel Representations},
  journal     = JMLR,
  volume      = {2},
  pages       = {243--264},
  year        = {2001}
}
```

Once L is computed, you can use L*L' as approximation to P'*A*P, where A is the matrix to be approximated (usually a kernel matrix), and P is a permutation matrix.

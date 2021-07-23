# SPOOF - SParse matrices in Object-Oriented Fortran

*Work in progress.*

This library is a modern, object-oriented Fortran library for sparse matrix computations. It largely draws on the [SPARSKIT](https://www-users.cse.umn.edu/~saad/software/SPARSKIT/) library, adding an object-oriented interface to improve usability.

Currently implemented:
- CSR matrix format: `type(CSRMatrix)`. Creation from and conversion to dense matrices.
- Diagonal matrix format: `type(DiagonalMatrix)`. Creation from and conversion to dense matrices.
- Matrix by vector multiplication for CSR and diagonal matrices: `out_matrix = matrix%multiply(vector)`.

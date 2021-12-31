# SPOOF - SParse matrices in Object-Oriented Fortran

*Work in progress.*

This library is a modern, object-oriented Fortran library for sparse matrix computations. It largely draws on the [SPARSKIT](https://www-users.cse.umn.edu/~saad/software/SPARSKIT/) library, adding an object-oriented interface to improve usability.

Currently implemented:
- CSR matrix format: `type(CSRMatrix)`. Creation from and conversion to dense matrices.
- Diagonal matrix format: `type(DiagonalMatrix)`. Creation from and conversion to dense matrices.
- Matrix by vector multiplication for CSR and diagonal matrices: `out_matrix = matrix%multiply(vector)`.
- Only double precision for the moment.

## How to use

### `fpm` - Fortran Package Manager

SPOOF can be built using the Fortran Package Manager ([fpm](https://github.com/fortran-lang/fpm)). To use in your fpm compatible project, simply include in your `fpm.toml` file as a dependency:

```toml
[dependencies]
spoof = { git = "https://github.com/samharrison7/spoof" }
```

Or you can clone the repo and build (`fpm build`) or install locally (`fpm install`) the library yourself using fpm.

### Grab the files

SPOOF is only a few files - just grab everything in `src/` and include in your own build process.

## Example

```fortran
program main
    use Spoof
    implicit none

    type(CSRMatrix)         :: csr_matrix
    type(DiagonalMatrix)    :: dia_matrix
    double precision        :: dense(10,10)
    double precision        :: vec(10)

    dense = 42.0D+00                    ! 10x10 matrix filled with 42
    vec = 10.0D+00                      ! 10-element vector filled with 10

    ! Create the matrices
    csr_matrix = CSRMatrix(dense)
    dia_matrix = DiagonalMatrix(dense)

    ! Multiply them by vec
    vec = csr_matrix%multiply(vec)      ! vec is filled with 4200
    vec = dia_matrix%multiply(vec)      ! vec is filled with 176400
end program
```

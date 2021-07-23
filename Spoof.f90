!-----------------------------------------------------------------------!
!> SPOOF - SParse matrices in Object-Oriented Fortran                   !
!> --------------------------------------------------                   !
!> A modern, object-oriented Fortran library for sparse                 !
!> matrix computations. Based on the SPARSKIT library:                  !
!> https://www-users.cse.umn.edu/~saad/software/SPARSKIT/               !
!>                                                                      !
!> Authors: Sam Harrison (sharrison@ceh.ac.uk)                          !
!> Repository: https://github.com/samharrison7/spoof                    !
!> License: MIT License,                                                !
!>          https://github.com/samharrison7/spoof/blob/develop/LICENSE  !
!-----------------------------------------------------------------------!
module Spoof
    use sparskit
    implicit none

    ! Double precision kind
    integer, parameter :: d = selected_real_kind(15, 307)

    !> Compressed Sparse Row matrix
    type, public :: CSRMatrix
        real(d), allocatable    :: data(:)          !! The non-zero elements of the matrix, stored as a contiguous array
        integer, allocatable    :: col_ind(:)       !! The indices of the columns of the values in the data array
        integer, allocatable    :: row_ptr(:)       !! The index in the data array and col_ind where a given row starts
        integer                 :: ncol             !! Number of columns in the matrix
        integer                 :: nrow             !! Number of rows in the matrix
      contains
        procedure               :: to_dense => to_dense_csr         !! Convert CSR to dense matrix
        procedure               :: to_array => to_dense_csr         !! Alias of to_dense, for SciPy-like interface
        procedure               :: multiply => multiply_by_vector_csr 
    end type

    !> Diagonal storage matrix
    type, public :: DiagonalMatrix
        real(d), allocatable    :: data(:,:)        !! data[:,k] stores the data for diagonal with offset[k]
        integer, allocatable    :: offsets(:)
        integer                 :: ncol             !! Number of columns in the matrix
        integer                 :: nrow             !! Number of rows in the matrix
        integer                 :: ndiag            !! The length of the main diagonal
        integer                 :: idiag            !! The number of diagonals stored
      contains
        procedure               :: to_dense => to_dense_dia         !! Convert diagonal to dense matrix
        procedure               :: to_array => to_dense_dia         !! Alias of to_dense, for SciPy-like interface
        procedure               :: multiply => multiply_by_vector_dia
    end type

    interface CSRMatrix
        procedure :: init_csr
        procedure :: init_csr_from_dense 
    end interface

    interface DiagonalMatrix
        procedure :: init_dia
        procedure :: init_dia_from_dense 
    end interface

  contains

!! -------------!
!! CSR MATRICES !
!!--------------!

    !> Create a CSR matrix from data, column index and row pointers,
    !! which is how the matrix is stored internally.
    !! See https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_row_(CSR,_CRS_or_Yale_format)
    function init_csr(data, col_ind, row_ptr, ncol) result(self)
        type(CSRMatrix)     :: self             !! This CSR matrix
        real(d)             :: data(:)          !! The non-zero elements of the matrix
        integer             :: col_ind(:)       !! The indices of the columns of the values in the data array
        integer             :: row_ptr(:)       !! The index in the data array and col_ind where a given row starts
        integer, optional   :: ncol             !! The number of columns in the matrix. If not present, we assume
                                                !! it to be the maximum value in the col_ind array (i.e. there are no
                                                !! all-zero columns after the final column in the col_ind array)
        allocate(self%data, source=data)
        allocate(self%col_ind, source=col_ind)
        allocate(self%row_ptr, source=row_ptr)
        self%nrow = size(self%row_ptr) - 1 
        if (present(ncol)) then
            self%ncol = ncol
        else
            self%ncol = maxval(self%col_ind)
        end if
    end function

    !> Create a CSR matrix from a dense matrix
    function init_csr_from_dense(dense) result(self)
        type(CSRMatrix)     :: self             !! This CSR matrix
        real(d)             :: dense(:,:)       !! The dense matrix to turn into the CSR matrix
        integer             :: nnz, ierr        ! Number of non-zero elements, error integer
        ! Store the parameters for the dnscsr function and allocate space
        self%nrow = size(dense, dim=1)
        self%ncol = size(dense, dim=2)
        nnz = count(dense /= 0.0_d)
        allocate(self%data(nnz))
        allocate(self%col_ind(nnz))
        allocate(self%row_ptr(self%nrow+1))
        ! Call the function to create the data, col_ind and row_ptr arrays
        call dnscsr(self%nrow, self%ncol, nnz, dense, self%nrow, &
                    self%data, self%col_ind, self%row_ptr, ierr)
    end function

    !> Convert this CSR matrix to a dense matrix
    function to_dense_csr(self) result(dense)
        class(CSRMatrix)        :: self         !! This CSR matrix
        real(d), allocatable    :: dense(:,:)   !! The output dense matrix
        integer                 :: ierr         ! Error integer
        ! Allocate space for the dense matrix
        allocate(dense(self%nrow,self%ncol))
        ! Convert the CSR matrix into dense matrix format
        call csrdns(self%nrow, self%ncol, self%data, self%col_ind, self%row_ptr, &
                    dense, self%nrow, ierr)
    end function

    !> Multiply this CSR matrix by a vector
    function multiply_by_vector_csr(self, vector) result(product)
        class(CSRMatrix)        :: self                 !! This CSR matrix  
        real(d)                 :: vector(self%ncol)    !! The vector to multiple the matrix by
        real(d)                 :: product(self%nrow)   !! The output vector
        ! Perform the multiplication
        call amux(self%nrow, vector, product, self%data, self%col_ind, self%row_ptr)
    end function

!! ------------------!
!! DIAGONAL MATRICES !
!!-------------------!

    function init_dia(data, offsets, shape) result(self)
        type(DiagonalMatrix)    :: self
        real(d)                 :: data(:,:)
        integer                 :: offsets(:)
        integer                 :: shape(2)
        allocate(self%data, source=data)
        allocate(self%offsets, source=offsets)
        self%nrow = shape(1)
        self%ncol = shape(2)
        self%ndiag = size(self%data, dim=1)
        self%idiag = size(self%offsets)
    end function

    function init_dia_from_dense(dense) result(self)
        type(DiagonalMatrix)    :: self             !! This diagonal matrix
        real(d)                 :: dense(:,:)       !! The dense matrix to convert to diagonal storage
        type(CSRMatrix)         :: csr_matrix       !! Temporary CSR matrix
        integer                 :: n                !! Order of the matrix
        ! First, convert to CSR
        csr_matrix = CSRMatrix(dense)
        n = minval(shape(dense))
        ! Store the shape
        self%nrow = size(dense, dim=1)
        self%ncol = size(dense, dim=2)
        ! Now convert to diagonal
        self%data = csrdia_alt(n, csr_matrix%data, csr_matrix%col_ind, csr_matrix%row_ptr, self%offsets)
        self%ndiag = size(self%data, dim=1)
        self%idiag = size(self%offsets)
    end function

    !> Convert from diagonal storage to dense matrix
    function to_dense_dia(self) result(dense)
        class(DiagonalMatrix)   :: self
        real(d), allocatable    :: dense(:,:)
        integer                 :: i, j, imin, imax
        ! Allocate the space for the dense matrix
        allocate(dense(self%nrow, self%ncol))
        dense = 0.0_d
        ! Loop over the diagonals and the elements in them and 
        ! fill the dense matrix
        do j = 1, size(self%data, dim=2)
            ! Figure out how to slice this diagonal's data in the self%data array
            ! (because offset diagonals won't use all the elements in the array)
            ! -ve offsets: data is stored at the end of the self%data array
            ! +ve offsets: data is stored at the start of the self%data array
            imin = 1 - min(self%offsets(j), 0)
            imax = size(self%data, dim=1) - max(self%offsets(j), 0)
            do i = imin, imax
                dense(i,i+self%offsets(j)) = self%data(i,j)
            end do
        end do
    end function

    !> Multiply this diagonal stored matrix by a vector
    function multiply_by_vector_dia(self, vector) result(product)
        class(DiagonalMatrix)   :: self                 !! This matrix
        real(d)                 :: vector(self%ncol)    !! The vector to multiple the matrix by
        real(d)                 :: product(self%nrow)   !! The output vector
        ! The amuxd function does the multiplication, returning the result in product
        call amuxd(self%ndiag, vector, product, self%data, self%ndiag, self%idiag, self%offsets)
    end function

end module
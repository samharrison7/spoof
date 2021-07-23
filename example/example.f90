program main
    use Spoof
    implicit none

    type(CSRMatrix)         :: csr_matrix
    type(DiagonalMatrix)    :: dia_matrix
    double precision        :: dense(10,10)
    double precision        :: vec(10)

    dense = 42.0D+00
    vec = 10.0D+00

    ! Create the matrices
    csr_matrix = CSRMatrix(dense)
    dia_matrix = DiagonalMatrix(dense)

    ! Multiply them by vec
    vec = csr_matrix%multiply(vec)              ! vec is filled with 4200
    vec = dia_matrix%multiply(vec)              ! vec is filled with 176400

end program


module sparskit

contains
!> These sparse matrix utilities have been copied from the SPARSKIT library, under
!! the GNU Lesser General Public License. Slight modifications have been made
!! here to ease integration into the NanoFASE model.
!!
!! Original author: Yousef Saad, University of Minnesota
!! URL: https://www-users.cse.umn.edu/~saad/software/SPARSKIT/
!!      and https://people.math.sc.edu/Burkardt/f_src/sparsekit/sparsekit.f90
!!
subroutine dnscsr ( nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr )
!*****************************************************************************80
!
!! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine does not check whether an element is small.  It considers 
!    that A(I,J) is zero only if it is exactly equal to zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NZMAX, the maximum number of nonzero elements 
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0 means normal return;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
    implicit none

    integer ( kind = 4 ) ncol
    integer ( kind = 4 ) ndns
    integer ( kind = 4 ) nrow

    real ( kind = 8 ) a(*)
    real ( kind = 8 ) dns(ndns,ncol)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(nrow+1)
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(*)
    integer ( kind = 4 ) next
    integer ( kind = 4 ) nzmax

    ierr = 0
    next = 1
    ia(1) = 1

    do i = 1, nrow

    do j = 1, ncol

    if ( dns(i,j) /= 0.0D+00 ) then

    if ( nzmax < next ) then
        ierr = i
        return
    end if

    ja(next) = j
    a(next) = dns(i,j)
    next = next + 1

    end if

    end do

    ia(i+1) = next

    end do

    return
end

subroutine csrdns ( nrow, ncol, a, ja, ia, dns, ndns, ierr )
!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer ( kind = 4 ) NDNS, the dimension of the DNS array.
!
!    Output, integer ( kind = 4 ) IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
    implicit none

    integer ( kind = 4 ) ncol
    integer ( kind = 4 ) ndns

    real ( kind = 8 ) a(*)
    real ( kind = 8 ) dns(ndns,ncol)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(*)
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(*)
    integer ( kind = 4 ) k
    integer ( kind = 4 ) nrow
    
    ierr = 0
    dns(1:nrow,1:ncol) = 0.0D+00

    do i = 1, nrow
    do k = ia(i), ia(i+1)-1
        j = ja(k)
        if ( ncol < j ) then
        ierr = i
        return
        end if
        dns(i,j) = a(k)
    end do
    end do

    return
end

subroutine amux ( n, x, y, a, ja, ia )
!*****************************************************************************80
!
!! AMUX multiplies a CSR matrix A times a vector.
!
!  Discussion:
!
!    This routine multiplies a matrix by a vector using the dot product form.
!    Matrix A is stored in compressed sparse row storage.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
!
!    Input, real X(*), and array of length equal to the column dimension 
!    of A.
!
!    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real Y(N), the product A * X.
!
    implicit none

    integer ( kind = 4 ) n

    real ( kind = 8 ) a(*)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(*)
    integer ( kind = 4 ) ja(*)
    integer ( kind = 4 ) k
    real ( kind = 8 ) t
    real ( kind = 8 ) x(*)
    real ( kind = 8 ) y(n)

    do i = 1, n
!
!  Compute the inner product of row I with vector X.
!
    t = 0.0D+00
    do k = ia(i), ia(i+1)-1
        t = t + a(k) * x(ja(k))
    end do

    y(i) = t

    end do

    return
end
subroutine csrdia ( n, idiag, job, a, ja, ia, ndiag, diag, ioff, ao, &
    jao, iao, ind )
  
  !*****************************************************************************80
  !
  !! CSRDIA converts Compressed Sparse Row to diagonal format.
  !
  !  Discussion:
  !
  !    This routine extracts IDIAG diagonals from the input matrix A,
  !    JA, IA, and puts the rest of the matrix in the output matrix AO,
  !    JAO, IAO.  The diagonals to be extracted depend on the value of JOB.
  !
  !    In  the first case, the diagonals to be
  !    extracted are simply identified by their offsets provided in ioff
  !    by the caller.  In the second case, the code internally determines
  !    the idiag most significant diagonals, i.e., those diagonals of the
  !    matrix which have the largest number of nonzero elements, and
  !    extracts them.
  !
  !    The algorithm is in place: ao, jao, iao can be overwritten on
  !    a, ja, ia if desired.
  !
  !    When the code is required to select the diagonals (job >= 10)
  !    the selection of the diagonals is done from left to right
  !    as a result if several diagonals have the same weight (number
  !    of nonzero elemnts) the leftmost one is selected first.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, integer ( kind = 4 ) IDIAG.  On intput, the number of diagonals 
  !    to be extracted.  On output, IDIAG may be modified to the
  !    actual number of diagonals found.
  !
  !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
  !    Compressed Sparse Row format.
  !
  ! job      = integer ( kind = 4 ). serves as a job indicator.  Job is better thought
  !         of as a two-digit number job=xy. If the first (x) digit
  !         is one on entry then the diagonals to be extracted are
  !         internally determined. In this case csrdia exctracts the
  !         idiag most important diagonals, i.e. those having the largest
  !         number on nonzero elements. If the first digit is zero
  !         then csrdia assumes that ioff(*) contains the offsets
  !         of the diagonals to be extracted. there is no verification
  !         that ioff(*) contains valid entries.
  !         The second (y) digit of job determines whether or not
  !         the remainder of the matrix is to be written on ao,jao,iao.
  !         If it is zero  then ao, jao, iao is not filled, i.e.,
  !         the diagonals are found  and put in array diag and the rest is
  !         is discarded. if it is one, ao, jao, iao contains matrix
  !         of the remaining elements.
  !         Thus:
  !         job= 0 means do not select diagonals internally (pick those
  !                defined by ioff) and do not fill ao,jao,iao
  !         job= 1 means do not select diagonals internally
  !                      and fill ao,jao,iao
  !         job=10 means  select diagonals internally
  !                      and do not fill ao,jao,iao
  !         job=11 means select diagonals internally
  !                      and fill ao,jao,iao
  !
  !  Input, integer ( kind = 4 ) NDIAG, the first dimension of array DIAG.
  !
  ! on return:
  !
  ! diag  = real array of size (ndiag x idiag) containing the diagonals
  !         of A on return
  !
  ! ioff  = integer ( kind = 4 ) array of length idiag, containing the offsets 
  ! of the diagonals to be extracted.
  !
  ! ao, jao
  !  iao  = remainder of the matrix in a, ja, ia format.
  !
  ! work arrays:
  !
  ! ind   = integer ( kind = 4 ) array of length 2*n-1 used as work space.
  !         needed only when job>=10 i.e., in case the diagonals are to
  !         be selected internally.
  !
    implicit none
  
    integer ( kind = 4 ) idiag
    integer ( kind = 4 ) ndiag
  
    real ( kind = 8 ) a(*)
    real ( kind = 8 ) ao(*)
    real ( kind = 8 ) diag(ndiag,idiag)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ia(*)
    integer ( kind = 4 ) iao(*)
    integer ( kind = 4 ) idum
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) ind(*)
    integer ( kind = 4 ) ioff(*)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) ja(*)
    integer ( kind = 4 ) jao(*)
    integer ( kind = 4 ) jmax
    integer ( kind = 4 ) job
    integer ( kind = 4 ) job1
    integer ( kind = 4 ) job2
    integer ( kind = 4 ) k
    integer ( kind = 4 ) ko
    integer ( kind = 4 ) l
    integer ( kind = 4 ) n
    integer ( kind = 4 ) n2

  
    job1 = job / 10
    job2 = job - job1 * 10
  
    if ( job1 /= 0 ) then
  
      n2 = n + n - 1
      call infdia ( n, ja, ia, ind, idum )

      idiag = idum

  !
  !  Determine the diagonals to accept.
  !
    !   ii = 0
  
    !   do
  
    !     ii = ii + 1
    !     jmax = 0
  
    !     do k = 1, n2
  
    !       j = ind(k)
  
    !       if ( jmax < j ) then
    !         i = k
    !         jmax = j
    !       end if
  
    !     end do
  
    !     if ( jmax <= 0 ) then
    !       ii = ii - 1
    !       exit
    !     end if
  
    !     ioff(ii) = i - n
    !     ind(i) = - jmax
  
    !     if ( idiag <= ii ) then
    !       exit
    !     end if
  
    !   end do
    !   idiag = ii
  
    end if
  !
  !  Initialize DIAG to zero.
  !
    diag(1:n,1:idiag) = 0.0D+00
  
    ko = 1
  !
  !  Extract diagonals and accumulate remaining matrix.
  !
    do i = 1, n
  
       do k = ia(i), ia(i+1)-1
  
          j = ja(k)
  
          do l = 1, idiag
             if ( j - i == ioff(l) ) then
               diag(i,l) = a(k)
               go to 51
             end if
          end do
  !
  !  Append element not in any diagonal to AO, JAO, IAO.
  !
          if ( job2 /= 0 ) then
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko + 1
          end if
  
  51      continue
  
       end do
  
       if ( job2 /= 0 ) then
         ind(i+1) = ko
       end if
  
    end do
  !
  !  Finish with IAO.
  !
    if ( job2 /= 0 ) then
      iao(1) = 1
      iao(2:n+1) = ind(2:n+1)
    end if
  
    return
  end
  subroutine infdia ( n, ja, ia, ind, idiags )

    !*****************************************************************************80
    !
    !! INFDIA obtains information on the diagonals of A.
    !
    !  Discussion:
    !
    !    This routine finds the lengths of each of the 2*N-1 diagonals of A
    !
    !    It also outputs the number of nonzero diagonals found.
    !
    !  Modified:
    !
    !    07 January 2004
    !
    !  Author:
    !
    !    Youcef Saad
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the matrix.
    !
    !    Input, integer ( kind = 4 ) JA(*), IA(N+1), the matrix information (but
    !    no values) in CSR Compressed Sparse Row format.
    !
    !    Output, integer ( kind = 4 ) IND(2*N-1); The K-th entry in IND contains the number 
    !    of nonzero elements in diagonal K, the numbering being from the 
    !    lowermost diagonal (bottom-left).  In other words IND(K) = length 
    !    of diagonal whose offset with respect to the main diagonal is = - N + K.
    !
    !    Output, integer ( kind = 4 ) IDIAG, the number of nonzero diagonals found.
    !
      implicit none
    
      integer ( kind = 4 ) n
    
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ia(n+1)
      integer ( kind = 4 ), intent(out) :: idiags
      integer ( kind = 4 ) ind(*)
      integer ( kind = 4 ) j
      integer ( kind = 4 ) ja(*)
      integer ( kind = 4 ) k
      integer ( kind = 4 ) n2
    
      n2 = n+n-1
      ind(1:n2) = 0
    
      do i = 1, n
         do k = ia(i), ia(i+1)-1
            j = ja(k)
            ind(n+j-i) = ind(n+j-i) + 1
         end do
      end do
    !
    !  Count the nonzero ones.
    !
      idiags = 0
    
      do k = 1, n2
        if ( ind(k) /= 0 ) then
          idiags = idiags + 1
        end if
      end do

    
      return
    end

function csrdia_alt(n, a, ja, ia, ioff) result(diag)
    ! A simpler alternative to csrdia, which assumes we want to return *all* diagonals
    ! with non-zero elements
    implicit none
    integer(kind=4) n                         ! Order of the matrix
    real (kind=8) a(*)                        ! Input CSR matrix data
    integer ( kind = 4 ) ja(*)                ! Input CSR matrix column indices
    integer ( kind = 4 ) ia(*)                ! Input CSR matrix row pointers
    real(kind=8), allocatable :: diag(:,:)    ! Output diagonal matrix data
    integer ( kind = 4 ), allocatable :: ioff(:)           ! Output diagonal matrix offsets
    integer ( kind = 4 ) idiag                ! Number of diagonals with non-zero elements
    integer ( kind = 4 ) ndiag                ! Max number of elements per diagonal
    integer ( kind = 4 ) n2                   ! The maximum number of diagonals
    integer ( kind = 4 ), allocatable :: ind(:)   ! The kth entry in ind contains the number of non-zero elements in diagonal,
                                                  ! with numbering from the lower most (bottom left)
    integer ( kind = 4 ) i, k, j, l           ! Iterators
    logical :: found_diag = .false.           ! Flag to exit loop

    ! Store ndiag (max number of elements per diagonal) and idiag (number of diagonals)
    ndiag = n         ! Assume this is the order of the matrix
    n2 = n + n - 1    ! Maximum number of diagonals to loop over
    allocate(ind(n2))
    ind(1:n2) = 0
    ! Loop over the diagonals and calculate the number of non-zero elements in each
    do i = 1, n
      do k = ia(i), ia(i+1) - 1
        j = ja(k)
        ind(n+j-i) = ind(n+j-i) + 1
      end do
    end do

    ! Use this to calculate the number of diagonals with non-zero elements (i.e. ind(k) > 0)
    idiag = count(ind /= 0)
    if (allocated(ioff)) deallocate(ioff)
    allocate(ioff(idiag))
    i = 0
    do k = 1, n2
      if (ind(k) /= 0) then
        i = i + 1
        ioff(i) = k - n           ! Set this offset around the main diagonal, k is ordered from bottom left
      end if
    end do

    ! Now we have ndiag and idiag, we can allocate space to diag and ioff
    allocate(diag(ndiag, idiag))
    diag = 0.0D+00

    ! Extract diagonals and accumulate remaining matrix
    do i = 1, n
      do k = ia(i), ia(i+1)-1
        j = ja(k)
        l = 1
        do while ((.not. found_diag) .and. (l <= idiag))
          if ( j - i == ioff(l) ) then
            diag(i,l) = a(k)
            found_diag = .true.
          end if
          l = l + 1
        end do
        found_diag = .false.
      end do
   end do
end function

subroutine amuxd ( n, x, y, diag, ndiag, idiag, ioff )
  !*****************************************************************************80
  !
  !! AMUXD multiplies a DIA matrix times a vector.
  !
  !  Discussion:
  !
  !    This routine multiplies a matrix by a vector when the original matrix 
  !    is stored in the DIA diagonal storage format.
  !
  !  Modified:
  !
  !    07 January 2004
  !
  !  Author:
  !
  !    Youcef Saad
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
  !
  !    Input, real X(*), array of length equal to the column dimension of
  !    the A matrix.
  !
  !    Output, real Y(N), the product A * X.
  !
  !    Input, real DIAG(NDIAG,IDIAG), the diagonals.
  !
  !    Input, integer ( kind = 4 ) NDIAG, the first dimension of array adiag as 
  !    declared in the calling program.
  !
  !    Input, integer ( kind = 4 ) IDIAG, the number of diagonals in the matrix.
  !
  !    Input, integer ( kind = 4 ) IOFF(IDIAG), the offsets of the diagonals of 
  !    the matrix: diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
  !
    implicit none
  
    integer ( kind = 4 ) idiag
    integer ( kind = 4 ) n
    integer ( kind = 4 ) ndiag
  
    real ( kind = 8 ) diag(ndiag,idiag)
    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) io
    integer ( kind = 4 ) ioff(idiag)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    real ( kind = 8 ) x(n)
    real ( kind = 8 ) y(n)
  
    y(1:n) = 0.0D+00
  
    do j = 1, idiag
      io = ioff(j)
      i1 = max ( 1, 1 - io )
      i2 = min ( n, n - io )
      do k = i1, i2
        y(k) = y(k) + diag(k,j) * x(k+io)
      end do
    end do
  
    return
  end

end module
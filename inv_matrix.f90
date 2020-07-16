! 逆行列をLAPACKで計算
! http://sak12.blogspot.com/2012/09/fortran.html
! 実行:
! gfortran -o test inv_matrix.f90 -llapack -lblas
! ./test
! NIFS Plasma Simulator
! nfort -o exec_inv_matrix inv_matrix.f90 -llapack -lblas_sequential && qsub script_inv_matrix.sh

! --------------------------------
! calculate inverse of matrix A(m,n)
! --------------------------------
function inv_lapack(A) result (AI)
    implicit none

    real(8), intent(in) :: A(:,:)
    real(8) :: AI(size(A,1),size(A,2))
    integer, parameter :: NB = 64
    integer :: info
    integer :: m, n

    ! for lapack
    integer :: LWORK
    integer, allocatable :: IPIV(:) ! dimension N
    real(8), allocatable :: WORK(:) ! dimension LWORK

    ! get size of A(m,n)
    m = size(A,1)
    n = size(A,2)

    LWORK = n*NB
    allocate(IPIV(n))
    allocate(WORK(LWORK))

    AI = A
    call DGETRF(m, n, AI, m, IPIV(1:min(m,n)), info)  ! factorize
    call DGETRI(n, AI, m, IPIV, WORK, LWORK, info)  ! inverse

    deallocate(IPIV)
    deallocate(WORK)

    if (INFO==0) then
        !write(*,*) "succeeded."
    endif

    return
end function

! --------------------------------
! show matrix A(m,n)
! --------------------------------
subroutine show_dmat(A,fmtf)
    implicit none
    real(8), intent(in) :: A(:,:)
    character(len=*), intent(in), optional :: fmtf
    character(len=25) :: fmtc
    integer :: m, n, i, j
    m = size(A,1); n = size(A,2)

    if (present(fmtf)) then
        fmtc = fmtf
    else
        fmtc = '(f13.4," ")'
    endif
    do i=1,m
        do j=1,n
            write(*,fmtc, advance='no') A(i,j)
        enddo
        write(*,*) ""
    enddo
    return
end subroutine show_dmat


program test
    implicit none

    interface
        function inv_lapack(A)
            real(8), intent(in) :: A(:,:)
            real(8) :: inv_lapack(size(A,2),size(A,1))
        end function inv_lapack

        subroutine show_dmat(A,fmtf)
            real(8), intent(in) :: A(:,:)
            character(len=*), intent(in), optional :: fmtf
        end subroutine show_dmat
    end interface

    real(8) :: A(2,2)
    real(8) :: AI(2,2)
    integer :: i, j

    A(1,1) = 1.d0
    A(1,2) = 2.d0
    A(2,1) = 4.5d0
    A(2,2) = 6.0d0

    write(*,'("A=")')
    call show_dmat(A, "(f9.3)")
    AI = inv_lapack(A)
    write(*,'("AI=")')
    call show_dmat(AI, "(f9.3)")

    stop
end program test

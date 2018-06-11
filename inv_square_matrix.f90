! 正方行列の逆行列をLAPACKで計算
! http://sak12.blogspot.com/2012/09/fortran.html
! 実行:
! gfortran -o test inv_matrix.f90 -llapack -lblas
! ./test


function inv_lapack(mat) result(mat_I)
  implicit none
  double precision, intent(in) :: mat(:,:)
  double precision :: mat_I(size(mat,1),size(mat,2))
  integer, parameter :: NB = 64 !このパラメータの意味?
  integer :: info
  integer :: n
  !For LAPACK
  integer :: LWORK
  integer, allocatable :: IPIV(:)
  double precision, allocatable :: WORK(:)

  ! 行列のサイズを取得
  n = size(mat,1)
  if ( n.ne.size(mat,2) ) then  !正方行列の判定
    stop 'Matrix is not square'
  end if

  LWORK = n*NB
  allocate(IPIV(n))
  allocate(WORK(LWORK))

  mat_I = mat
  call DGETRF(n, n, mat_I, n, IPIV(1:n), info) ! 対角化
  call DGETRI(n, mat_I, n, IPIV, WORK, LWORK, info) ! 逆行列

  deallocate(IPIV)
  deallocate(WORK)

  return

end function

subroutine show_mat(mat)
    implicit none
    double precision, intent(in) :: mat(:,:)
    integer :: n, i, j
    n = size(mat,1)

    do i = 1, n
      do j = 1, n
        write(*,'(f9.3)', advance='no') mat(i,j)
      end do
      write(*,*) ''
    end do

end subroutine show_mat

program main
  implicit none

  interface
      function inv_lapack(mat)
          double precision, intent(in) :: mat(:,:)
          double precision :: inv_lapack(size(mat,2),size(mat,1))
      end function inv_lapack

      subroutine show_mat(A)
          double precision, intent(in) :: A(:,:)
      end subroutine show_mat
  end interface

  double precision :: mat(2,2)
  double precision :: mat_I(2,2)
  integer i, j

  mat(1,1) = 1.d0
  mat(1,2) = 2.d0
  mat(2,1) = 4.5d0
  mat(2,2) = 6.0d0

  call show_mat(mat)
  mat_I = inv_lapack(mat)
  call show_mat(mat_I)


end program main

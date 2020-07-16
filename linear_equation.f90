!対称行列を係数行列とする連立一次方程式を解く
! 実行:
! gfortran -o test linear_equation.f90 -llapack -lblas
! ./test
! NIFS Plasma Simulator
! nfort -o exec_linear_equation linear_equation.f90 -llapack -lblas_sequential && qsub script_linear_equation.sh

subroutine set_A(mat)
  double precision, intent(out) :: mat(:,:)
  mat = reshape( (/1.80d0, 5.25d0, 1.58d0,-1.11d0, &
                 2.88d0,-2.95d0,-2.69d0,-0.66d0, &
                 2.05d0,-0.95d0,-2.90d0,-0.59d0, &
                -0.89d0,-3.80d0,-1.04d0, 0.80d0/), &
              (/4,4/) )
end subroutine set_A

subroutine set_B(mat)
  double precision, intent(out) :: mat(:)
  mat = (/ 9.52d0, 24.35d0, 0.77d0, -6.22d0 /)
end subroutine set_B

function dgesv_lapack(A,B) result(X)
  implicit none
  double precision, intent(in) :: A(:,:)
  double precision, intent(in) :: B(:)
  double precision :: X(size(A,1))
  character(len=1), parameter :: UPLO = 'U' ! 上三角要素を使用
  integer :: N
  integer, allocatable :: IPIV(:)
  double precision :: LW
  double precision, allocatable :: WORK(:)
  integer :: LWORK
  integer :: INFO

  ! 係数/右辺行列のサイズ
  N = size(A,1)

  allocate(IPIV(N))

  !計算
  call dgesv(N, 1, A, N, IPIV, B, N, INFO)

  X = B

  return
end function dgesv_lapack

subroutine show_mat(mat)
    implicit none
    double precision, intent(in) :: mat(:,:)
    integer :: m, n, i, j
    m = size(mat,1)
    n = size(mat,2)

    do i = 1, m
      do j = 1, n
        write(*,'(f9.3)', advance='no') mat(i,j)
      end do
      write(*,*) ''
    end do

end subroutine show_mat

subroutine show_vec(vec)
    implicit none
    double precision, intent(in) :: vec(:)
    integer :: m, i
    m = size(vec,1)

    do i = 1, m
        write(*,'(f9.3)') vec(i)
      write(*,*) ''
    end do

end subroutine show_vec


program main
  implicit none
  interface
      function dgesv_lapack(A, B)
          double precision, intent(in) :: A(:,:)
          double precision, intent(in) :: B(:)
          double precision :: dgesv_lapack(size(A,1))
      end function dgesv_lapack

      subroutine set_A(A)
          double precision, intent(out) :: A(:,:)
      end subroutine set_A

      subroutine set_B(A)
          double precision, intent(out) :: A(:)
      end subroutine set_B

      subroutine show_mat(A)
          double precision, intent(in) :: A(:,:)
      end subroutine show_mat

      subroutine show_vec(A)
          double precision, intent(in) :: A(:)
      end subroutine show_vec
  end interface

  double precision, dimension(4,4) :: A
  double precision, dimension(4) :: B
  double precision, dimension(4) :: X

  ! 係数行列の定義
  call set_A(A)
  call show_mat(A)

  !右辺行列の定義
  call set_B(B)
  call show_vec(B)

  !行列方程式を解く
  X = dgesv_lapack(A,B)
  call show_vec(X)


end program main

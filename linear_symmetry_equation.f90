!対称行列を係数行列とする連立一次方程式を解く
! 実行:
! gfortran -o test linear_symmetry_equation.f90 -llapack -lblas
! ./test
! NIFS Plasma Simulator
! nfort -o exec_linear_symmetry_equation linear_symmetry_equation.f90 -llapack -lblas_sequential && qsub script_linear_symmetry_equation.sh

subroutine set_A(mat)
  double precision, intent(out) :: mat(:,:)
  mat(1,1:5) = (/-5.86d0, 3.99d0, -5.93d0,-2.82d0, 7.69d0/)
  mat(2,1:5) = (/3.99d0, 4.46d0, 2.58d0, 4.42d0, 4.61d0/)
  mat(3,1:5) = (/-5.93d0, 2.58d0, -8.52d0, 8.57d0, 7.69d0/)
  mat(4,1:5) = (/-2.82d0, 4.42d0, 8.57d0, 3.72d0, 8.07d0/)
  mat(5,1:5) = (/7.69d0, 4.61d0, 7.69d0, 8.07d0, 9.83d0/)
end subroutine set_A

subroutine set_B(mat)
  double precision, intent(out) :: mat(:,:)
  mat(1,1:3) = (/1.32d0, -6.33d0, -8.77d0/)
  mat(2,1:3) = (/2.22d0, 1.69d0, -8.33d0/)
  mat(3,1:3) = (/0.12d0, -1.56d0, 9.54d0/)
  mat(4,1:3) = (/-6.41d0, -9.49d0, 9.56d0/)
  mat(5,1:3) = (/6.33d0, -3.67d0, 7.48d0/)
end subroutine set_B

function dsysv_lapack(A,B) result(X)
  implicit none
  double precision, intent(in) :: A(:,:)
  double precision, intent(in) :: B(:,:)
  double precision :: X(size(A,1),size(B,2))
  character(len=1), parameter :: UPLO = 'U' ! 上三角要素を使用
  integer :: N, NRHS
  integer, allocatable :: IPIV(:)
  double precision :: LW
  double precision, allocatable :: WORK(:)
  integer :: LWORK
  integer :: INFO

  ! 係数/右辺行列のサイズ
  N = size(A,1)
  NRHS = size(B,2)

  allocate(IPIV(N))

  ! 最適なワークスペース用配列を確保
  call DSYSV(UPLO, N, NRHS, A, N, IPIV, B, N, LW, -1, INFO)
  LWORK = int(LW)
  allocate(WORK(1:LWORK))

  !計算
  call DSYSV(UPLO, N, NRHS, A, N, IPIV, B, N, WORK, LWORK, INFO)
  deallocate(WORK)

  X = B

  return
end function dsysv_lapack

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


program main
  implicit none
  interface
      function dsysv_lapack(A, B)
          double precision, intent(in) :: A(:,:)
          double precision, intent(in) :: B(:,:)
          double precision :: dsysv_lapack(size(A,1),size(B,2))
      end function dsysv_lapack

      subroutine set_A(A)
          double precision, intent(out) :: A(:,:)
      end subroutine set_A

      subroutine set_B(A)
          double precision, intent(out) :: A(:,:)
      end subroutine set_B

      subroutine show_mat(A)
          double precision, intent(in) :: A(:,:)
      end subroutine show_mat
  end interface

  double precision, dimension(5,5) :: A
  double precision, dimension(5,3) :: B
  double precision, dimension(5,3) :: X

  ! 係数行列の定義
  call set_A(A)
  call show_mat(A)

  !右辺行列の定義
  call set_B(B)
  call show_mat(B)

  !行列方程式を解く
  X = dsysv_lapack(A,B)
  call show_mat(X)


end program main

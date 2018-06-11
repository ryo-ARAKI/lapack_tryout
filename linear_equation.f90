!対称行列を係数行列とする連立一次方程式を解く

subroutine set_A(mat)
  double precision, intent(out) :: mat(:,:)
  mat(1,1) = -5.86d0
  mat(1,2) = 3.99d0
  mat(1,3) = -5.93d0
  mat(1,4) = -2.82d0
  mat(1,5) = 7.69d0
  mat(2,1) = 3.99d0
  mat(2,2) = 4.46d0
  mat(2,3) = 2.58d0
  mat(2,4) = 4.42d0
  mat(2,5) = 4.61d0
  mat(3,1) = -5.93d0
  mat(3,2) = 2.58d0
  mat(3,3) = -8.52d0
  mat(3,4) = 8.57d0
  mat(3,5) = 7.69d0
  mat(4,1) = -2.82d0
  mat(4,2) = 4.42d0
  mat(4,3) = 8.57d0
  mat(4,4) = 3.72d0
  mat(4,5) = 8.07d0
  mat(5,1) = 7.69d0
  mat(5,2) = 4.61d0
  mat(5,3) = 7.69d0
  mat(5,4) = 8.07d0
  mat(5,5) = 9.83d0
end subroutine set_A

subroutine set_B(mat)
  double precision, intent(out) :: mat(:,:)
  mat(1,1) = 1.32d0
  mat(1,2) = -6.33d0
  mat(1,3) = -8.77d0
  mat(2,1) = 2.22d0
  mat(2,2) = 1.69d0
  mat(2,3) = -8.33d0
  mat(3,1) = 0.12d0
  mat(3,2) = -1.56d0
  mat(3,3) = 9.54d0
  mat(4,1) = -6.41d0
  mat(4,2) = -9.49d0
  mat(4,3) = 9.56d0
  mat(5,1) = 6.33d0
  mat(5,2) = -3.67d0
  mat(5,3) = 7.48d0
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
  end interface

  double precision, dimension(5,5) :: A
  double precision, dimension(5,3) :: B
  double precision, dimension(5,3) :: X
  integer m, n, o
  integer i, j

  m = size(A,1)
  n = size(A,2)
  o = size(B,2)

  ! 係数行列の定義
  call set_A(A)
  do i = 1, m
    do j = 1, n
      write(*,'(f9.3)', advance='no') A(i,j)
    end do
    write(*,*) ''
  end do

  !右辺行列の定義
  call set_B(B)
  do i = 1, n
    do j = 1, o
      write(*,'(f9.3)', advance='no') B(i,j)
    end do
    write(*,*) ''
  end do

  !行列方程式を解く
  X = dsysv_lapack(A,B)
  do i = 1, m
    do j = 1, o
      write(*,'(f9.3)', advance='no') X(i,j)
    end do
    write(*,*) ''
  end do


end program main

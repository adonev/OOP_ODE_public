module MatrixExponential
    ! Code adopted by Aleks Donev from 
    ! https://people.sc.fsu.edu/~jburkardt/f_src/matrix_exponential/matrix_exponential.html
    use Precision
    implicit none
    private
    
    integer, parameter :: rk = wp ! Working precision
    public :: mat_expm_matlab, mat_expm_taylor

contains    

subroutine mat_expm_matlab ( n, a, e )

!*****************************************************************************80
!
!! R8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Input:
!
!    integer N, the dimension of the matrix.
!
!    real ( kind = rk ) A(N,N), the matrix.
!
!  Output:
!
!    real ( kind = rk ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) a_norm
  real ( kind = rk ) c
  real ( kind = rk ) d(n,n)
  real ( kind = rk ) e(n,n)
  integer ee
  integer k
  logical p
  integer , parameter :: q = 6
  real ( kind = rk ) r8_log_2
  real ( kind = rk ) r8mat_norm_li
  integer s
  real ( kind = rk ) x(n,n)

  a2(1:n,1:n) = a(1:n,1:n)

  a_norm = r8mat_norm_li ( n, n, a2 )

  ee = int ( r8_log_2 ( a_norm ) ) + 1

  s = max ( 0, ee + 1 )

  a2(1:n,1:n) = a2(1:n,1:n) / 2.0D+00 ** s

  x(1:n,1:n) = a2(1:n,1:n)

  c = 0.5D+00

  call r8mat_identity ( n, e )
  e(1:n,1:n) = e(1:n,1:n) + c * a2(1:n,1:n)

  call r8mat_identity ( n, d )
  d(1:n,1:n) = d(1:n,1:n) - c * a2(1:n,1:n)

  p = .true.

  do k = 2, q

    c = c * real ( q - k + 1, kind = rk ) &
      / real ( k * ( 2 * q - k + 1 ), kind = rk )

    x(1:n,1:n) = matmul ( a2(1:n,1:n), x(1:n,1:n) )

    e(1:n,1:n) = e(1:n,1:n) + c * x(1:n,1:n)

    if ( p ) then
      d(1:n,1:n) = d(1:n,1:n) + c * x(1:n,1:n)
    else
      d(1:n,1:n) = d(1:n,1:n) - c * x(1:n,1:n)
    end if

    p = .not. p

  end do
!
!  E -> inverse(D) * E
!
  call r8mat_minvm ( n, n, d, e, e )
!
!  E -> E^(2*S)
!
  do k = 1, s
    e(1:n,1:n) = matmul ( e(1:n,1:n), e(1:n,1:n) )
  end do

  return
end


subroutine mat_expm_taylor ( n, a, e )

!*****************************************************************************80
!
!! R8MAT_EXPM2 uses the Taylor series for the matrix exponential.
!
!  Discussion:
!
!    Formally,
!
!      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
!
!    This function sums the series until a tolerance is satisfied.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 November 2011
!
!  Author:
!
!    Cleve Moler, Charles Van Loan
!
!  Reference:
!
!    Cleve Moler, Charles VanLoan,
!    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
!    Twenty-Five Years Later,
!    SIAM Review,
!    Volume 45, Number 1, March 2003, pages 3-49.
!
!  Input:
!
!    integer N, the dimension of the matrix.
!
!    real ( kind = rk ) A(N,N), the matrix.
!
!  Output:
!
!    real ( kind = rk ) E(N,N), the estimate for exp(A).
!
  implicit none

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) e(n,n)
  real ( kind = rk ) f(n,n)
  integer k
  logical r8mat_is_insignificant

  e(1:n,1:n) = 0.0D+00

  call r8mat_identity ( n, f )

  k = 1

  do

    if ( r8mat_is_insignificant ( n, n, e, f ) ) then
      exit
    end if

    e(1:n,1:n) = e(1:n,1:n) + f(1:n,1:n)

    f(1:n,1:n) = matmul ( a(1:n,1:n), f(1:n,1:n) ) / real ( k, kind = rk )
    k = k + 1

  end do

  return
end

end module MatrixExponential


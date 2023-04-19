! Compile with
! gfortran -o TestODEIntegrator.x TestODEIntegrator.f90  MatrixExponential.o MatrixInverse.o Precision.o -llapack
! and then execute ./TestODEIntegrator.x
program TestODEIntegrator
    use Precision
    use LinearOperator
    use BandedMatrix
    use ODEIntegrand
    implicit none
                    
    call TestLinearIntegrand(n=2) ! Tests linear algebra in BandedMatrix class
    
contains

function A_lambda(n,lambda) result(A)
    integer, intent(in) :: n
    real(wp), intent(in) :: lambda
    real(wp) :: A(n,n) 
    
    A(1,1)=0.222222e-1_wp+0.888889_wp*lambda
    A(1,2)=0.628539e-1_wp-0.314269_wp*lambda
    A(2,1)=0.628536e-1_wp-0.314269_wp*lambda
    A(2,2)=0.177777_wp+0.111111_wp*lambda
    !write(*,*) "A=", A
    
end function

subroutine TestLinearIntegrand(n)
    integer, intent(in) :: n
 
    integer :: i
    real(wp), dimension(n,n) :: A, invA, A_alpha
    real(wp), dimension(n) :: x, b, b_tmp
    real(wp) :: lambda, alpha
 
    type(BandMat), target :: jacobian
    type(RHSFunc) :: rhs 
    
    lambda=1.0_wp           
    A = A_lambda(n,lambda)
    
    ! From Maple:
    invA(1,1) = 0.144444337752028318e1_wp
    invA(1,2) = 0.125707793940648971e1_wp
    invA(2,1) = 0.125707793940648971e1_wp
    invA(2,2) = 0.455556565221083165e1_wp
    
    jacobian = BandMat(n=n, m=n, n_bands=[-1,-1], id=1)
    call jacobian%create(n=n)
    call jacobian%init(matrix=A, tol=1e-6_wp, alpha=0.0_wp)

    rhs = RHSFunc(n=n, is_linear=.true., jacobian=jacobian, id=1)

    !--------------------
    ! Solvers
    !--------------------           
    call random_number(b)
    call rhs%solve(b=b, x=x)
    write(*,*) "Error in A^(-1)*b = ", x-matmul(invA,b)
    call rhs%eval(x=x, f=b_tmp)
    write(*,*) "Error in eval A*x =", b_tmp-matmul(A,x)
    write(*,*) "Error in residual A*x-b = ", b-b_tmp

    alpha=1.0_wp ! Try nonzero alpha
    A_alpha=A
    do i=1,n
       A_alpha(i,i)=A(i,i)+alpha
    end do     
    call rhs%solve(b=b, x=x, alpha=1.0_wp)
    call rhs%eval(x=x, f=b_tmp)
    write(*,*) "Error in eval (alpha*I+A)*x =", b_tmp-matmul(A_alpha,x)
    write(*,*) "Error in residual = (alpha*I+A)*x-b = ", b-b_tmp
    
end subroutine       

end program


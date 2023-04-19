! TODO: Make t a vector, dont assume A is invertible, if alpha=0 reuse factorization twice

module BandedMatrix
    use, intrinsic :: ISO_C_BINDING ! For interfacing with C
    use Precision
    use LinearOperator
    use MatrixInverse
    use TridiagonalSolver
    use MatrixExponential
    implicit none
    private

    ! The BandedMatrix class implements a LinearOperator interface for banded matrices
    ! In this example code only fully dense, tridiagonal, and diagonal matrices are supported,
    ! but full support for all banded matrices is easy to add using an interface to LAPACK
    
    ! Abstract representation of a linear map/matrix/operator A from R^n->R^m
    type, extends(LinOp), public :: BandMat
        integer  :: n_bands(2)=-1 ! Range of band below and above the diagonal, negative for fully dense
            ! [-1,-1] fully dense, [-1:0] lower triangular, [0:-1] upper triagonal, 
            ! [0,0] means diagonal, [1:1] means tridiagonal, [2:2] pentadiagonal, etc.
        real(wp), dimension(:,:), allocatable :: A, A_p_aI, invA, expA 
            ! Matrix A, A+alpha*Id, (A+alpha*Id)^(-1), and exp(A*t)
            ! In more serious codes invA would be instead an LU/Cholesky factorization or alike
        logical(c_bool) :: periodic=.false. ! Does this banded matrix come from periodic BCs in 1D    
        logical, private :: ready=.false. ! Must call constructor first
    contains
        ! Overloaded procedures from parent class LinOp
        procedure, pass(this) :: create => CreateBandMat
        procedure, pass(this) :: init => InitBandMat
        procedure, pass(this) :: destroy => DestroyBandMat
        procedure, pass(this) :: apply => ApplyBandMat
        procedure, pass(this) :: solve => SolveBandSys
        procedure, pass(this) :: exp_t => ApplyBandExp
    end type

    ! Add constructor of class with same name as type (convention)
    interface BandMat
        module procedure ConstructBandMat
    end interface
    
contains
    
function ConstructBandMat(n,m,n_bands,periodic,id) result(this) 
   ! This constructor is required because LinOp is an abstract type
   type(BandMat) :: this ! Constructed object
   integer, intent(in) :: n ! Must be supplied
   integer, intent(in), optional :: m ! Defaults to n (square matrix)
   integer, intent(in), optional :: n_bands(2) ! Number of bands, defaults to [-1,-1] 
   logical, intent(in), optional :: periodic          
   integer, intent(in), optional :: id
   
   if(present(m)) then
       if(.not.(m==n)) stop "Only square matrices (n=m) supported in CreateBandMat"
   end if
   
   ! call this%LinOp%create(n=n,m=n,id=id) ! This seems to be illegal in Fortran since type is abstract
   call CreateLinOp(this,n=n,m=n,id=id) ! Use the utility routine

   if(present(n_bands)) this%n_bands=n_bands 
   if(present(periodic)) this%periodic=periodic       
   this%ready=.true.
         
end function

subroutine CreateBandMat(this,n,m,id) ! Create (allocate memory, open files, etc.)
   class(BandMat), intent(inout) :: this
   integer, intent(in) :: n ! Must be supplied
   integer, intent(in), optional :: m ! Defaults to n (square matrix)            
   integer, intent(in), optional :: id
   
   integer :: i
    
   if(.not.this%ready) then
       stop "Banded matrix must be created using the constructor for BandMat or a type extension of BandMat"
   end if
   
   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       ! We only need to store the diagonal
       allocate(this%A(n,1), this%A_p_aI(n,1), this%invA(n,1), this%expA(n,1))
       this%A(:,1)=1.0_wp ! Default to identity matrix
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       allocate(this%A(n,3), this%A_p_aI(n,3)) ! Only need to store 3 diagonals
       ! Default to identity matrix:
       this%A(:,1)=0.0_wp
       this%A(:,2)=1.0_wp
       this%A(:,3)=0.0_wp
   else ! dense matrix by default
       allocate(this%A(n,n), this%A_p_aI(n,n), this%invA(n,n), this%expA(n,n))
       ! Default to identity matrix:
       this%A=0.0_wp
       do i=1,n
          this%A(i,i)=1.0_wp
       end do   
   end if

   write(*,*) "Created banded matrix with id=", this%id
  
end subroutine

subroutine InitBandMat(this,matrix,t,alpha,dt,tol) ! Initialize (set values, factorize, etc.)
   class(BandMat), intent(inout) :: this
   real(wp), intent(in), optional :: matrix(:,:) ! Data used to create the matrix   
   real(wp), intent(in), optional :: t ! If present, this is a parameter of the matrix, A<-A(t), e.g., time in an ODE
   real(wp), intent(in), optional :: dt ! If present, assume ApplyExpOp is called with the same dt, and precompute exp(A*dt)
   real(wp), intent(in), optional :: alpha ! Solve is (alpha*Id+A)*x=b; default is alpha=0
   real(wp), intent(in), optional :: tol ! Maximum relative error tolerance
   
   integer :: i, n   
   n=this%n

   call InitLinOp(this, matrix=matrix, t=t, alpha=alpha, dt=dt, tol=tol)

   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       if(.not.(size(matrix,2)==1)) write(*,*) &
          "WARNING: matrix should be of size (n,1) for diagonal matrices", size(matrix)
       if(present(matrix)) this%A(:,1) = matrix(:,1)
       this%A_p_aI=this%alpha+this%A
       this%invA(:,1)=1.0_wp/this%A_p_aI(:,1)
       if(present(dt)) this%expA(:,1)=exp(this%A(:,1)*dt)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       if(.not.(size(matrix,2)==3)) write(*,*) &
          "WARNING: matrix should be of size (n,3) for diagonal matrices", size(matrix)
       if(present(matrix)) this%A(:,1:3) = matrix(:,1:3)
       this%A_p_aI=this%A
       this%A_p_aI(:,2) = this%A(:,2) + this%alpha ! Add a multiple of identity to diagonal
       ! No need to precompute inverse for tridiagonal
       ! For now we don't have a special function for exponential of a tridiagonal matrix
   else ! dense matrix by default
       if(.false.) then ! Introduce a BUG on purpose for MSG debugging demo
           if(.not.(size(matrix,2)==n)) write(*,*) &
              "WARNING: matrix should be of size (n,n) for dense matrices", size(matrix) 
           if(present(matrix)) this%A = matrix
       else ! No BUG
           if(present(matrix)) then
              if(.not.(size(matrix,2)==n)) write(*,*) &
                 "WARNING: matrix should be of size (n,n) for dense matrices", size(matrix)
              this%A = matrix
           end if
       end if         
       this%A_p_aI=this%A
       do i=1,n
          this%A_p_aI(i,i)=this%A(i,i)+this%alpha
       end do   
       this%invA=mat_inv(this%A_p_aI)
       if(present(dt)) call mat_expm_matlab(n=n, a=this%A*this%dt, e=this%expA)
   end if
   
end subroutine

subroutine ApplyBandMat(this,x,b) ! Compute b=A*x
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(out) :: b(:) ! Actual argument: b(m)

   integer :: n   
   
   n=this%n

   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       b=this%A(:,1)*x(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       stop "tridiagonal matrix times vector not yet implemented"
   else ! dense matrix by default
       b=matmul(this%A,x) ! Use built in intrinsic function matmul
   end if
   
end subroutine

! Child classes are free to restrict the A^(-1) and exp(A*t) operations
! to n=m, or, to compute something like a least squares type solution

subroutine SolveBandSys(this,b,x,alpha) ! Solve (alpha*Id+A)*x=b to tolerance
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: b(:) ! Actual argument: b(m)
   real(wp), intent(out) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(in), optional :: alpha
   
   real(wp) :: alpha_
   integer :: n   
   
   n=this%n
   if(.not.(this%m==n)) stop "SolveBandSys only works for square matrices"
   
   alpha_=this%alpha
   if(present(alpha)) alpha_=alpha
   
   if(.not.(alpha_==this%alpha)) then
       write(*,*) "Recomputing inverse for BandedMatrix with id=", this%id
       call this%init(alpha=alpha_)
   end if            
      
   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       x=this%invA(:,1)*b(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       call solve_tridiag(a=this%A_p_aI(:,1), b=this%A_p_aI(:,2), c=this%A_p_aI(:,3), &
                          r=b, u=x, n=n, periodic=this%periodic)
   else ! dense matrix by default
       x=matmul(this%invA,b) ! x=(alpha*Id+A)^(-1)*b
   end if
   
end subroutine

subroutine ApplyBandExp(this,x,b,dt) ! Compute b=exp(A*dt)*x to tolerance
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(out) :: b(:) ! Actual argument: b(m)
   real(wp), intent(in), optional :: dt ! If not present, use the t set when initialized

   integer :: n   
   real(wp) :: dt_ ! We want exp(A*dt_)*x
   
   n=this%n
   if(.not.(this%m==n)) stop "ApplyBandExp only works for square matrices"

   dt_=this%dt
   if(present(dt)) dt_=dt 

   if(.not.(dt_==this%dt)) then
      write(*,*) "Recomputing expA for BandedMatrix with id=", this%id
      call this%init(dt=dt_)
   end if
   
   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       b=this%expA(:,1)*x(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       stop "Exponential of a tridiagonal matrix not yet implemented"
   else ! dense matrix by default
       b=matmul(this%expA,x) ! b=exp(A*t)*x               
   end if    
   
end subroutine
        
subroutine DestroyBandMat(this) ! Does nothing but print a message to confirm destruction
  class(BandMat), intent(inout) :: this

  write(*,*) "Destroying banded matrix with id=", this%id
  
  if(allocated(this%A)) deallocate(this%A)
  if(allocated(this%invA)) deallocate(this%invA)
  if(allocated(this%expA)) deallocate(this%expA)
  
  this%ready = .false.
  
  !call this%LinOp%destroy() ! This seems to illegal in Fortran, since LinOp is abstract
  call DestroyLinOp(this)

end subroutine

end module BandedMatrix


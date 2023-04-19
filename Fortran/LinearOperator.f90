! TODO: Make t a vector, dont assume A is invertible, if alpha=0 reuse factorization twice

module LinearOperator
    use Precision
    implicit none
    private    

    ! The LinearOperator class contains methods for computing b=A*x, solving the linear system Ax=b, computing exp(A*t)*x
    ! an initialize method that precomputes factorizations, and a finalize method that frees memory
    ! The matrix may optionally depend on some parameter, A<-A(t), e.g., time in an ODE
   
    public :: CreateLinOp, InitLinOp, DestroyLinOp ! base implementations to be used by parent classes if desired
    
    ! Abstract representation of a linear map/matrix/operator A from R^n->R^m
    type, abstract, public :: LinOp 
        integer  :: n=1, m=1 ! Op: R^n->R^m
        integer  :: id=0 ! An integer used to identify instances of operator
        real(wp) :: tol=-1.0_wp ! Maximum relative tolerance required for operations on A
                ! Here we will assume same accuracy is needed for all operations
        real(wp) :: dt=0.0_wp ! If ApplyExpOp is called with the same t, store it  
        real(wp) :: t=0.0_wp ! An internal parameter of A(t)
        real(wp) :: alpha=0.0_wp ! Solve is (alpha*Id+A)*x=b
    contains
        procedure(CreateLinOp), deferred, pass(this) :: create
        procedure(InitLinOp), deferred, pass(this) :: init
        procedure(DestroyLinOp), deferred, pass(this) :: destroy
        procedure(ApplyLinOp), deferred, pass(this) :: apply
        procedure(SolveLinSys), deferred, pass(this) :: solve
        procedure(ApplyExpOp), deferred, pass(this) :: exp_t
    end type

    abstract interface
        
        subroutine ApplyLinOp(this,x,b) ! Compute b=A*x
            import ! Import all variables from this module
            class(LinOp), intent(inout) :: this
            real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
            real(wp), intent(out) :: b(:) ! Actual argument: b(m)
        end subroutine

        ! Child classes are free to restrict the A^(-1) and exp(A*t) operations
        ! to n=m, or, to compute something like a least squares type solution

        subroutine SolveLinSys(this,b,x,alpha) ! Solve b=(alpha*Id+A)*x to tolerance
            import ! Import all variables from this module
            class(LinOp), intent(inout) :: this
            real(wp), intent(in) :: b(:) ! Actual argument: b(m)
            real(wp), intent(out) :: x(:)  ! Actual argument: x(n)
            real(wp), intent(in), optional :: alpha
        end subroutine

        subroutine ApplyExpOp(this,x,b,dt) ! Compute b=exp(A*dt)*x to tolerance
            import ! Import all variables from this module
            class(LinOp), intent(inout) :: this
            real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
            real(wp), intent(out) :: b(:) ! Actual argument: b(m)
            real(wp), intent(in), optional :: dt ! If not present, use the dt given during initialization
        end subroutine
        
    end interface    

contains

    subroutine CreateLinOp(this,n,m,id) ! Create (allocate memory, open files, etc.)
        class(LinOp), intent(inout) :: this
        integer, intent(in) :: n ! Must be supplied
        integer, intent(in), optional :: m ! Defaults to n (square matrix)            
        integer, intent(in), optional :: id  

        if(present(id)) this%id=id
        write(*,*) "Creating linear operator with id=", this%id

        this%n=n
        if(present(m)) then
            this%m=m
        else
            this%m=n
        end if
        
    end subroutine

    subroutine InitLinOp(this,matrix,t,alpha,dt,tol) ! Initialize/update (set values, factorize, etc.)
        class(LinOp), intent(inout) :: this
        real(wp), intent(in), optional :: matrix(:,:) ! Data used to create the matrix
        real(wp), intent(in), optional :: t ! If present, this is a parameter of the matrix, A<-A(t), e.g., time in an ODE
        real(wp), intent(in), optional :: alpha ! Solve is (alpha*Id+A)*x=b; default is alpha=0
        real(wp), intent(in), optional :: dt ! If present, assume ApplyExpOp is called with the same dt, and precompute exp(A*dt)
        real(wp), intent(in), optional :: tol ! Maximum relative error tolerance

        if(present(t)) this%t=t
        if(present(alpha)) this%alpha=alpha
        if(present(dt)) this%dt=dt
        if(present(tol)) this%tol=tol

    end subroutine
        
    subroutine DestroyLinOp(this) ! Does nothing but print a message to confirm destruction
        class(LinOp), intent(inout) :: this
        
        write(*,*) "Destroying linear operator with id=", this%id
        
    end subroutine
        
end module LinearOperator


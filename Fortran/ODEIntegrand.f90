module ODEIntegrand
    use Precision
    use LinearOperator
    implicit none
    private

    ! The ODEIntegrand class represents a r.h.s. of an ODE dx(t)/dt=f(x,t)
    ! This class is not abstract to avoid limitations on the use of abstract classes
    ! The default implementation is to either use a function pointer (similar to function handles in Matlab)
    ! or, if no function pointer provided, assume the r.h.s. is linear in x, f(x,t)=A(t)*x
    
    ! Abstract representation of a (generally) nonlinear function f(x,t) where x is in R^n.
    type, public :: RHSFunc
        integer  :: n=1 ! Number of x variables
        integer  :: id=0 ! An integer used to identify instances of operator
        logical, private :: ready=.false. ! Must call constructor first
        real(wp) :: tol=-1.0_wp ! Maximum relative tolerance required for operations on f
                ! Here we will assume same accuracy is needed for all operations
        procedure(EvalFunc), nopass, pointer :: eval_f=>null() ! If supplied, f(x,t) is nonlinear, use pointer to evaluate f(x,t)
        class(LinOp), pointer :: jacobian=>null() ! This is either a linear function f(x)=A(t)*x, or A=df/dx
        real(wp) :: alpha=0.0_wp ! Solve is alpha*x+f(x)=b
        logical, private :: is_linear=.false. ! Is f(x,t) linear in x
    contains
        ! Overloaded procedures from parent class LinOp
        procedure, pass(this) :: create => CreateRHS ! Create (allocate memory etc.)
        procedure, pass(this) :: init => InitRHS ! Initialize/update
        procedure, pass(this) :: destroy => DestroyRHS ! Destroy (deallocate)
        procedure, pass(this) :: eval => EvalRHS ! Evaluate f(x,t)
        procedure, pass(this) :: solve => SolveRHS ! Solve implicit system alpha*x+f(x,t)=b to tolerance
    end type

    interface
        subroutine EvalFunc(n,f,x,t,rhs)
            import ! Import definition of RHSFunc from host module
            integer, intent(in) :: n ! Number of variables
            real(wp), intent(out) :: f(:) ! Actual argument f(n)
            real(wp), intent(in) :: x(:) ! Actual argument x(n)
            real(wp), intent(in), optional :: t
            ! Internal subroutines can be used to pass extra parameters to the EvalFunc
            ! Or, better, use an object of class(RHSFunc)
            class(RHSFunc), intent(inout), optional :: rhs
        end subroutine
    end interface

    ! Add constructor of class with same name as type (convention)
    interface RHSFunc
        module procedure ConstructRHS
    end interface
    
contains
    
function ConstructRHS(n,is_linear,eval_f,jacobian,id) result(this) ! Constructor
   type(RHSFunc) :: this ! Constructed object
   integer, intent(in) :: n ! Number of variables
   logical, intent(in), optional :: is_linear ! Is f(x,t) linear in x
   procedure(EvalFunc), pointer, optional :: eval_f
   class(LinOp), intent(inout), target, optional :: jacobian
   integer, intent(in), optional :: id
   
   call this%create(n,is_linear,eval_f,jacobian,id)   
         
end function

subroutine CreateRHS(this,n,is_linear,eval_f,jacobian,id) ! Create (allocate memory, open files, etc.)
   class(RHSFunc), intent(inout) :: this ! Passed object
   integer, intent(in) :: n ! Number of variables
   logical, intent(in), optional :: is_linear ! Is f(x,t) linear in x
   procedure(EvalFunc), pointer, optional :: eval_f
   class(LinOp), intent(inout), target, optional :: jacobian
   integer, intent(in), optional :: id

   this%n = n
   if(present(is_linear)) this%is_linear=is_linear
   if(present(eval_f)) this%eval_f=>eval_f
   if(present(jacobian)) this%jacobian=>jacobian
   if(this%is_linear .and. (.not.associated(this%jacobian))) stop "Must supply the Jacobian A(t)=df/dx to Construct/CreateRHS"
   if(present(id)) this%id = id

   this%ready=.true.    

   write(*,*) "Created RHSFunc with id=", this%id
  
end subroutine

subroutine InitRHS(this,matrix,t,alpha,dt,tol) ! Initialize (set values, factorize, etc.)
   class(RHSFunc), intent(inout) :: this ! Passed object
   real(wp), intent(in), optional :: matrix(:,:) ! Data used to create the matrix   
   real(wp), intent(in), optional :: t ! If present, this is a parameter of the matrix, A<-A(t), e.g., time in an ODE
   real(wp), intent(in), optional :: alpha ! Solve is (alpha*Id+A)*x=b; default is alpha=0
   real(wp), intent(in), optional :: dt ! If present, assume ApplyExpOp is called with the same dt, and precompute exp(A*dt)
   real(wp), intent(in), optional :: tol ! Maximum relative error tolerance
   
   integer :: n
   n=this%n

   if(present(tol)) this%tol=tol
   if(present(alpha)) this%alpha=alpha 
   if(associated(this%jacobian)) call this%jacobian%init(matrix=matrix, t=t, alpha=this%alpha, dt=dt, tol=this%tol)
   
end subroutine

subroutine EvalRHS(this,f,x,t) ! Compute b=A*x
   class(RHSFunc), intent(inout) :: this
   real(wp), intent(out) :: f(:) ! Actual argument f(n)
   real(wp), intent(in) :: x(:) ! Actual argument x(n)
   real(wp), intent(in), optional :: t

   if(this%is_linear) then
      if(present(t)) then
         if(.not.(t==this%jacobian%t)) then
            write(*,*) "Re-initializing/updating Jacobian with id=", this%jacobian%id, &
                       " of RHSFunc with id=", this%id, " since t changed"
            call this%jacobian%init(t=t,alpha=this%alpha,tol=this%tol) ! Re-initialize/update matrix
         end if   
      end if   
      call this%jacobian%apply(x=x, b=f)
   else if(associated(this%eval_f)) then
      call this%eval_f(n=this%n, f=f, x=x, t=t, rhs=this)
   else
      stop "Do not know how to evaluate f(x,t)"
   end if        
      
end subroutine

subroutine SolveRHS(this,b,x,alpha) ! Solve implicit system alpha*x+f(x,t)=b to tolerance
   class(RHSFunc), intent(inout) :: this
   real(wp), intent(in) :: b(:) ! Actual argument: b(n)
   real(wp), intent(out) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(in), optional :: alpha

   real(wp) :: alpha_
   integer :: n   
   
   n=this%n
    
   alpha_=this%alpha
   if(present(alpha)) alpha_=alpha
   
   if(.not.(alpha_==this%alpha)) then
       write(*,*) "Re-initializing/updating RHSFunc with id=", this%id, " since alpha changed"
       call this%init(alpha=alpha_)
   end if            

   if(this%is_linear) then  
      call this%jacobian%solve(b=b, x=x, alpha=alpha)
   else
      stop "Do not know how to solve (alpha*I+A(t))*x=b"
   end if        
   
end subroutine

subroutine DestroyRHS(this) ! Does nothing but print a message to confirm destruction
  class(RHSFunc), intent(inout) :: this

  write(*,*) "Destroying RHSFunc with id=", this%id
    
  this%ready = .false.  
  if(associated(this%jacobian)) call this%jacobian%destroy()

end subroutine

end module ODEIntegrand


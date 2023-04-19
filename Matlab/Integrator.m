classdef Integrator < handle
%{
Integrator class solves temporal ODEs of the form x'(t) = L*x(t)+f(x,t)
Basic implementation: 
forward Euler for explicit part, backward Euler for implicit
(x^{n+1}-x^{n})/dt = L*x^{n+1}+f(x^n,t^n)
%}

 properties (Access = public)
  n   {mustBeNumeric}   % number of unknowns
  t0   {mustBeNumeric}   % initial time t0
  InitCond   (:,:) {mustBeNumeric}  % initial condition at t=t0
  LinearPart       % linear operator L (LinearOperator class)
  NonlinearPart      % integrand f  (Integrand class)
  dt   {mustBeNumeric}   % timestep

  % These properties need to be public in order for the subclasses to have access to them
  x   (:,:) {mustBeNumeric}  % current x value (starts as initial condition InitCond)
  t   {mustBeNumeric}   % current t value (starts as initial time t0)
  f_x  {mustBeNumeric}   % f(x,t) at x at time t
 end
 
 methods
 
%------------------------------------------------------------------------------------------
% Main Methods
%------------------------------------------------------------------------------------------

  function obj = Integrator(varargin)
   % Constructor: inputs n, L (LinearOperator object), f(x,t) (Integrand object), dt, t0
   
   Defaults = {[],[],[],1.0,0.0};
   Defaults(1:nargin) = varargin;
   
   if isempty(Defaults{1})
    error('Please provide number of unknowns');
   end
   
   obj.n = Defaults{1};
   obj.LinearPart = Defaults{2}; % in set function, checks class
   obj.NonlinearPart = Defaults{3}; % in set function, checks class
   
   % Set initial time, if given
   obj.t0 = Defaults{5};
   obj.t = Defaults{5};
   
   % A child class that requires exp(A*dt) must set LinearPart.dt via LinearPart.Update_dt
   % This tells the LinearOperator that the action of exp(A*dt) will be needed
   obj.dt=Defaults{4};
   
   if ~isempty(obj.LinearPart)
     obj.LinearPart.Update_alpha(-1/obj.dt);
   end
     
  end
  
  function SetInitCond(obj,x0,varargin)
   % Sets initial condition x0 at t=t0, and initializes x and t. 
   % Takes inputs x0 and t0 (optional, defaults to 0)
   
   % validate x0 is the correct size (n by 1)
   validateattributes(x0,{'double'},{'size',[obj.n,1]});
   obj.InitCond = x0;
   
   if nargin > 2 % if initial time provided, set it.
    obj.t0 = varargin{1};
   end  
   
   % Initialize internal state:  
   obj.x = x0;
   obj.t = obj.t0;
   if ~isempty(obj.NonlinearPart)    
    obj.f_x = obj.NonlinearPart.Evaluate(obj.x,obj.t);
   end
      
  end  

  function Update(obj,varargin)
   % Updates dt and/or matrix A in linear part.
   % This will recompute any necessary matrix factorizations/exponential. 
   % Otherwise, matrix factorizatins/exponentials will be reused. 
   % Optional inputs are dt and LinearPart; 
   % if LinearPart is changed, dt must be given even if it has not changed.
   if ~isempty(varargin{1})
    obj.dt = varargin{1};
    if nargin==2   % If new dt is given but not a new LinearPart
     if ~isempty(obj.LinearPart)
      obj.LinearPart.Update_alpha(-1/obj.dt);
     end
    end
   end
   if nargin>2
    obj.LinearPart.Update_LinOpA(varargin{2},-1/obj.dt);
   end
  end
  
  function Reset(obj,destroy)
   % Resets Integrator; user must call SetInitCond again before integrating
   % If destroy is true, free memory to garbage collector
   % If destroy is false, go back to the initial condition
   if(destroy) % Destroy object (makes it unusable)
    obj.t = [];
    obj.x = [];
    obj.f_x = [];   
   else % Reset back to initial condition
    obj.t = obj.t0;
    obj.x = obj.InitCond; 
    if ~isempty(obj.NonlinearPart)    
     obj.f_x = obj.NonlinearPart.Evaluate(obj.x,obj.t);
    end
   end
  end
   
%------------------------------------------------------------------------------------------
% Timestepping Methods
%------------------------------------------------------------------------------------------

  function AdvanceOneStep(obj)
   %{
   Advance by 1 step with Backward Euler for implicit part,
   Forward Euler for explicit part. Inputs are
   x^n and t^n ('previous' x and t values).
   Also returns the value f(x^n) since this may be useful for multistep methods
   %}
   
   RHS = -1/obj.dt*(obj.x);
   
   if ~isempty(obj.NonlinearPart)    
    RHS = RHS-obj.f_x;
   end
   
   if ~isempty(obj.LinearPart)
    obj.x = obj.LinearPart.SolveLinSys(RHS);
   else
    obj.x = -obj.dt*RHS;
   end
   obj.t = obj.t+obj.dt;
   if ~isempty(obj.NonlinearPart)    
    obj.f_x = obj.NonlinearPart.Evaluate(obj.x,obj.t);
   end   
   
  end % AdvanceOneStep
  
  function [X,T] = Advance(obj,Nsteps,nstep_output)
   % Advance forward by Nsteps steps, and return the solution every nstep_output steps
   % If nstep_output>0, don't output the first value (useful for running in batches of steps)
   % If nstep_output>0, also output the first value (useful if running in one batch)
   % If nstep_output==0, an empty array will be returned
   
   if(abs(nstep_output)>0)
    n_outputs=floor(Nsteps/abs(nstep_output));
    X=zeros(obj.n,n_outputs);
    T=zeros(1,n_outputs);
   else
    X=[];
    T=[];
   end

   step_out=0;
   if(nstep_output<0)
     step_out=step_out+1;
     X(:,step_out)=obj.x;
     T(:,step_out)=obj.t;
   end
     
   for step = 1:Nsteps % time loop
    obj.AdvanceOneStep;
    if(abs(nstep_output)>0) % output
     if (mod(step,abs(nstep_output))==0)
      step_out=step_out+1;
      X(:,step_out)=obj.x;
      T(:,step_out)=obj.t;
     end
    end % if output
   end % time loop
   
  end % Advance
  
%------------------------------------------------------------------------------------------
% Set Methods -- automatically called by Matlab when the value of a class property is changed
%------------------------------------------------------------------------------------------
 
  function set.LinearPart(obj,LinearPart)
   % Sets the LinearPart property with input value. 
   % Checks whether input is a LinearOperator object, unless it is empty
   if ~isempty(LinearPart)
    if ~strcmpi(class(LinearPart),'LinearOperator')
     error('Only LinearOperator objects are accepted for the LinearPart property');
    end
    obj.LinearPart = LinearPart;
   end
  end
  
  function set.NonlinearPart(obj,NonlinearPart)
   % Sets the LinearPart property with input value.
   % Checks whether input is an Integrand object, unless it is empty
   if ~isempty(NonlinearPart)
    if ~strcmpi(class(NonlinearPart),'Integrand')
    error('Only Integrand objects are accepted for the NonlinearPart property');
    end
    obj.NonlinearPart = NonlinearPart;
   end
  end  
  
 end % methods
 
end % class Integrator


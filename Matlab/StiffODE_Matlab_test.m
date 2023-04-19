clear all
close all
format long
format compact

%% Setup
nvars = 2; % This is A_eig simple 2x2 test with exact solution
exact_sol = @(t) [sin(t); cos(t)];
nstep_output=-1; % Output the solution every so many time steps (keep sign negative)

% Explicit part
B_eig = @(x,t,eig) [0.2e0 * x(2) * x(1) + ...
   cos(t) + (0.222222e-1 + 0.888889e0 * eig) * sin(t) + (0.628539e-1 - 0.314269e0 * eig) * cos(t) - 0.2e0 * sin(t) * cos(t) ; ...
                    0.3e0 * x(1) ^ 2 - 0.1e0 * x(2) ^ 2 - ...
   sin(t) + (0.628536e-1 - 0.314269e0 * eig) * sin(t) + (0.177777e0 + 0.111111e0 * eig) * cos(t) - 0.3e0 * sin(t) ^ 2 + 0.1e0 * cos(t) ^ 2];
   
% Jacobian:
dB_dx = @(x,t) [0.2e0 * x(2), 0.2e0 * x(1); 0.6e0 * x(1), -0.2e0 * x(2);];
  
% Matrix in linear (implicit) part, with largest eigenvalue eig>1
% Problem becomes more stiff as eig grows:
A_eig = @(eig) -[0.222222e-1+0.888889e0*eig,0.628539e-1-0.314269e0*eig;...
 0.628536e-1-0.314269e0*eig,0.177777e0+0.111111e0*eig];

StiffODE = @(t,x,A_eig,B,eig) A_eig(eig)*x+B(x,t,eig);

% For fully implicit methods we need the full Jacobian:
df_dx = @(eig,x,t) A_eig(eig)+dB_dx(x,t);  

%% Basic test for dt=~0.01
max_eig = 1.0; % Moderate stiffness
tfinal = 2.0*pi;
dt = tfinal/40;
t0 = 0.0;

IC = exact_sol(t0); % Exact solution at t=0
nsteps = round(tfinal/dt);
tfinal = nsteps*dt; % Make sure all integrators end at the same time

if(0) % Test out the ODE rhs
%% Solve it using ODE15s ("true solution")
dxdt = @(t,x) StiffODE(t,x,A_eig,B_eig,max_eig);
opts = odeset('RelTol',1e-11, 'AbsTol',1e-11);
[T,X] = ode15s(dxdt, [0 tfinal], IC, opts);
X=X.';
T=T.';
X_exact=exact_sol(T);

figure(1); clf
plot(T,X(1,:),'ro--', T,X_exact(1,:),'r-', ...
     T,X(2,:),'ks--', T,X_exact(2,:),'k-', ...
     'MarkerSize',4, 'DisplayName','ode15s tol=1e-11');

end % Test

%% Solve it using our own ODE solvers:
solve_flag = 0;
 % = 0  use matinv for linear solves (great for small number of variables)
 % = 1  use Cholesky for linear solves, where  (alpha*I+A_eig) is SPD (Not recommended -- in this case, losing accuracy)
 % = -1 use Cholesky for linear solves, where (-1)*(alpha*I+A_eig) is SPD (Not recommended -- in this case, losing accuracy)
 % = 2  use LU for linear solves

A = A_eig(max_eig);
stiff_part = LinearOperator(A,solve_flag);
B = @(x,t) B_eig(x,t,max_eig);
nonstiff_part = Integrand(B,nvars);

% Constructor: inputs n, L (LinearOperator object), f(x,t) (Integrand object), dt, t0
myIntegrator = Integrator(nvars, stiff_part, nonstiff_part, dt, t0);

%--- Subclass options, uncomment as necessary
% myIntegrator = ETDRK2Integrator(nvars, stiff_part, nonstiff_part, dt, t0);
% myIntegrator = AB2BDF2Integrator(nvars, stiff_part, nonstiff_part, dt, t0);
% myIntegrator = ImplicitRK2AB2Integrator(nvars, stiff_part, nonstiff_part, dt, t0);
 
%--- choices for ImplicitRKAB2Integrator, uncomment as necessary
% myIntegrator.explicit_flag = 1;  % if FE for explicit part
% myIntegrator.c = 0;    % if BE for implicit part
% myIntegrator.c = 1;    % if FE for implicit part

% Solve ODE
%------------------------

% Set initial x0 at t0
myIntegrator.SetInitCond(IC);

name=['FE+BE, nsteps=',int2str(nsteps)];
[X,T]=myIntegrator.Advance(nsteps,nstep_output);
X_exact=exact_sol(T);

% Plot solution:
%------------------------
figure(2); clf
plot(T,X(1,:),'ro--', T,X_exact(1,:),'r-', ...
     T,X(2,:),'ks--', T,X_exact(2,:),'k-', ...
     'MarkerSize',4, 'DisplayName',name);

% Now plot errors:
%------------------------
figure(3); clf
plot(T,X(1,:)-X_exact(1,:),'ro--', ...
     T,X(2,:)-X_exact(2,:),'ks--');
hold on;

% Confirm order of accuracy by resolving with half the time step size
order=1; % Should be first-order accurate for FE/BE
nsteps=nsteps*2;
dt=dt/2;
myIntegrator.Reset(0); % Reset back to initial condition
myIntegrator.Update(dt);
[X,T]=myIntegrator.Advance(nsteps,2*nstep_output);
X_exact=exact_sol(T);
plot(T,2^order*(X(1,:)-X_exact(1,:)),'rd--', ...
     T,2^order*(X(2,:)-X_exact(2,:)),'k*--');
hold on;


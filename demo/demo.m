clear all
close all
clc

addpath([pwd '/submodel']);
model_file_name = 'smallDrivelineSimScape_oscSpring_backlash';
open_system(model_file_name);
smallDrivelineSimScape_config


%% Open simulink model

Tend = 1;
Ts = 0.001;

cs = getActiveConfigSet(model_file_name);
cs.set_param('SimscapeLogDataHistory', 5000);    
cs.set_param('SimscapeLogDecimation', 1);    
cs.set_param('SimscapeLogLimitData', 'off');   % Limit data points 
cs.set_param('SimscapeLogName', 'simlog');    
cs.set_param('SimscapeLogType', 'all');
cs.set_param('SolverType', 'Fixed-step');   % Type
cs.set_param('SystemTargetFile', 'grt.tlc');   % System target file
cs.set_param('StartTime', '0.0');   % Start time 
cs.set_param('StopTime', num2str(Tend));   % Stop time 
cs.set_param('FixedStep', num2str(Ts));   % Fixed-step size (fundamental sample time) 
cs.set_param('Solver', 'ode4');   % Solver 
cs.set_param('MaxDataPoints', num2str(Tend/Ts));

set_param([model_file_name '/Solver Configuration'],'UseLocalSolver','on');
set_param([model_file_name '/Solver Configuration'],'LocalSolverSampleTime',num2str(Ts));
set_param([model_file_name '/Solver Configuration'],'DoFixedCost','on');
set_param([model_file_name '/Solver Configuration'],'MaxNonlinIter','10');
set_param([model_file_name '/Solver Configuration'],'LocalSolverChoice','NE_TRAPEZOIDAL_ADVANCER');
set_param([model_file_name '/Solver Configuration'],'ResidualTolerance','1e-9');

%% Run conversion script

% hit codegen button
out = system(['python ../run.py ' model_file_name]);
assert(out==0);
clc

rehash

%%


model = Model;

nx = model.nx;
nz = model.nz;
nu = model.nu;

disp('Unreduced DAE')
fprintf('#diff states %d\n', nx);
fprintf('#algebraic states %d\n', nz);
fprintf('#inputs states %d\n', nu);

%%
model.variable_names'

%%

figure()
spy(sparse(model.m))
figure()
spy(sparse(model.a))

%%
addpath('/home/jgillis/programs/casadi/matlab_install/casadi')
import casadi.*

x = SX.sym('x',nx);
z = SX.sym('z',nz);
u = SX.sym('u',nu);

model.f([x;z])


%% Reduce model
[Fr,xr,zr] = model.Fr;

nxr = numel(xr);
nzr = numel(zr);

disp('Reduced DAE')
fprintf('#diff states %d\n', nxr);
fprintf('#algebraic states %d\n', nzr);

%% Inspect reduced model

[Mr, rhs] = Fr(x(xr), z(zr), u);

g = rhs(nxr+1:end)

dg_dz = jacobian(g,z(zr))

jacobian(dg_dz,z(zr))

%% Create ODE model

f_ode = Mr(1:nxr,1:nxr)\rhs(1:nxr)

zsol = dg_dz\substitute(g,z(zr),0);

f_fully_explicit_ode = substitute(f_ode,z(zr),zsol)

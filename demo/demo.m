clear all
close all
clc

addpath('..')
sdl_suspension

model_file_name = 'sdl_suspension';
open_system(model_file_name);
%double_pendulum_config

%%
model_file_name = 'sdl_flexible_shaft';

%% Open simulink model

cs = getActiveConfigSet(model_file_name);

Tend = eval(cs.get_param('StopTime'));
Ts = 0.01;

N = ceil(Tend/Ts);

cs.set_param('SimscapeLogDataHistory', N);    
cs.set_param('SimscapeLogDecimation', 1);    
cs.set_param('SimscapeLogLimitData', 'off'); 
cs.set_param('SimscapeLogName', 'simlog');    
cs.set_param('SimscapeLogType', 'all');
cs.set_param('SolverType', 'Fixed-step');
cs.set_param('SystemTargetFile', 'grt.tlc');
cs.set_param('StartTime', '0.0');
cs.set_param('FixedStep', num2str(Ts)); 
cs.set_param('Solver', 'ode4');
cs.set_param('MaxDataPoints', num2str(N));

set_param([model_file_name '/Solver Configuration'],'UseLocalSolver','on');
set_param([model_file_name '/Solver Configuration'],'LocalSolverSampleTime',num2str(Ts));
set_param([model_file_name '/Solver Configuration'],'DoFixedCost','on');
set_param([model_file_name '/Solver Configuration'],'MaxNonlinIter','10');
set_param([model_file_name '/Solver Configuration'],'LocalSolverChoice','NE_TRAPEZOIDAL_ADVANCER');
set_param([model_file_name '/Solver Configuration'],'ResidualTolerance','1e-9');

simOut = sim(model_file_name,cs);

simlog = simOut.get('simlog');

t = simOut.get('tout');
if t(1)~=0
   t = [0;t]; 
end
dt = t(2)-t(1);
assert(all(diff(dt)==dt));


%% Run conversion script

rtwbuild(model_file_name)
%%
out = system(['python ../run.py ' model_file_name]);
assert(out==0);

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

model.f([x;z],u)


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

%%
X = zeros(nx,numel(t));
for i=1:nx
   data = eval(['simlog.' model.variable_names{i}]);
   values = data.series.values;
   X(i,:) = values;
end
%Z = zeros(nz,numel(t));
%for i=1:nz
%   data = eval(['simlog.' model.variable_names{nx+i}]);
%   values = data.series.values;
%   Z(i,:) = values;
%end
% Control vector
U = zeros(model.nu,numel(t));
for i=1:model.nu
   data = eval(['simlog.' model.input_names{i}]);
   values = data.series.values;
   U(i,:) = values;
end

%%
[ode_r_expl,xr,zr] = model.ode_r_expl;
nxr = numel(xr);
nzr = numel(zr);

rhsr_model = zeros(nxr,numel(t)-1);

FD = (X(:,2:end)-X(:,1:end-1))/dt;

for i=1:numel(t)-1
    [rhs] = ode_r_expl(X(xr,i),U(:,i));
    rhsr_model(:,i) = full(rhs);
end

delta_oder = FD(xr,1:end-1)-(rhsr_model(:,1:end-1)+rhsr_model(:,2:end))/2;


%%
dae_expl = model.dae_expl;

nxr = numel(xr);
nzr = numel(zr);
[dae_r_expl,xr,zr] = model.dae_r_expl;

delta_dae = zeros(model.nz,numel(t)-1);
delta_daer = zeros(nzr,numel(t)-1);

rhs_model = zeros(nx,numel(t)-1);
rhsr_model = zeros(nxr,numel(t)-1);

FD = (X(:,2:end)-X(:,1:end-1))/dt;

for i=1:numel(t)-1
    [rhs,res] = dae_expl(X(:,i),Z(:,i),U(:,i));
    delta_dae(:,i) = full(res);
    rhs_model(:,i) = full(rhs);
    
    
    [rhs,res] = dae_r_expl(X(xr,i),Z(zr,i),U(:,i));
    delta_daer(:,i) = full(res);
    rhsr_model(:,i) = full(rhs);
end

delta_ode = FD(:,1:end-1)-(rhs_model(:,1:end-1)+rhs_model(:,2:end))/2;
delta_oder = FD(:,1:end-1)-(rhsr_model(:,1:end-1)+rhsr_model(:,2:end))/2;


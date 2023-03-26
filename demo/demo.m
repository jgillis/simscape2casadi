clear all
close all
clc

addpath('../models/proprietary/R2016b/FM')
addpath('..')

smallDrivelineSimScape_oscSpring_flex_config

model_file_path = 'smallDrivelineSimScape_oscSpring_flex';
open_system(model_file_path);

%% Run conversion script

rtwbuild(model_file_path)
%%
out = system(['python ../run.py ' model_file_path]);
assert(out==0);

rehash


%%

addpath('/home/jgillis/programs/casadi/matlab_install/casadi')
model = Model;

model.input_names

%%
model.variable_names
%%
model.parameter_names
%%
model.nx
model.np
%%

import casadi.*

x = SX.sym('x',model.nx);
z = SX.sym('z',model.nz);
u = SX.sym('u',model.nu);
p = SX.sym('p',model.np);
t = SX.sym('t');

dae = model.dae_expl;

[ode, alg] = dae(x,z,u,p,t,[],[],[]);

size(ode)
size(alg)

%%
ode
alg
%%
jacobian(ode,[x;z])

%%
jacobian(alg,z)
celldisp(symvar(ans))
%%
model.variable_names{4}
%%
zsol = -jacobian(alg,z)\substitute(alg,z,0)
size(zsol)


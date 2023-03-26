clear
close all
bdclose all
clc

% addpath('../models/proprietary/R2022b/FM')
% addpath('..')

mdl_name = 'driveline_springdamper';

if ~bdIsLoaded(mdl_name)
    load_system(mdl_name);
    open_system(mdl_name);
end

%% Run conversion script

rtwbuild(mdl_name);
%%

% get the current Matlab version number
prj = currentProject;
% Get the project root folder:
projectRoot = prj.RootFolder;

model_path = get_param(mdl_name, 'FileName');
[mdl_folder,mdl_name,~] = fileparts(model_path);

syscmd = ['python ',...
              fullfile(convertStringsToChars(projectRoot),'run.py'),...
              ' ', fullfile(mdl_folder,mdl_name)];

out = system(syscmd);
assert(out==0);

rehash


%%

% addpath('/home/jgillis/programs/casadi/matlab_install/casadi')
eval(['model = ',mdl_name,'_DAE;']);
% model = Model;

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


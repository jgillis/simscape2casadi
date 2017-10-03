% Parametric: open global simscape settings (Home->Preferences),
%  and tick 'show run-time parameter settings'

clear all
close all
clc
addpath('/home/jgillis/programs/casadi/matlab_install/casadi')
import casadi.*


warning('off','symbolic:sym:sym:DeprecateExpressions')
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid')
warning('off','physmod:simscape:compiler:sli:logging:CodeGenNotSupported')
addpath('..')

% Cleanup if test folder is in a dirty state
if exist('slprj','dir')
  rmdir('slprj','s')
end
d=dir;
for d={d.name}
   d=d{1};
   if ~isempty(strfind(d,'grt'))
       rmdir(d,'s')
   end
end

models = {'../models/R2014b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper_param'};

for model_file_c=models
    model_file = model_file_c{1};
    disp(['model: ' model_file])
    [path,model_file_name,ext] = fileparts(model_file);
    AAA =1.7;
    BBB=10;

    load_system(model_file);

    % Open simulink model

    cs = getActiveConfigSet(gcs);

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
    cs.set_param('LimitDataPoints', 'off');
    cs.set_param('EnableMemcpy','off')
    cs.set_param('RTWVerbose','off');
    
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


    % Run conversion script
    disp('codegen')
    rtwbuild(model_file_name)
    disp('simscape2casadi')
    out = system(['python ../run.py ' model_file_name]);
    assert(out==0);

    rehash

    %
    model = Model;

    nx = model.nx;
    nz = model.nz;
    nu = model.nu;
    np = model.np;

    disp('Unreduced DAE')
    fprintf('#diff states %d\n', nx);
    fprintf('#algebraic states %d\n', nz);
    fprintf('#inputs %d\n', nu);
    fprintf('#parameters %d\n', np);
    
    x = SX.sym('x',nx);
    z = SX.sym('z',nz);
    u = SX.sym('u',nu);
    p = SX.sym('u',np);
    
    model.f([x;z],u,p)

    % Reduce model
    disp('reduce DAE')
    [Fr,xr,zr] = model.Fr;

    nxr = numel(xr);
    nzr = numel(zr);

    disp('Reduced DAE')
    fprintf('#diff states %d\n', nxr);
    fprintf('#algebraic states %d\n', nzr);
    %
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
    P = zeros(model.np,1);
    for i=1:model.np
       P(i) = eval(model.parameter_names{i});
    end

    [ode_r_expl,xr,zr] = model.ode_r_expl;
    nxr = numel(xr);
    nzr = numel(zr);

    rhsr_model = zeros(nxr,numel(t)-1);

    FD = (X(:,2:end)-X(:,1:end-1))/dt;

    for i=1:numel(t)-1
        [rhs] = ode_r_expl(X(xr,i),U(:,i),P);
        rhsr_model(:,i) = full(rhs);
    end

    delta_oder = FD(xr,1:end-1)-full((rhsr_model(:,1:end-1)+rhsr_model(:,2:end))/2);


    assert(max(max(abs(delta_oder)))<1e-10)
    %
    close_system(model_file_name, 0);

    rmdir([model_file_name '_grt_rtw'],'s')

    delete Model.m
    delete temp.m
    delete(model_file_name)
    rmdir('slprj','s')

end
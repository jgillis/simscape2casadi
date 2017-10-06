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

models = {...
          '../models/R2014b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper_param',...
          '../models/R2016b/driveline_springdamper_gearbox',...
          '../models/R2016b/driveline_springdamper_timedep',...
          '../models/R2016b/driveline_springdamper_backlash',...
          '../models/R2016b/driveline_springdamper_flex',...
          '../models/R2016b/driveline_ICE',...
          '../models/R2016b/fail_driveline_mass3',...
          '../models/R2016b/fail_driveline_mass'}; % Mass matrix has entries outside of 1:nxr
          %../models/R2016b/fail_driveline_springdamper_flex',...
          %'../models/R2016b/driveline_EM',...
          %'../models/R2016b/proprietary/FM/DCT_v_0_1'};%,...
          %'../models/R2016b/driveline_springdamper_LUT'};

for model_file_c=models
    model_file = model_file_c{1};
    disp(['model: ' model_file])
    [path,model_file_name,ext] = fileparts(model_file);
    
    if exist([model_file '_config.m'])
        run([model_file '_config.m'])
    end

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
    max_iter = 10;
    set_param([model_file_name '/Solver Configuration'],'MaxNonlinIter',num2str(max_iter));
    
    % According to Simscape documentation, backward Euler is activated at
    % the initial time, and at discrete changes
    %set_param([model_file_name '/Solver Configuration'],'LocalSolverChoice','NE_TRAPEZOIDAL_ADVANCER');
    set_param([model_file_name '/Solver Configuration'],'LocalSolverChoice','NE_BACKWARD_EULER_ADVANCER');
    set_param([model_file_name '/Solver Configuration'],'ResidualTolerance','1e-9');

    simOut = sim(model_file_name,cs);

    simlog = simOut.get('simlog');

    ts = simOut.get('tout');
    if ts(1)~=0
       ts = [0;ts]; 
    end
    dt = ts(2)-ts(1);
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
    nq = model.nq;

    x = SX.sym('x',nx);
    z = SX.sym('z',nz);
    u = SX.sym('u',nu);
    p = SX.sym('u',np);
    q = SX.sym('u',nq);
    t = SX.sym('t');
    
    %
    X = zeros(nx,numel(ts));
    for i=1:nx
       data = eval(['simlog.' model.variable_names{i}]);
       values = data.series.values;
       X(i,:) = values;
    end
    Z = zeros(nz,numel(ts));
    for i=1:nz
       data = eval(['simlog.' model.variable_names{nx+i}]);
       values = data.series.values;
       Z(i,:) = values;
    end
    % Control vector
    U = zeros(model.nu,numel(ts));
    for i=1:model.nu
       data = eval(['simlog.' model.input_names{i}]);
       values = data.series.values;
       U(i,:) = values;
    end
    P = zeros(model.np,1);
    for i=1:model.np
       P(i) = eval(model.parameter_names{i});
    end


    [dae_r_expl,xr,zr] = model.dae_r_expl;
    nxr = numel(xr);
    nzr = numel(zr);

    Xtraj = zeros(nxr,numel(ts));
    Ztraj = zeros(nzr,numel(ts));
        
    Xtraj(:,1) = X(xr,1);
    Ztraj(:,1) = Z(zr,1);

    %xn = x+dt*f(xn,zn)
    
    xnext = SX.sym('xnext',nxr);
    znext = SX.sym('znext',nzr);
    tnext = SX.sym('t');
    
    [rhs,res] = dae_r_expl(xnext,znext,u,p,q,tnext);
    rf = Function('rf',{[xnext;znext],x(xr),u,p,q,tnext},{[xnext-(x(xr)+dt*rhs);res]});
    rf = rootfinder('rf','newton',rf,struct('max_iter',max_iter));
    
    for i=1:numel(ts)-1
      % This is odd: Simscape seems to take U(:,i+1) here instead of U(:,i)
      sol = rf([Xtraj(:,i);Ztraj(:,i)],Xtraj(:,i),U(:,i+1),P,0,ts(i+1));
      sol = full(sol);
      Xtraj(:,i+1) = sol(1:nxr);
      Ztraj(:,i+1) = sol(nxr+1:end);
    end
    
    assert(max(max(abs(X(xr,:)-Xtraj)))<1e-10)
    
    %
    close_system(model_file_name, 0);

    rmdir([model_file_name '_grt_rtw'],'s')

    delete Model.m
    delete temp.m
    if exist(model_file_name)
        delete(model_file_name)
    end
    if exist([model_file_name '.exe'])
        delete([model_file_name '.exe'])
    end
    rmdir('slprj','s')

end

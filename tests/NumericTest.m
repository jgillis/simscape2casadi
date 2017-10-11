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

has_proprietary = exist('../models/proprietary','dir');


models = {...
          '../models/R2014b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper',...
          '../models/R2016b/driveline_springdamper_param',...
          '../models/R2016b/driveline_springdamper_gearbox',...
          '../models/R2016b/driveline_springdamper_timedep',...
          '../models/R2016b/driveline_springdamper_backlash',...
          '../models/R2016b/driveline_springdamper_flex',...
          '../models/R2016b/driveline_EM',...
          '../models/R2016b/fail_driveline_mass3',...
          '../models/R2016b/fail_driveline_mass',...
          '../models/R2016b/driveline_spring2',...
          '../models/R2016b/driveline_ICE',...
          '../models/R2016b/fail_driveline_springdamper_flex',...
          '../models/R2016b/driveline_springdamper_LUT',...
          '../models/R2016b/driveline_springdamper_constant',...
          '../models/proprietary/R2016b/PMA/DiskClutch_Inertia_New'};
       %   '../models/R2016b/proprietary/FM/smallDrivelineSimScape_oscSpring_flex'};
          %'../models/R2016b/driveline_springdamper_LUT'};
          % '../models/R2016b/proprietary/PMA/DiskClutch_Inertia_New'};
         % %'../models/R2016b/driveline_springdamper_clutch'};%,...
          %'../models/R2016b/proprietary/FM/DCT_v_0_1'};%,...

for model_file_c=models
    
    % Skip priorietary models when not available
    if ~isempty(strfind(path,'proprietary')) && ~has_proprietary
        warning(['Skipping proprietary model ' model_file])
        continue
    end
    
    model_file = model_file_c{1};
    disp(['model: ' model_file])
    [path,model_file_name,ext] = fileparts(model_file);
    
    load_system(model_file);

    % Open simulink model

    cs = getActiveConfigSet(gcs);

    Tend = eval(cs.get_param('StopTime'));
    eval_check = true;
    integration_check = true;
    skip_ode = false;
    check_ode = true;
    Ts = 0.01;    
    if exist([model_file '_config.m'])
        run([model_file '_config.m'])
    end
    assert(integration_check || eval_check);

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
    cs.set_param('Toolchain','Automatically locate an installed toolchain');
    cs.set_param('GenCodeOnly',true);
    cs.set_param('SaveTime',true);
    
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

    ts = simOut.get(cs.get_param('TimeSaveName'));
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
    disp('Unreduced DAE')
    fprintf('#diff states %d\n', nx);
    fprintf('#algebraic states %d\n', nz);
    fprintf('#inputs %d\n', nu);
    fprintf('#parameters %d\n', np);
    fprintf('#major modes %d\n', nq);
    x = SX.sym('x',nx);
    z = SX.sym('z',nz);
    u = SX.sym('u',nu);
    p = SX.sym('u',np);
    q = SX.sym('q',nq);
    t = SX.sym('t');
    
    model.f([x;z],u,p,q,t)

    % Reduce model
    disp('reduce DAE')
    [Fr,xr,zr] = model.Fr;

    nxr = numel(xr);
    nzr = numel(zr);

    disp('Reduced DAE')
    fprintf('#diff states %d\n', nxr);
    fprintf('#algebraic states %d\n', nzr);
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

    if eval_check
        dae_expl = model.dae_expl;

        [dae_r_expl,xr,zr] = model.dae_r_expl;

        delta_dae = zeros(model.nz,numel(ts)-1);
        delta_daer = zeros(nzr,numel(ts)-1);

        rhs_model = zeros(nx,numel(ts)-1);
        rhsr_model = zeros(nxr,numel(ts)-1);

        FD = (X(:,2:end)-X(:,1:end-1))/dt;

        for i=1:numel(ts)-1
            [rhs,res] = dae_expl(X(:,i),Z(:,i),U(:,i),P,0,ts(i));
            delta_dae(:,i) = full(res);
            rhs_model(:,i) = full(rhs);


            [rhs,res] = dae_r_expl(X(xr,i),Z(zr,i),U(:,i),P,0,ts(i));
            delta_daer(:,i) = full(res);
            rhsr_model(:,i) = full(rhs);
        end

        %if ~isempty(delta_dae)
        %  assert(max(max(abs(delta_dae)))<1e-10)
        %end
        %delta_ode = FD(:,1:end-1)-rhs_model(:,2:end);
        %assert(max(max(abs(delta_ode)))<1e-10)

        if ~isempty(delta_daer)
          assert(max(max(abs(delta_daer)))<1e-10)
        end

        delta_oder = FD(xr,1:end-1)-rhsr_model(:,2:end);
        assert(max(max(abs(delta_oder)))<1e-10)

        if check_ode
            model.ode_expl;


            [ode_r_expl,xr,zr,unsafe] = model.ode_r_expl;

            ode_r_expl

            nxr = numel(xr);
            nzr = numel(zr);

            rhsr_model = zeros(nxr,numel(ts)-1);

            FD = (X(:,2:end)-X(:,1:end-1))/dt;

            for i=1:numel(ts)-1
                [rhs] = ode_r_expl(X(xr,i),U(:,i),P,0,ts(i));
                rhsr_model(:,i) = full(rhs);
            end

            % trapezoidal
            % delta_oder = FD(xr,1:end-1)-full((rhsr_model(:,1:end-1)+rhsr_model(:,2:end))/2);
            delta_oder = FD(xr,1:end-1)-rhsr_model(:,2:end);

            if unsafe
                disp('Unsafe ODE; ignoring non-regular values in numeric test')
                delta_oder(isnan(delta_oder)) = 0;
                delta_oder(isinf(delta_oder)) = 0;
            end

            assert(max(max(abs(delta_oder)))<1e-10)
        end
    end
    
    if integration_check
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

        assert(max(max(abs(X(xr,:)-Xtraj)))<1e-9)
    end
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
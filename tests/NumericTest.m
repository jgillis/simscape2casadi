% Parametric: open global simscape settings (Home->Preferences),
%  and tick 'show run-time parameter settings'

clear
close all
clc
% addpath('/home/jgillis/programs/casadi/matlab_install/casadi')
import casadi.*


warning('off','symbolic:sym:sym:DeprecateExpressions')
warning('off','symbolic:generate:FunctionNotVerifiedToBeValid')
warning('off','physmod:simscape:compiler:sli:logging:CodeGenNotSupported')
% addpath('..')

% get the current Matlab version number
p = currentProject;
% Get the project root folder:
projectRoot = p.RootFolder;

% Cleanup if test folder is in a dirty state
if exist('slprj','dir')==1
  rmdir('slprj','s')
end
d=dir;
for d={d.name}
   d=d{1};
   if contains(d,'grt')
       rmdir(d,'s')
   end
end

has_proprietary = exist('../models/proprietary','dir');

model_folder = '../models/R2022b/';

models = {'driveline_springdamper',...
          'driveline_springdamper_param',...
          'driveline_springdamper_params',...
          'driveline_springdamper_gearbox',...
          'driveline_springdamper_timedep',...
          'driveline_springdamper_backlash',...
          'driveline_springdamper_flex',...
          'driveline_EM',...
          'fail_driveline_mass3',...
          'fail_driveline_mass',...
          'driveline_spring2',...
          'driveline_springdamper_constant'};
          ...%'driveline_ICE',...
          ...%'driveline_ICE_fail2',...
          ...%'fail_driveline_springdamper_flex',...
          ...%'driveline_springdamper_LUT',...
          ...%'driveline_springdamper_torqueconv',...
          ...%'driveline_springdamper_LUT2'}; 


for model_file_c=models
    
    % Skip priorietary models when not available
    if contains(path,'proprietary') && ~has_proprietary
        warning(['Skipping proprietary model ' model_file])
        continue
    end
    
    model_file = model_file_c{1};
    disp(['model: ' model_file])
    [path,model_file_path,ext] = fileparts(model_file);
%     addpath(path);
    
    eval_check = true;
    integration_check = true;
    skip_ode = false;
    check_ode = true;
    output_check = true;
    Ts = 0.01;    
    if exist([model_file '_config.m'],'file')
        run([model_file '_config.m'])
    end
    assert(integration_check || eval_check);
    
    load_system(model_file);

    % Open simulink model

    cs = getActiveConfigSet(gcs);

    Tend = eval(cs.get_param('StopTime'));


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
    
    set_param([model_file_path '/Solver Configuration'],'UseLocalSolver','on');
    set_param([model_file_path '/Solver Configuration'],'LocalSolverSampleTime',num2str(Ts));
    set_param([model_file_path '/Solver Configuration'],'DoFixedCost','on');
    max_iter = 10;
    set_param([model_file_path '/Solver Configuration'],'MaxNonlinIter',num2str(max_iter));
    
    % According to Simscape documentation, backward Euler is activated at
    % the initial time, and at discrete changes
    %set_param([model_file_path '/Solver Configuration'],'LocalSolverChoice','NE_TRAPEZOIDAL_ADVANCER');
    set_param([model_file_path '/Solver Configuration'],'LocalSolverChoice','NE_BACKWARD_EULER_ADVANCER');
    set_param([model_file_path '/Solver Configuration'],'ResidualTolerance','1e-9');

    simOut = sim(model_file_path,cs);

    simlog = simOut.get('simlog');

    logsout = simOut.get('logsout');

    ts = simOut.get(cs.get_param('TimeSaveName'));
    if ts(1)~=0
       ts = [0;ts]; 
    end
    dt = ts(2)-ts(1);
    assert(all(diff(dt)==dt));


    % Run conversion script
    disp('codegen')
    rtwbuild(model_file_path)
    disp('simscape2casadi')
    model_path = get_param(model_file_path, 'FileName');
    [mdl_folder,mdl_name,~] = fileparts(model_path);
    syscmd = ['python ',...
              fullfile(convertStringsToChars(projectRoot),'run.py'),...
              ' ', fullfile(mdl_folder,mdl_name)];
    out = system(syscmd);
    assert(out==0);

    rehash

    %
    eval(['model = ',model_file_c{1},'_DAE;']);

    nx = model.nx;
    nz = model.nz;
    nu = model.nu;
    np = model.np;
    nq = model.nq;
    nw = model.nw;
    ns = model.ns;
    ny = model.ny;
    disp('Unreduced DAE')
    fprintf('#diff states %d\n', nx);
    fprintf('#algebraic states %d\n', nz);
    fprintf('#inputs %d\n', nu);
    fprintf('#parameters %d\n', np);
    fprintf('#major modes %d\n', nq);
    fprintf('#delays %d\n', nw);
    x = SX.sym('x',nx);
    z = SX.sym('z',nz);
    u = SX.sym('u',nu);
    p = SX.sym('u',np);
    t = SX.sym('t');
    q = SX.sym('q',nq);
    w = SX.sym('w',nw);
    s = SX.sym('s',ns);
    
    model.f([x;z],u,p,t,q,w,s)

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
        %         data = eval(['simlog.' model.input_names{i}]);
        %         values = data.series.values;
        %         U(i,:) = values;
        dataset = get(logsout,['In',num2str(i)]);
        U(i,:) = dataset.Values.Data;
    end
    P = zeros(model.np,1);
    j=0; para_path_old = ''; para_name_old = '';
    for i=1:model.np
        %        P(i) = eval(model.parameter_names{i});
        paraVal = eval([get_param(model.parameter_paths{i},model.parameter_names{i})]);
        if ~strcmp(model.parameter_paths{i},para_path_old) ||...
           ~strcmp(model.parameter_names{i},para_name_old)
            j=i;
        end
        if numel(paraVal)==1
            P(i) = paraVal;
        else
            [row,col] = ind2sub(size(paraVal),i-j+1);
            P(i) = paraVal(row,col);
        end
        para_path_old = model.parameter_paths{i};
        para_name_old = model.parameter_names{i};
    end
    % Outputs
    Y = zeros(model.ny,numel(ts));
    for i=1:model.ny
        %        data = eval(['simlog.' model.output_names{i}]);
        %        values = data.series.values;
        %        Y(i,:) = values;
        dataset = get(logsout,['Out',num2str(i)]);
        Y(i,:) = dataset.Values.Data;
    end
    if eval_check
        dae_expl = model.dae_expl;

        [dae_r_expl,xr,zr] = model.dae_r_expl;

        delta_dae = zeros(model.nz,numel(ts)-1);
        delta_daer = zeros(nzr,numel(ts)-1);

        rhs_model = zeros(nx,numel(ts)-1);
        rhsr_model = zeros(nxr,numel(ts)-1);

        y_model = zeros(ny,numel(ts)-1);
        yr_model = zeros(ny,numel(ts)-1);

        FD = (X(:,2:end)-X(:,1:end-1))/dt;

        for i=1:numel(ts)-1
            [rhs,res,y] = dae_expl(X(:,i),Z(:,i),U(:,i),P,ts(i),0,0,0);
            delta_dae(:,i) = full(res);
            rhs_model(:,i) = full(rhs);
            y_model(:,i) = full(y);

            [rhs,res,y] = dae_r_expl(X(xr,i),Z(zr,i),U(:,i),P,ts(i),0,0,0);
            delta_daer(:,i) = full(res);
            rhsr_model(:,i) = full(rhs);
            yr_model(:,i) = full(y);
        end

        if ~isempty(delta_dae)
          assert(max(max(abs(delta_dae)))<1e-10)
        end
        delta_ode = FD(:,1:end-1)-rhs_model(:,2:end);
        assert(max(max(abs(delta_ode)))<1e-10)

        if ~isempty(delta_daer)
          assert(max(max(abs(delta_daer)))<1e-10)
        end

        delta_oder = FD(xr,1:end-1)-rhsr_model(:,2:end);
        assert(max(max(abs(delta_oder)))<1e-10)
        if ny>0
            assert(max(max(abs(Y(:,1:end-1)-y_model)))<1e-8)
            assert(max(max(abs(Y(:,1:end-1)-yr_model)))<1e-8)
        end
        if check_ode
            model.ode_expl;


            [ode_r_expl,xr,zr,unsafe] = model.ode_r_expl;

            nxr = numel(xr);
            nzr = numel(zr);

            rhsr_model = zeros(nxr,numel(ts)-1);
            yr_model = zeros(ny,numel(ts)-1);

            FD = (X(:,2:end)-X(:,1:end-1))/dt;

            for i=1:numel(ts)-1
                [rhs,y] = ode_r_expl(X(xr,i),U(:,i),P,ts(i),0,0,0);
                rhsr_model(:,i) = full(rhs);
                yr_model(:,i) = full(y);
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
            DY = Y(:,1:end-1)-yr_model;
            DY(isnan(DY) | isinf(DY)) = 0;
            if output_check && ny>0
              assert(max(max(abs(DY)))<1e-8)
            end
        end
    end
    
    if integration_check
        [dae_r_expl,xr,zr] = model.dae_r_expl;
        nxr = numel(xr);
        nzr = numel(zr);

        Xtraj = zeros(nxr,numel(ts));
        Ztraj = zeros(nzr,numel(ts));
        Ytraj = zeros(ny,numel(ts));

        Xtraj(:,1) = X(xr,1);
        Ztraj(:,1) = Z(zr,1);


        %xn = x+dt*f(xn,zn)

        xnext = SX.sym('xnext',nxr);
        znext = SX.sym('znext',nzr);
        tnext = SX.sym('t');

        [rhs,res,y] = dae_r_expl(xnext,znext,u,p,tnext,q,w,s);
        rf_res = Function('rf',{[xnext;znext],x(xr),u,p,tnext,q,w,s},{[xnext-(x(xr)+dt*rhs);res]});
        outputf = Function('outputf',{xnext,znext,u,p,tnext,q,w,s},{y});
        rf = rootfinder('rf','newton',rf_res,struct('max_iter',max_iter));
        %rf = rootfinder('rf','nlpsol',rf_res,struct('nlpsol','ipopt','nlpsol_options',struct('ipopt',struct('print_level',0),'print_time',false)));

        disp('integration_check')
        for i=1:numel(ts)-1
          % This is odd: Simscape seems to take U(:,i+1) here instead of U(:,i)
          [sol] = rf([Xtraj(:,i);Ztraj(:,i)],Xtraj(:,i),U(:,i+1),P,ts(i+1),0,0,0);
          
          %err = norm(full(rf_res(sol, Xtraj(:,i),U(:,i+1),P,ts(i+1),0,0,0)));
          % rf_res([Xtraj_orig(xr,i+1);Ztraj_orig(zr,i+1)], Xtraj_orig(:,i),U(:,i+1),P,ts(i+1),0,0,0)
          %assert(err<1e-8);
          sol = full(sol);
          Xtraj(:,i+1) = sol(1:nxr);
          Ztraj(:,i+1) = sol(nxr+1:end);
        end
        
        for i=1:numel(ts)
           Ytraj(:,i) = full(outputf(Xtraj(:,i),Ztraj(:,i),U(:,i),P,ts(i),0,0,0));
        end

        assert(max(max(abs(X(xr,1:i)-Xtraj(:,1:i))))<1e-8)
        if ny>0 && output_check
          assert(max(max(abs(Y(:,1:i)-Ytraj(:,1:i))))<1e-8)
        end
    end
    
    %
    close_system(model_file_path, 0);

%     rmdir([model_file_path '_grt_rtw'],'s')
% 
%     delete Model.m
    delete temp.m
%     if exist(model_file_path,'dir')
%         delete(model_file_path)
%     end
%     if exist([model_file_path '.exe'],'file')
%         delete([model_file_path '.exe'])
%     end
%     rmdir('slprj','s')
% 
%     rmpath(path);
end

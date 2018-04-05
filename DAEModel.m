%{
MIT License

Copyright (c) 2017 Joris Gillis, KU Leuven

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%}
classdef DAEModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %
    % variables:
    %
    % x -- differential states
    % z -- algebraic states
    % u -- controls
    % p -- parameters
    % t -- time
    % q -- currently not in use; pass 0
    % w -- delay pipeline
    % s -- start time of simulation?
    %
    %
    % Delay:
    %
    % Delays are describe by a nw-vector-valued delay pipeline
    %
    % Each entry in the delay pipeline is treated separately.
    % More information about each entry is obtained from a couple of helper
    % functions, which are also nw-vector-valued:
    %
    %  del_t(p)        delay length [s]
    %  del_tmax(p)     maximum delay length [s] 
    %  del_v(x,u,p,t)  input of the delay pipeline
    %  del_v0(p)       output of the delay pipeline when the input hasn't
    %                  reached the output yet
    
    properties
    end
    
    methods
        function self = DAEModel()
        end
        function out = n(self)
          out = size(self.sp_a,1);
        end
        function out = ny(self)
          out = numel(self.output_names);
        end
        function out = np(self)
          out = numel(self.parameter_names);
        end
        function out = nq(self)
          out = numel(self.mmode_names);
        end
        function out = nz(self)
          out = self.n - self.nx;
        end
        function out = nx(self)
          out = nnz(sum(casadi.DM(self.sp_m,1),2));
        end
        function out = nu(self)
          out = size(self.sp_b,2);
        end
        function out = dae_expl(self)
          import casadi.*
          x = SX.sym('x',self.nx);
          z = SX.sym('z',self.nz);
          u = SX.sym('u',self.nu);
          p = SX.sym('p',self.np);
          t = SX.sym('t');
          q = SX.sym('q',self.nq);
          w = SX.sym('w',self.nw);
          s = SX.sym('s',self.ns);
          
          % M [dot(x);dot(z)] = E(x,z,u,p,t,q,w,s)

          E = self.a(p) * [x;z] + self.b(p)*u + self.f([x;z],u,p,t,q,w,s);
          M = self.m(p);
          M = M(1:self.nx,1:self.nx);
          Y = self.y([x;z],u,p,t);
          
          out = Function('E',{x,z,u,p,t,q,w,s},{M\E(1:self.nx),E(self.nx+1:end),Y},{'x','z','u','p','t','q','w','s'},{'ode','alg','y'});
        end
        function [out,xr,zr] = dae_r_expl(self)
          import casadi.*
          [Fun, xr, zr] = self.Fr();
          x = SX.sym('x',self.nx);
          z = SX.sym('z',self.nz);
          u = SX.sym('u',self.nu);
          p = SX.sym('p',self.np);
          t = SX.sym('t');
          q = SX.sym('q',self.nq);
          w = SX.sym('w',self.nw);
          s = SX.sym('s',self.ns);
          
          [M,rhs,y] = Fun(x(xr),z(zr),u,p,t,q,w,s);
          

          nxr = numel(xr);
          M = M(1:nxr,1:nxr);

          res = rhs(nxr+1:end);
          rhs = rhs(1:nxr);
          
          r = M\rhs;
          for i=1:length(r)
             if is_constant(r(i))
                 assert(is_regular(r(i)), 'M not invertible; may happen when you connect flexible element in series without inertia inbetween')
             end
          end
          out = Function('E',{x(xr),z(zr),u,p,t,q,w,s},{r,res,y},{'xr','zr','u','p','t','q','w','s'},{'ode','alg','y'});
        end
        function [out,unsafe] = ode_expl(self)
          import casadi.*
          dae = self.dae_expl();
          
          dae_in = sx_in(dae);
          x = dae_in{1};
          z = dae_in{2};
          u = dae_in{3};
          p = dae_in{4};
          t = dae_in{5};
          q = dae_in{6};
          w = dae_in{7};
          s = dae_in{8};
          [rhs,res,y] = dae(dae_in{:});
          
          % Check if res in linear in z
          if is_empty(res)
              unsafe = 0;
          else
            assert(~any(which_depends(res, z, 2, false)))
              J = jacobian(res,z);

              unsafe = ~isempty(strfind(str(J),'?'));

              zsol = -J\substitute(res,z,0);
              rhs = substitute(rhs,z,zsol);
              y = substitute(y,z,zsol);
          end
          
          out = Function('E',{x,u,p,t,q,w,s},{rhs,y},char('x','u','p','t','q','w','s'),{'rhs','y'});
        end
        function [out,xr,zr,unsafe] = ode_r_expl(self)
          import casadi.*
          [dae,xr,zr] = self.dae_r_expl();
          
          dae_in = sx_in(dae);
          x = dae_in{1};
          z = dae_in{2};
          u = dae_in{3};
          p = dae_in{4};
          t = dae_in{5};
          q = dae_in{6};
          w = dae_in{7};
          s = dae_in{8};
          
          [rhs,res,y] = dae(dae_in{:});
          
          % Check if res in linear in z
          size(res)
          size(z)
          if is_empty(res)
              unsafe = 0;
          else
            assert(~any(which_depends(res, z, 2, false)))
            J = jacobian(res,z);
            unsafe = ~isempty(strfind(str(J),'?'));
            zsol = -J\substitute(res,z,0);
            rhs = substitute(rhs,z,zsol);
            y = substitute(y,z,zsol);
          end
          
          out = Function('E',{x,u,p,t,q,w,s},{rhs,y},{'x','u','p','t','q','w','s'},{'rhs','y'});
        end
        function [Fun, xr, zr] = Fr(self)
            model = Model;

            nx = self.nx;
            nz = self.nz;
            nu = self.nu;
            np = self.np;
            nq = self.nq;
            nw = self.nw;
            ns = self.ns;

            import casadi.*
            x = SX.sym('x',nx);
            x_orig = x;
            z = SX.sym('z',nz);
            z_orig = z;
            u = SX.sym('u',nu);
            p = SX.sym('p',np);
            t = SX.sym('t');
            q = SX.sym('q',nq);
            w = SX.sym('w',nw);
            s = SX.sym('s',ns);
            
            x_old = SX(x);
            z_old = SX(z);

            pp = p;

            M = self.m(p);

            M = M(1:nx,1:nx);

            A = self.a(p);
            B = self.b(p);

            f_expr = self.f([x;z],u,p,t,q,w,s);

            patt = @(A) full(DM(sparsity(sparsify(SX(A))),1));
            
            % M [dot(x);dot(z)] = E

            % Find trivial algebraic variables: E_(nx+i) = z_i + b_j u + f_(nx+i)(not on z_i)

            no_dependence = ~which_depends(f_expr,z,1,true)';
            assert(~depends_on(f_expr(find(no_dependence)),z));

            candidates = sum(patt(A),2)==1 & sum(patt(A(:,1:nx)),2)==0;
            candidates = candidates & no_dependence;

            % z-indices of 
            e_triv = find(candidates(nx+1:end));
            if ~isempty(e_triv)
                disp('Eliminating trivial algebraic variables')
                e_trivi = find(~candidates(nx+1:end));

                [col,~] = find(patt(A(nx+e_triv,:))');
                z_triv = col-nx;
                z_val = -f_expr(nx+e_triv)-B(nx+e_triv,:)*u;
                
                for i=1:numel(e_triv)
                   z_val(i) = z_val(i)/A(nx+e_triv(i),nx+z_triv(i)); 
                end

                z_old(z_triv) = z_val;
                                
                c = (1:nx+nz);
                c(col) = [];
                A = [A(1:nx,:);A(nx+e_trivi,:)];
                f_expr = substitute([f_expr(1:nx);f_expr(nx+e_trivi)],z(z_triv),z_val) + A(:,col)*z_val;
                A = A(:,c);
                B = [B(1:nx,:);B(nx+e_trivi,:)];

                zsel = c(nx+1:end)-nx;
                z = z(zsel);
                nz = numel(z);
            end
          
            t_mupad = sym('t');
            x_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(x),'uni',false);
            dx_mupad = cell(nx,1);
            for i=1:nx
               dx_mupad{i} = diff(x_mupad{i},t_mupad);
            end
            z_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(z),'uni',false);
            u_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(u),'uni',false);
            p_mupad = cellfun(@(e) sym(name(e)),vertsplit(p),'uni',false);
            q_mupad = cellfun(@(e) sym(name(e)),vertsplit(q),'uni',false);
            w_mupad = cellfun(@(e) sym(name(e)),vertsplit(w),'uni',false);
            s_mupad = cellfun(@(e) sym(name(e)),vertsplit(s),'uni',false);

            [A_mupad, fA_mupad, fA_casadi] = DAEModel.to_mupad(A, 'A');
            [B_mupad, fB_mupad, fB_casadi] = DAEModel.to_mupad(B, 'B');
            [M_mupad, fM_mupad, fM_casadi] = DAEModel.to_mupad(M, 'M');
            [f_expr_mupad, fF_mupad, fF_casadi] = DAEModel.to_mupad(f_expr, 'F');
            
            vars = [vertcat(x_mupad{:});vertcat(z_mupad{:})];

            eqs = [M_mupad*vertcat(dx_mupad{:});zeros(nz,1)]==A_mupad*vars+B_mupad*vertcat(u_mupad{:})+f_expr_mupad;



            %vars = vertcat(x_mupad{:});
            %[newEqs,newVars,~] = reduceRedundancies(eqs,vars);
            
            %assert(isLowIndexDAE(newEqs,newVars))
            

            [newEqs,newVars,R] = reduceDAEIndex(eqs,vars);
            %assert(numel(R)==0)
            [newEqs,newVars,R] = reduceRedundancies(newEqs,newVars);
            
            replaced_variables = [R.constantVariables;R.replacedVariables];
            
            
            replacement_x = [];
            replacement_z = [];
            
            replaced_variables_x = [];
            replaced_variables_z = [];

            replaced_labels = replaced_variables(:,1);
            for i=1:numel(replaced_labels)
              label = char(replaced_labels(i));
              label = label(1:end-3);
              var_type = label(1);
              var_index = str2num(label(3:end));
              if (var_type=='x')
                replacement_x = [replacement_x;var_index+1];
                replaced_variables_x = [replaced_variables_x;replaced_variables(i,2)];
              elseif (var_type=='z')
                replacement_z = [replacement_z;var_index+1];
                replaced_variables_z = [replaced_variables_z;replaced_variables(i,2)];
              end
            end
            
            replaced_variables = [replaced_variables_x;replaced_variables_z];

            fE_mupad=[];
            fD_mupad = [];
            deriv_info = DAEModel.collect_derivatives([newEqs;replaced_variables]);
            fE_casadi=cell(length(deriv_info),1);
            for i=1:length(deriv_info)
               funname = deriv_info{i}{1};
               arg = deriv_info{i}{2};
               f = eval(['f' funname(3) '_casadi{' funname(4:end) '}']);
               f_mupad = eval(['f' funname(3) '_mupad{' funname(4:end) '}']);
               
               args = arrayfun(@char,children(f_mupad),'uni',false);
               
               n={};
               for j=1:n_in(f)
                 n{j} = name_in(f,j-1);
               end

               if numel(args)==1
                 fD_mupad=[fD_mupad;sym(['D(' funname ')(' strjoin(args,',') ')'])];  
               else
                 fD_mupad=[fD_mupad;sym(['D([' num2str(arg) '], ' funname ')(' strjoin(args,',') ')'])];
               end
               fE_casadi{i}=f.factory(['f_E' num2str(i)],n,{['jac:out:' n{arg}]});
               fE_mupad = [fE_mupad;sym(['f_E' num2str(i) '(' strjoin(args,',') ')'])];
            end
            newEqs = subs(newEqs,fD_mupad,fE_mupad);


            assert(isLowIndexDAE(newEqs,newVars))
            
            [M,F] = massMatrixForm(newEqs,newVars);

            XZ = sym('XZ',[length(newVars),1]);
            
            U = sym('U',[nu,1]);
            P = sym('P',[np,1]);
            T = sym('T',[1,1]);
            Q = sym('Q',[nq,1]);
            W = sym('W',[nw,1]);
            S = sym('S',[ns,1]);  
            unsafe = arrayfun(@(x) ~isempty(strfind(char(x),'diff')),F);
            if any(unsafe)
                warning('diff ignored')
            end
            F = subs(F,[newVars;vertcat(u_mupad{:});vertcat(p_mupad{:});vertcat(q_mupad{:});vertcat(w_mupad{:});vertcat(s_mupad{:})],[XZ;U;P;Q;W;S]);
            F(unsafe) = NaN;
            unsafe = arrayfun(@(x) ~isempty(strfind(char(x),'diff')),M);
            if any(unsafe)
                warning('diff ignored')
            end
            M = subs(M,[newVars;vertcat(u_mupad{:});vertcat(p_mupad{:});vertcat(q_mupad{:});vertcat(w_mupad{:});vertcat(s_mupad{:})],[XZ;U;P;Q;W;S]);
            M(unsafe) = NaN;
            unsafe = arrayfun(@(x) ~isempty(strfind(char(x),'diff')),replaced_variables);
            if any(unsafe)
                warning('diff ignored')
            end
            replaced_variables = subs(replaced_variables,[newVars;vertcat(u_mupad{:});vertcat(p_mupad{:});vertcat(q_mupad{:});vertcat(w_mupad{:});vertcat(s_mupad{:})],[XZ;U;P;Q;W;S]);
            replaced_variables(unsafe) = NaN;

            names_newVars = cellfun(@char,num2cell(newVars),'uni',false);
            names_x = cellfun(@char,x_mupad,'uni',false);
            names_z = cellfun(@char,z_mupad,'uni',false);
            
            xr = [];
            zr = [];
            
            for i=1:length(names_newVars)
                found = false;
                for j=1:length(names_x)
                    if strcmp(names_newVars{i},names_x{j})
                       xr = [xr;j]; 
                       found=true;
                    end
                end
                for j=1:length(names_z)
                    if strcmp(names_newVars{i},names_z{j})
                       zr = [zr;j]; 
                       found=true;
                    end
                end
                if ~found
                    zr = [zr; nz+1];
                    n = names_newVars{i};
                    z = [z;SX.sym(n(1:end-3))];
                    nz = numel(z);
                end
            end

            F = matlabFunction(F,'Vars',{XZ,U,P,t_mupad,Q,W,S},'File','temp');
            F = fileread('temp.m');
            i = strfind(F,'%');
            F = F(i(1):end);

            fileID = fopen('temp.m','w');
            f_mupad  = [fA_mupad;fB_mupad;fM_mupad;fF_mupad;fE_mupad];
            f_casadi = [fA_casadi;fB_casadi;fM_casadi;fF_casadi;fE_casadi];
            for i=1:length(f_mupad)
                fprintf(fileID,'%s=f_casadi{%d};',name(f_casadi{i}),i);
            end
            
            fwrite(fileID,F);
            fclose(fileID);
            rehash;

            in1 = [x(xr);z(zr)];
            in2 = u;
            in3 = p;
            U1 = u;
            P1 = p;
            W1 = w;
            S1 = s;
            in6 = w;
            temp;
            
            M = matlabFunction(M,'Vars',{XZ,U,P,t_mupad,Q,W,S},'File','temp');
            M = fileread('temp.m');
            i = strfind(M,'%');
            M = M(i(1):end);

            fileID = fopen('temp.m','w');
            for i=1:length(f_mupad)
                fprintf(fileID,'f%d=f_casadi{%d};',i,i);
            end
            
            fwrite(fileID,M);
            fclose(fileID);
            rehash;
            
            in1 = [x(xr);z(zr)];
            in2 = u;
            in3 = p;
            U1 = u;
            P1 = p;
            W1 = w;
            S1 = s;
            in6 = w;
            temp;

            R = matlabFunction(replaced_variables,'Vars',{XZ,U,P,t_mupad,Q,W,S},'File','temp');
            R = fileread('temp.m');
            i = strfind(R,'%');
            R = R(i(1):end);

            fileID = fopen('temp.m','w');
            for i=1:length(f_mupad)
                fprintf(fileID,'f%d=f_casadi{%d};',i,i);
            end
            
            fwrite(fileID,R);
            fclose(fileID);
            rehash;
            
            in1 = [x(xr);z(zr)];
            in2 = u;
            in3 = p;
            U1 = u;
            P1 = p;
            W1 = w;
            S1 = s;
            in6 = w;
            temp;
            
            if ~isempty(replacement_x)
              x_old(replacement_x) = replaced_variables(1:numel(replacement_x));
            end
            if ~isempty(replacement_z)
              z_old(replacement_z) = replaced_variables(numel(replacement_x)+1:end);
            end
            %Mcopy = SX(M);
            %Mcopy(1:length(xr),1:length(xr))=0;
            %assert(nnz(sparsify(Mcopy))==0,'Mass matrix has entries outside of 1:nxr')
            
            x = x(xr(:));
            dx = SX.sym('dx',nx);
            dx = dx(xr(:));
            z = z(zr(:));
            
            while true
 
                M = sparsify(SX(M));
            
                nx = numel(x);
                nz = numel(z);

                Mpat = patt(M);

                %
                assert(nnz(Mpat(nx+1:end,nx+1:end))==0)
                assert(nnz(Mpat(1:nx,nx+1:end))==0)

                if nnz(Mpat(nx+1:end,1:nx))~=0
                   warning('extra elimination needed'); 
                   
                   % Transformation: order rows such that
                   % equations without z go first
                   %
                   % Needed because we assume a partition of M
                   priority = find(~sum(patt(jacobian(F,z)),2));
                   priorityi = find(sum(patt(jacobian(F,z)),2));
                   M = [M(priority,:);M(priorityi,:)];
                   F = [F(priority,:);F(priorityi,:)];
                   
                   % Transformation: non-zero rows in algebraic part of
                   % mass matrix can (in fact must!) be eliminated
                   elr = find(sum(Mpat(nx+1:end,1:nx),2));
                   elri = find(~sum(Mpat(nx+1:end,1:nx),2));
                   

                   A = jacobian(F(nx+elr),[z]);
                   
                   elc = find(sum(patt(A),1))';
                   elci = find(~sum(patt(A),1))';
                   A = A(:,elc);
                   
                   assert(size(A,1)==size(A,2));
 
                   zsol = A\(substitute(-F(nx+elr),z(elc),0)+M(nx+elr,1:nx)*dx);
                   e = substitute(M*[dx;zeros(nz,1)]-F,z(elc),zsol);

                   M = jacobian(e,dx);
                   F = -substitute(e,dx,0);
                   
                   zsol_no_dx = substitute(zsol,dx,M(1:nx,1:nx)\F(1:nx));
                   
                   x_old = substitute(x_old,z(elc),zsol_no_dx);
                   z_old = substitute(z_old,z(elc),zsol_no_dx);

                   M = M([(1:nx)';(nx+elri)],:);
                   M = [M zeros(size(M,1),size(M,1)-size(M,2))];
                   F = F([(1:nx)';(nx+elri)],:);

                   z = z(elci);

                end

                Mcopy = SX(M);
                Mcopy(1:nx,1:nx)=0;
                if nnz(sparsify(Mcopy))==0
                    break
                end
            end
            
            Y = self.y([x_old;z_old],u,p,t);

            Fun = Function('F',{x,z,u,p,t,q,w,s},{M,F,Y},{'xr','zr','u','p','t','q','w','s'},{'M','F','y'});
              
            [~,zr] = find(patt(jacobian(z,z_orig)));
            [~,xr] = find(patt(jacobian(x,x_orig)));
            

            
        end
    end
    methods(Static)
        function out = collect_derivatives(e)
            s = char(e);
            matches = regexp(s,'D\(.*?\)','match');
            
            out = {};
 
            for i=1:numel(matches)
               m = matches{i};
               if m(3)=='['
                   m2 = regexp(m,'[\d+\]','match');
                   m2 = m2{1};
                   arg = str2num(m2(2:end-1));
                   
                   m3 = strsplit(m,',');
                   m3 = m3{2};
                   funname = m3(2:end-1);
               else
                   arg = 1;
                   funname = m(3:end-1);
               end
               
               % Why strings, because unqique works for string cell 
               out = {out{:} [funname ',' num2str(arg)]};
               
            end
            out = unique(out);
            
            for i=1:numel(out)
               out{i} = strsplit(out{i},',');
               out{i}{2} = str2num(out{i}{2});
            end
            
        end
        function [E_mupad, f_mupad, f_casadi] = to_mupad(E, label)
           import casadi.*
           args = symvar(E);
           v = vertcat(args{:});
           
           Enum = full(DM(E));
           
           
           
           s = DM(sparsity(jacobian(E,v)),1);
           N = size(s,1);
           fcount = 1;
           f_casadi = cell(N,1);
           for i=1:N
             deps = {};
             deps_names = {};
             deps_i = find(full(s(i,:)==1));
             for d=deps_i
                deps = {deps{:} v(d)};
                deps_names = {deps_names{:} name(v(d))};
             end
             if ~isempty(deps_i)
               f_casadi{fcount} = Function(['f_' label num2str(fcount)],deps,{E(i)},deps_names,{'out'});
               fcount = fcount+1;
             end
           end

           V_mupad = cell(N,1);
           for i=1:numel(v)
             V_mupad{i} = sym(name(v(i)));
           end

           E_mupad = sym(zeros(size(E)));

           fcount = 1;
           f_mupad = cell(N,1);
           for i=1:N
              deps_i = find(full(s(i,:)==1));
              deps = '';
              for d=deps_i
                  n = name(v(d));
                  if n(1)=='p' || n(1)=='t' || n(1)=='q' || n(1)=='s' 
                    deps = [deps  n ','];  
                  else
                    deps = [deps  n '(t),'];
                  end
              end
              
              if isempty(deps_i)
                E_mupad(i) = Enum(i);
              else
                f = sym(['f_' label num2str(fcount) '(' deps(1:end-1) ')']);
                f_mupad{fcount} = f;
                E_mupad(i) = f;
                fcount = fcount+1;
              end
            end
            f_mupad  = f_mupad(1:fcount-1);
            f_casadi = f_casadi(1:fcount-1);
        end
    end
    
end


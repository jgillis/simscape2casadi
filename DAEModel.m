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

          E = self.a(p) * [x;z] + self.b(p)*u + self.f([x;z],u,p,t);
          M = self.m(p);
          M = M(1:self.nx,1:self.nx);
          out = Function('E',{x,z,u,p,t},{M\E(1:self.nx),E(self.nx+1:end)},{'x','z','u','p','t'},{'ode','alg'});
        end
        function [out,xr,zr] = dae_r_expl(self)
          import casadi.*
          [Fun, xr, zr] = self.Fr();
          x = SX.sym('x',self.nx);
          z = SX.sym('z',self.nz);
          u = SX.sym('u',self.nu);
          p = SX.sym('p',self.np);
          t = SX.sym('t');
          
          [M,rhs] = Fun(x(xr),z(zr),u,p,t);
          nxr = numel(xr);
          M = M(1:nxr,1:nxr);

          res = rhs(nxr+1:end);
          rhs = rhs(1:nxr);
          
          out = Function('E',{x(xr),z(zr),u,p,t},{M\rhs,res},{'xr','zr','u','p','t'},{'ode','alg'});
        end
        function out = ode_expl(self)
          import casadi.*
          dae = self.dae_expl();
          
          dae_in = sx_in(dae);
          x = dae_in{1};
          z = dae_in{2};
          u = dae_in{3};
          p = dae_in{4};
          t = dae_in{5};
          [rhs,res] = dae(dae_in{:});
          
          % Check if res in linear in z
          assert(~any(cell2mat(which_depends(res, z, 2, false))))
          J = jacobian(res,z);
          
          res

          zsol = -J\substitute(res,z,0);
          rhs = substitute(rhs,z,zsol);
          
          out = Function('E',{x,u,p,t},{rhs},char('x','u','p','t'),char('rhs'));
        end
        function [out,xr,zr] = ode_r_expl(self)
          import casadi.*
          [dae,xr,zr] = self.dae_r_expl();
          
          dae_in = sx_in(dae);
          x = dae_in{1};
          z = dae_in{2};
          u = dae_in{3};
          p = dae_in{4};
          t = dae_in{5};
          
          [rhs,res] = dae(dae_in{:});
          
          % Check if res in linear in z
          assert(~any(cell2mat(which_depends(res, z, 2, false))))
          J = jacobian(res,z);

          zsol = -J\substitute(res,z,0);
          rhs = substitute(rhs,z,zsol);
          
          out = Function('E',{x,u,p,t},{rhs},{'x','u','p','t'},{'rhs'});
        end
        function [Fun, xr, zr] = Fr(self)
            model = Model;

            nx = self.nx;
            nz = self.nz;
            nu = self.nu;
            np = self.np;

            import casadi.*
            x = SX.sym('x',nx);
            z = SX.sym('z',nz);
            u = SX.sym('u',nu);
            p = SX.sym('p',np);
            t = SX.sym('t');
            
            pp = p;

            M = self.m(p);

            M = M(1:nx,1:nx);

            A = self.a(p);
            B = self.b(p);

            f_expr = self.f([x;z],u,p,t);

            
            t_mupad = sym('t');
            x_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(x),'uni',false);
            dx_mupad = cell(nx,1);
            for i=1:nx
               dx_mupad{i} = diff(x_mupad{i},t_mupad);
            end
            z_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(z),'uni',false);
            u_mupad = cellfun(@(e) sym([name(e) '(t)']),vertsplit(u),'uni',false);
            p_mupad = cellfun(@(e) sym(name(e)),vertsplit(p),'uni',false);

            [A_mupad, fA_mupad, fA_casadi] = DAEModel.to_mupad(A, 'A');
            [B_mupad, fB_mupad, fB_casadi] = DAEModel.to_mupad(B, 'B');
            [M_mupad, fM_mupad, fM_casadi] = DAEModel.to_mupad(M, 'M');
            [f_expr_mupad, fF_mupad, fF_casadi] = DAEModel.to_mupad(f_expr, 'F');
            
            vars = [vertcat(x_mupad{:});vertcat(z_mupad{:})];

            eqs = [M_mupad*vertcat(dx_mupad{:});zeros(nz,1)]==A_mupad*vars+B_mupad*vertcat(u_mupad{:})+f_expr_mupad;

            [newEqs,newVars,~] = reduceRedundancies(eqs,vars);
            size(newEqs);
            size(newVars);

            assert(isLowIndexDAE(newEqs,newVars))

            [M,F] = massMatrixForm(newEqs,newVars);

            X = sym('X',[nx,1]);
            Z = sym('Z',[nz,1]);
            U = sym('U',[nu,1]);
            P = sym('P',[np,1]);
            T = sym('T',[1,1]);
            
            F = subs(F,[vars;vertcat(u_mupad{:});vertcat(p_mupad{:})],[X;Z;U;P]);
            M = subs(M,[vertcat(p_mupad{:})],[P]); 

            newVars = subs(newVars,vars,[X;Z]);

            [~,xr] = find(jacobian(newVars,[X]));
            [~,zr] = find(jacobian(newVars,[Z]));

            F = matlabFunction(F,'Vars',{[X;Z],U,P,t_mupad},'File','temp');
            F = fileread('temp.m');
            i = strfind(F,'%');
            F = F(i(1):end);

            fileID = fopen('temp.m','w');
            f_mupad  = [fA_mupad;fB_mupad;fM_mupad;fF_mupad];
            f_casadi = [fA_casadi;fB_casadi;fM_casadi;fF_casadi];
            for i=1:length(f_mupad)
                fprintf(fileID,'%s=f_casadi{%d};',name(f_casadi{i}),i);
            end
            
            fwrite(fileID,F);
            fclose(fileID)
            rehash

            in1 = [x;z];
            in3 = p;
            U1 = u;
            temp
            
            M = matlabFunction(M,'Vars',{P},'File','temp');
            M = fileread('temp.m');
            i = strfind(M,'%');
            M = M(i(1):end);

            fileID = fopen('temp.m','w');
            for i=1:length(f_mupad)
                fprintf(fileID,'f%d=f_casadi{%d};',i,i);
            end
            
            fwrite(fileID,M);
            fclose(fileID)
            rehash
            
            in1 = p;
            temp

            Fun = Function('F',{x(xr),z(zr),u,p,t},{M,F},{'xr','zr','u','p','t'},{'M','F'});
        end
    end
    methods(Static)
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
                  if n(1)=='p' || n(1)=='t'
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


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
        function [Fun, xr, zr] = Fr(self)
            model = Model;

            nx = self.nx;
            nz = self.nz;
            nu = self.nu;

            import casadi.*
            x = SX.sym('x',nx);
            z = SX.sym('z',nz);
            u = SX.sym('u',nu);
            t = SX.sym('t');

            M = self.m;

            M = M(1:nx,1:nx);

            A = self.a;
            B = self.b;
            r = self.r';

            f_expr = self.f([x;z]);

            F_struct_x = DM(sparsity(jacobian(f_expr,x)),1);
            F_struct_z = DM(sparsity(jacobian(f_expr,z)),1);

            f_casadis = cell(nx+nz,1);
            for i=1:nx+nz
              deps_i = find(full(F_struct_x(i,:)==1));
              deps = {t};
              for d=deps_i
                  deps = {deps{:} x(d)};
              end
              deps_i = find(full(F_struct_z(i,:)==1));
              for d=deps_i
                  deps = {deps{:} z(d)};
              end
              f_casadis{i} = Function('f',deps,{f_expr(i)});
            end

            syms t
            x = cell(nx,1);
            dx = cell(nx,1);
            for i=1:nx
              x{i} = sym(['x' num2str(i) '(t)']);
              dx{i} = diff(x{i},t);
            end

            z = cell(nz,1);
            for i=1:nz
              z{i} = sym(['z' num2str(i) '(t)']);
            end

            u = cell(nu,1);
            for i=1:nu
              u{i} = sym(['u' num2str(i) '(t)']);
            end

            f = cell(nx+nz,1);
            for i=1:nx+nz
              deps_i = find(full(F_struct_x(i,:)==1));
              deps = 't';
              for d=deps_i
                  deps = [deps ',x' num2str(d) '(t)'];
              end
              deps_i = find(full(F_struct_z(i,:)==1));
              for d=deps_i
                  deps = [deps ',z' num2str(d) '(t)'];
              end
              if ~is_zero(f_expr(i))
                f{i} = sym(['f' num2str(i) '(' deps ')']);
              else
                f{i} = 0;
              end
            end


            A = full(model.a);
            B = full(model.b);
            M = full(model.m);

            vars = [vertcat(x{:});vertcat(z{:})];

            eqs = M*[vertcat(dx{:});zeros(nz,1)]==A*vars+B*vertcat(u{:})+vertcat(f{:});

            [newEqs,newVars,~] = reduceRedundancies(eqs,vars);
            size(newEqs);
            size(newVars);

            assert(isLowIndexDAE(newEqs,newVars))

            [M,F] = massMatrixForm(newEqs,newVars);

            Mf = matlabFunction(M);
            M = Mf();

            X = sym('x',[nx,1]);
            Z = sym('Z',[nz,1]);
            U = sym('U',[nu,1]);

            F = subs(F,[vars;vertcat(u{:})],[X;Z;U]);

            newVars = subs(newVars,vars,[X;Z]);

            [~,xr] = find(jacobian(newVars,[X]));
            [~,zr] = find(jacobian(newVars,[Z]));

            F = matlabFunction(F,'Vars',{[X;Z],U,t},'File','temp');
            F = fileread('temp.m');
            i = strfind(F,'%');
            F = F(i(1):end);

            fileID = fopen('temp.m','w');
            for i=1:nx+nz
              if ~is_zero(f_expr(i))
                fprintf(fileID,'f%d=f_casadis{%d};',i,i);
              end
            end
            
            fwrite(fileID,F);
            fclose(fileID)
            rehash

            x = SX.sym('x',model.nx);
            z = SX.sym('z',model.nz);
            u = SX.sym('u',model.nu);

            in1 = [x;z];
            U1 = u;
            t = 0;
            temp

            Fun = Function('F',{x(xr),z(zr),u},{M,F});
        end
    end
    
end


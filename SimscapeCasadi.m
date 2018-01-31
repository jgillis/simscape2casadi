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
classdef SimscapeCasadi
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
      function r = if_else(cond,iftrue,iffalse)
        if islogical(cond)
          cond = double(cond);  
        end
        r = if_else(casadi.SX(cond),iftrue,iffalse);
      end
      function r = tlu2_1d_linear_linear_value(in)
        r=SimscapeCasadi.pw_lin(in{3},in{1},in{2});
      end
      function r = tlu2_2d_linear_linear_value(in)
        % The SX implementation
        M = reshape(in{3},numel(in{1}),numel(in{2}));
        res = cell(1,numel(in{2}));
        for i=1:numel(in{2})
            res{i} = SimscapeCasadi.pw_lin(in{4},in{1},M(:,i));
        end
        res = [res{:}];
        r = SimscapeCasadi.pw_lin(in{5},in{2},res);
      end
    
      function y=pw_lin(t, tval, val)
        import casadi.*
        N = numel(tval);
        assert(N>=2);
        assert(numel(val)==N);
        tval = SX(tval);

        % Gradient for each line segment
        g = SX(1, N-1);
        for i=1:N-1
          g(i) = (val(i+1)- val(i))/(tval(i+1)-tval(i));
        end
  
        lseg = SX(1, N-1);
        for i=1:N-1
          lseg(i) = val(i) + g(i)*(t-tval(i));
        end
        y = pw_const(t, vec(tval(2:N-1))', vec(lseg)');
      end
  end
    
end


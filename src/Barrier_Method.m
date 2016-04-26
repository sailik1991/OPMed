% The Barrier Method for minimizing optimizations problems with only
% inequality constraints is as follows:
%
% min k * f_0(x) + sum_{i = 1 to m} \phi_i
%
% where, \phi_i is the i-th inequality constraint

% Example functions defined with symbols
% min f_0(x) = (1/x1 + 1/x2 + 1/x3)
% s.t.  x1 + x2 <= 2,
%       x1 + x3 <= 2,
%       x2 + x3 <= 2,
%       x1 >= 0,
%       x2 >= 0,
%       x3 >= 0,

syms x1 x2 x3 k

symbols = [x1 x2 x3 k];
f = k * (1/x1 + 1/x2 + 1/x3) - (log(2-x1-x2) + log(2-x2-x3) + log(2-x1-x3) + log(x1) + log(x2) + log(x3));

g = gradient(f, [x1 x2 x3]);
H = hessian(f, [x1 x2 x3]);

eps=1e-05;
a = 0.5;
b = 0.9;
mew = 12;

x_0 = [0.5 0.5 0.5];
k_0 = 1;
i=0;
while 1
    i = i + 1;
    
    %fprintf('Iteration no. (i) & t & Backtracking t & Norm & Step Size & x_i \\\\ \n')
    %fprintf('\\hline\n')
    
    % Solve for x(k)
    while 1
        g_x = eval(subs(g, [x1 x2 x3 k], [x_0 k_0]));
        h_x = eval(subs(H, [x1 x2 x3 k], [x_0 k_0]));
        dir = h_x \ g_x;
        lambda_2 = transpose(g_x) * dir; 
        dir = -1 * dir;

        if lambda_2/2 <= eps
            break 
        end
        t = arjimo(f, g, x_0, symbols, k_0, dir, a, b);
        x_0 = x_0 + t * transpose(dir);
        %fprintf('%d & %d & %0.6f & %0.6f & %0.6f & (%0.6f, %0.6f, %0.6f) \\\\ \n', j, k_0, t, double(norm(eval(subs(g, [x1 x2 x3 k], [x_0 k_0])))), double(t), double(x_0) );
    end
    
    %fprintf('\\hline\n')
    
    if 6/k_0 < eps
       break 
    end
    k_0 = mew*k_0;
end
disp(double(x_0))
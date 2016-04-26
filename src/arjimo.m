% This function implements Armijo's Backtracking line search
% @Inputs
%   f(x)
%   g(x) = f'(x)
%   x_0 point for calculating descent
%   symbols = list of symbols to represent f(x) in abstract format
%   values = values to be assigned to the symbols (other than 'x_0')
%   dir = Descent direction
%   a, b = Parameters for Arjimo's backtracking line search
%
%   @Returns
%       't' value for step size
function a_x = arjimo(f, g, x_0, symbols, other_values, dir, a, b)
    t = 1;
    val = [x_0 other_values];
    f_x = eval(subs(f, symbols, val));
    while 1
        x_new = x_0 + t * transpose(dir);
        val = [x_new other_values];
        f_x_new = eval(subs(f, symbols, val));
        val = [x_0 other_values];
        chng = a*t*transpose(eval(subs(g, symbols, val))) * dir;
        if isreal(f_x_new) == 0
            t = b * t;
            continue
        end
        if f_x_new <= f_x + chng
            break
        end
        t = b * t;
    end
    a_x = t;
end
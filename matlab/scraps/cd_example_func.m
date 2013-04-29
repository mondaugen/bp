function [ f, g ] = cd_example_func( x )
%CD_EXAMPLE_FUNC Trying out gradient descent optimixation
f = (3/2)*x(1)^2 + 2*x(2)^2 + (3/2)*x(3)^2 + x(1)*x(3) + 2*x(2)*x(3) - 3*x(1) - x(3);
g = [3*x(1) + x(3) - 3
     4*x(2) + 2*x(3)
     x(1) + 2*x(2) + 3*x(3) - 1];
end


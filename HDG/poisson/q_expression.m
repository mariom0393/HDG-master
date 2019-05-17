syms x y lambda 

u = 4*y^2-4*lambda^2*y*exp(-lambda*y)*cos(6*pi*x)+lambda*(-2*lambda*y);

delta_u = gradient(u, [x y]);
n = [0 -1];

q = -delta_u(2);

q_func = matlabFunction(q)
syms lambda x y
u = 4*y^2-4*lambda^2*y*exp(-lambda*y)*cos(6*pi*x)+lambda*exp(-2*lambda*y);

n = [0 -1];

u_xy = gradient(u,[x y]);

t = dot(u_xy,n)
t_func = matlabFunction(t)
syms x y k a b lambda

%Analytical symbolic expression
u = sinh(a*sin(k*pi*x)+b*cos(pi*(x+lambda*y)));
%u = k*log((x+a)^2)+cos(b*(y-lambda)^3); %agus
%u = cos(a*pi*(k*(x^2-b)))*sinh(-lambda*y); % eugenio

%Derived expressions
delta_u = gradient(u, [x, y]); 
s = - k*divergence(delta_u,[x, y]); %source

%Function u ---> Dirichlet 
FUNC_u = matlabFunction(u);
u_D = FUNC_u(1.2,0.75,5,3,x,y); %(a,b,k,lambda,x,y)

%Neumann (y=0)
n_N = [0 -1]; %check normal vector according to boundary
FUNC_t = matlabFunction(dot(k*delta_u,n_N));
t = FUNC_t(1.2,0.75,5,3,x,y);

%Robin (x=1)
n_R = [1 0];
FUNC_g = matlabFunction(simplify(dot(k*delta_u,n_R)+lambda*u));
g = FUNC_g(1.2,0.75,5,3,x,y);

LaTeX_Expr = latex(g)
clearvars
syms x y k a b gamma

%Analytical symbolic expression
u = sinh(a*sin(k*pi*x)+b*cos(pi*(x+gamma*y)));

%u = k*log((x+a)^2)+cos(b*(y-gamma)^3); %agus
%u = cos(a*pi*(k*(x^2-b)))*sinh(-gamma*y); % eugenio

%Derived expressions
delta_u = gradient(u, [x, y]); 
s = - k*divergence(delta_u,[x, y]); %source
s_F= matlabFunction(s);

%Function u ---> Dirichlet 
FUNC_u = matlabFunction(u);
u_D = FUNC_u(1.2,0.75,5,3,x,y); %(a,b,k,gamma,x,y)

x_vec0 = linspace(0,1,200);
y_vec0 = linspace(0,1,200);

[x_vec, y_vec] = meshgrid(x_vec0,y_vec0);

z_vec = FUNC_u(1.2,0.75,5,3,x_vec,y_vec);

s=surf(x_vec,y_vec,z_vec);
axis equal
colormap(jet)
s.EdgeColor = 'none';


%Neumann (y=0)
n_N = [0 -1]; %check normal vector according to boundary
FUNC_t = matlabFunction(dot(k*delta_u,n_N));


%Robin (x=1)
n_R = [1 0];
FUNC_g = matlabFunction(dot(k*delta_u,n_R)+gamma*u);
g = FUNC_g(1.2,0.75,3,5,x,y)

%LaTeX_Expr = latex(u_D)
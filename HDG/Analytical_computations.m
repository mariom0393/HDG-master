clearvars
syms x y k a b gamma

%Analytical symbolic expression
%u = sinh(a*sin(k*pi*x)+b*cos(pi*(x+gamma*y)));

%u = k*log((x+a)^2)+cos(b*(y-gamma)^3); %agus
u = cos(a*pi*(k*(x^2-b)))*sinh(-gamma*y); % eugenio
us = cos(4*pi*(2.5*(x^2-0.3)))*sinh(-1.1*y); % eugenio
um= matlabFunction(u);
latexum = latex(u);

%Derived expressions
delta_u = gradient(u, [x, y]); 
delta_us = gradient(us, [x, y]); 
%Source Term
s = - k*divergence(delta_u,[x, y]); 
sm= matlabFunction(s);
latexsm = latex(s);

%Function u ---> Dirichlet 
u_D = um(a,b,k,gamma,x,y);
u_Dx=um(4,0.3,2.5,1.1,0,y);%x=0
u_Dy=um(4,0.3,2.5,1.1,x,1);%y=1

%Plot of the function in 3D
x_vec0 = linspace(0,1,100);
y_vec0 = linspace(0,1,100);

[x_vec, y_vec] = meshgrid(x_vec0,y_vec0);

z_vec = um(4,0.3,2.5,1.1,x_vec,y_vec);

surf(x_vec,y_vec,z_vec);


%Neumann (x=1)
n_N = [1 0]; %check normal vector according to boundary
tm = matlabFunction(dot(k*delta_u,n_N));
tms = matlabFunction(dot(k*delta_us,n_N));
latextm= latex(dot(k*delta_u,n_N));
latextms= latex(dot(2.5*delta_us,n_N));

%Robin (y=0)
n_R = [0 -1];
FUNC_g = matlabFunction(simplify(dot(k*delta_u,n_R)+gamma*u));
gm = FUNC_g(a,b,k,gamma,x,y);
latexgm = latex(simplify(dot(k*delta_u,n_R)+gamma*u));
latexgms = latex(simplify(dot(2.5*delta_us,n_R)+1.1*us));
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
qxm=matlabFunction(delta_u(1));
qym=matlabFunction(delta_u(2));

%Plot of the function in 3D
x_vec0 = linspace(0,1,500);
y_vec0 = linspace(0,1,500);

[x_vec, y_vec] = meshgrid(x_vec0,y_vec0);

qx_vec = qxm(4,0.3,1.1,2.5,x_vec,y_vec);
qy_vec = qym(4,0.3,1.1,2.5,x_vec,y_vec);

figure(8)
sqx=surf(x_vec,y_vec,qx_vec);
%axis equal
colormap(jet)
sqx.EdgeColor = 'None';

figure(9)
sqy=surf(x_vec,y_vec,qy_vec);
%axis equal
colormap(jet)
sqy.EdgeColor = 'None';

%Source Term
s = - k*divergence(delta_u,[x, y]); 
sm= matlabFunction(s);
latexsm = latex(s);

%Function u ---> Dirichlet 
u_D = um(a,b,k,gamma,x,y);
u_Dx=um(4,0.3,2.5,1.1,0,y);%x=0
u_Dy=um(4,0.3,2.5,1.1,x,1);%y=1

z_vec = sm(4,0.3,1.1,2.5,x_vec,y_vec);
figure(10)
surf(x_vec,y_vec,z_vec);
s=surf(x_vec,y_vec,z_vec);
%axis equal
colormap(jet)
s.EdgeColor = 'None';


%Neumann (x=1)
n_N = [1 0]; %check normal vector according to boundary
tm = matlabFunction(dot(k*delta_u,n_N));
tms = matlabFunction(dot(k*delta_us,n_N));
latextm= latex(dot(k*delta_u,n_N));
latextms= latex(dot(2.5*delta_us,n_N));

z_vec = tm(4,0.3,1.1,2.5,x_vec,y_vec);
figure(11)
s2=surf(x_vec,y_vec,z_vec);
%axis equal
colormap(jet)
s2.EdgeColor = 'None';

%Robin (y=0)
n_R = [0 -1];
FUNC_g = matlabFunction(simplify(dot(k*delta_u,n_R)+gamma*u));
gm = FUNC_g(a,b,gamma,k,x,y);
latexgm = latex(simplify(dot(k*delta_u,n_R)+gamma*u));
latexgms = latex(simplify(dot(2.5*delta_us,n_R)+1.1*us));

z_vec = FUNC_g(4,0.3,1.1,2.5,x_vec,y_vec);
figure(12)
s3=surf(x_vec,y_vec,z_vec);
%axis equal
colormap(jet)
s3.EdgeColor = 'None';
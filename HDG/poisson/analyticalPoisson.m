function u = analyticalPoisson(X)

% Parameters
k = 5;
gamma = 3;
a=1.2;
b=0.75;

% Points
x = X(:,1);
y = X(:,2);

% Solution
u = sinh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y)));



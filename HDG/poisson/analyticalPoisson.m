function u = analyticalPoisson(X)

% Parameters
lambda = 4;
f = 6;

% Points
x = X(:,1);
y = X(:,2);

% Solution
u = 4*y.^2 - 4*lambda^2*y.*exp(-lambda*y).*cos(f*pi*x) + lambda*exp(-2*lambda*y);



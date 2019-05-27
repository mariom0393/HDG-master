function qy = analyticalPoissonqy(X)

% Parameters
k=2.5;
gamma=1.1;
a=4;
b=0.3;


% Points
x = X(:,1);
y = X(:,2);

% Solution
qy=-gamma.*cosh(gamma.*y).*cos(a.*k.*pi.*(b-x.^2));


function qx = analyticalPoissonqx(X)

% Parameters
k=2.5;
gamma=1.1;
a=4;
b=0.3;


% Points
x = X(:,1);
y = X(:,2);

% Solution
qx=a.*k.*x.*pi.*sinh(gamma.*y).*sin(a.*k.*pi.*(b-x.^2)).*-2.0;


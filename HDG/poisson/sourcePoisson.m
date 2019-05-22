function s = sourcePoisson(X, mu)

% Parameters
k=2.5;
gamma=1.1;
a=4;
b=0.3;

% Points
x = X(:,1);
y = X(:,2);

s = (1/k)*(k.*(gamma.^2.*sinh(gamma.*y).*cos(a.*k.*pi.*(b-x.^2))+a.*k.*pi.*sinh(gamma.*y).*sin(a.*k.*pi.*(b-x.^2)).*2.0-a.^2.*k.^2.*x.^2.*pi.^2.*sinh(gamma.*y).*cos(a.*k.*pi.*(b-x.^2)).*4.0));

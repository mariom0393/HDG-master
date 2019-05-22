function s = sourcePoisson(X, mu)

% Parameters
k = 5;
gamma = 3;
a=1.2;
b=0.75;

% Points
x = X(:,1);
y = X(:,2);

s = mu/k*(k.*(cosh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y))).*(b.*pi.^2.*cos(pi.*(x+gamma.*y))+a.*k.^2.*pi.^2.*sin(k.*x.*pi))-sinh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y))).*(b.*pi.*sin(pi.*(x+gamma.*y))-a.*k.*pi.*cos(k.*x.*pi)).^2-b.^2.*gamma.^2.*pi.^2.*sin(pi.*(x+gamma.*y)).^2.*sinh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y)))+b.*gamma.^2.*pi.^2.*cos(pi.*(x+gamma.*y)).*cosh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y)))));
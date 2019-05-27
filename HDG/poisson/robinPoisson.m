function g_R = robinPoisson(X)

% Parameters
k=2.5;
gamma=1.1;
a=4;
b=0.3;

% Points
x = X(:,1);
y = X(:,2);

%g_R=(1/k)*(-k.*sinh(k.*y).*cos(a.*gamma.*pi.*(b-x.^2))+cos(pi.*conj(a).*conj(gamma).*(conj(b)-conj(x).^2)).*cosh(conj(k).*conj(y)).*conj(gamma).*conj(k));
g_R=(1/k)*(-gamma.*sinh(gamma.*y).*cos(a.*k.*pi.*(b-x.^2))+cos(pi.*conj(a).*conj(k).*(conj(b)-conj(x).^2)).*cosh(conj(gamma).*conj(y)).*conj(gamma).*conj(k));
g_R=(1/k)*(-gamma.*sinh(gamma.*y).*cos(a.*k.*pi.*(b-x.^2))+cos(pi.*conj(a).*conj(k).*(conj(b)-conj(x).^2)).*cosh(conj(gamma).*conj(y)).*conj(gamma).*conj(k));

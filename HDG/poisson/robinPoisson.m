function q_R = robinPoisson(X)

% Parameters
k = 5;
gamma = 3;
a=1.2;
b=0.75;

% Points
x = X(:,1);
y = X(:,2);

q_R = (1/k)*(gamma.*sinh(a.*sin(k.*x.*pi)+b.*cos(pi.*(x+gamma.*y)))-cosh(cos(pi.*conj(x)+pi.*conj(gamma).*conj(y)).*conj(b)+conj(a).*sin(pi.*conj(k).*conj(x))).*conj(k).*(pi.*sin(pi.*conj(x)+pi.*conj(gamma).*conj(y)).*conj(b)-pi.*cos(pi.*conj(k).*conj(x)).*conj(a).*conj(k)));
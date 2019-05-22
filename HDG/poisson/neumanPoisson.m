function q_N = neumanPoisson(X)

% Parameters
k = 5;
gamma = 3;
a=1.2;
b=0.75;

% Points
x = X(:,1);
y = X(:,2);

q_N = (1/k)*(pi.*cosh(cos(pi.*(conj(x)+conj(gamma).*conj(y))).*conj(b)+conj(a).*sin(pi.*conj(k).*conj(x))).*sin(pi.*(conj(x)+conj(gamma).*conj(y))).*conj(b).*conj(k).*conj(gamma));


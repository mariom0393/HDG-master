function q_N = neumanPoisson(X)

% Parameters
k=2.5;
gamma=1.1;
a=4;
b=0.3;

% Points
x = X(:,1);
y = X(:,2);

q_N=(1/k)*(pi.*sin(pi.*conj(a).*conj(gamma).*(conj(b)-conj(x).^2)).*sinh(conj(k).*conj(y)).*conj(a).*conj(gamma).^2.*conj(x).*-2.0);

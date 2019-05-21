function q_N = neumanPoisson(X)

% Parameters
lambda = 4;

% Points
x = X(:,1);
y = X(:,2);

% Minus Laplacian
%q_N= y.*-8.0+lambda.^2.*2.0+lambda.^2.*exp(-lambda.*y).*cos(x.*pi.*6.0).*4.0...
    -lambda.^3.*y.*exp(-lambda.*y).*cos(x.*pi.*6.0).*4.0;

q_N = conj(y).*-8.0+exp(conj(lambda).*conj(y).*-2.0).*conj(lambda).^2.*2.0+exp(-conj(lambda).*conj(y)).*conj(lambda).^2.*cos(pi.*conj(x).*6.0).*4.0-...
    exp(-conj(lambda).*conj(y)).*conj(lambda).^3.*conj(y).*cos(pi.*conj(x).*6.0).*4.0;


% for i=1:length(q_N)
%     if y(i) == 0
%        q_N(i) = q_N(i);
%     else
%        q_N(i) = 0;
%     end
% end
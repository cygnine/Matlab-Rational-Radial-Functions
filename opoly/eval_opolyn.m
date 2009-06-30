function[p] = eval_opolyn(x,alpha,beta,n);

% function[p] = eval_opolyn(x,alpha,beta,n);
% Evaluates the normalized orthogonal polynomials defined by the recurrence
% coefficients alpha and beta. Assumes alpha and beta are long enough as
% necessary to evaluate the max(n)'th polynomial. 
% 
% Supports vectorization in n and x.
%
% Monic:
% p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}
% Normalized:
% sqrt(b_{n+1}) p_{n+1} = (x-a_n) p_n - sqrt(b_n) p_{n-1}

% 20080522: acn

% Pre-processing:
x = x(:);
n = n(:);
N = max(n);

p = zeros([length(x) N+1]);

p(:,1) = 1/sqrt(beta(1));
if N==0; 
  p = p(:,n+1);
  return;
end

p(:,2) = p(:,1).*(x-alpha(1));

for q=1:N;
  % Normalization of previous polynomial:
  p(:,q+1) = p(:,q+1)/sqrt(beta(q+1));
  % Computation of next orthogonal polynomial:
  p(:,q+2) = (x-alpha(q+1)).*p(:,q+1) - sqrt(beta(q+1))*p(:,q);
end

p(:,end) = p(:,end)/sqrt(beta(N+2));
p = p(:,n+1);

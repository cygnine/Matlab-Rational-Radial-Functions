function[x,w] = opoly_glq(as,bs,N,a,b);

% function[x,w] = opoly_glq(as,bs,N);
% Returns the N-point Gauss-Lobatto polynomial quadrature for the orthogonal
% polynomials corresponding to the recurrence coefficients as and bs.
% The fixed Lobatto points are a,b
%
% 20080522: acn

as = as(1:N);
bs = bs(1:N);

temp = eval_opoly([a; b],as,bs,[N-1,N-2]);
modif = inv(temp)*[a*temp(1,1); b*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
as(N) = modif(1);
bs(N) = modif(2);

as = as(:);
bs = bs(:);

J = spdiags([sqrt([bs(2:end);0]) as sqrt(bs)], -1:1, N, N);

[v,d] = eig(full(J));
x = diag(d);
w = bs(1)*v(1,:).^2.';

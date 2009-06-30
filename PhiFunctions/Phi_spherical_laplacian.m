function[laplacian,diff_mat] = Phi_spherical_laplacian(x,w,s,alpha,scale)

% Returns the nodal representation of the 2D spherical laplacian, f'' + 1/z*f'.
% Assumes the spectral representation is bijective. 
% In practice, this is really bad form...using collocation + weird
% multiplication by x^2...can code the right way with more time.

x = x(:);
w = w(:);
assert(size(x,1)==size(w,1), '(x,w) pair must be a valid quadrature rule');

N = size(x,1);
pl = Phi_eval(x,0:(N-1),s,alpha,scale);
dpl = dPhi_divided_x_eval(x,0:(N-1),s,alpha,scale);

diff_mat_x = dpl*(spdiags(w,0,N,N)*pl).';

% Use: d^2 f/dx^2 = 1/x*d/dx ( 1/x*df/dx ) * x^2 + 1/x*df/dx
% In practice, the *x.^2 part is well-conditioned: look at dr_dx_divided_x
laplacian = spdiags(x.^2,0,N,N)*diff_mat_x^2 + 2*diff_mat_x;

% Compute real differentiation matrix: is not that accurate for this basis set
dpl = dPhi_eval(x,0:(N-1),s,alpha,scale);
diff_mat = dpl*(spdiags(w,0,N,N)*pl).';

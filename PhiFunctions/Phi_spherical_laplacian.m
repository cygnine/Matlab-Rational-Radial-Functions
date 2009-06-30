function[laplacian,diff_mat] = PL_spherical_laplacian(x,w,s,alpha,scale)

% Returns the nodal representation of the 2D spherical laplacian, f'' + 1/z*f'.
% Assumes the spectral representation is bijective.

x = x(:);
w = w(:);
assert(size(x,1)==size(w,1), '(x,w) pair must be a valid quadrature rule');
assert(all(x~=0), 'If you give me a node at x=0, this construction will produce Infs');

N = size(x,1);
pl = Phi_eval(x,0:(N-1),s,alpha,scale);
dpl = dPhi_eval(x,0:(N-1),s,alpha,scale);

diff_mat = dpl*(spdiags(w,0,N,N)*pl).';

laplacian = diff_mat^2 + spdiags(1./x,0,N,N)*diff_mat;

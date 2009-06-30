%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute 2D radial rational chebshev Laplacian - Neumann bcs on r=[0,Lrh]
% Inputs: Lrh - length of domain, Nrh - number of mesh points
% Outputs: r - mesh, L - Laplacian, Dx - 1D differentiation matrix,
% w - integration weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rh,L,Dr,wh] = Compute_2D_radial_Laplacian_PL(Nrh,Lrh,s,alpha)

pl_parameters;

% s: default = 3/4, 
% noninclusive bounds: min: 1/2, max: infinity
% Tunable parameter: Chebyshev = 3/4
%                    Higher s -> less `emphasis' near x=infinity, but gives less
%                                numerical stability

% alpha: default = -1/2
% noninclusive bounds: min: -1, max: infinity
% Tunable parameter: Chebyshev = -1/2
%                    Higher alpha -> less `emphasis' near x=0, but gives less
%                                    numerical stability

% Default method: just let scale = Lrh
scale = Lrh;

% If you want Lrh to *really* be the physical nodal domain, uncomment out the
% following line of code.
% From N + length of domain, derive affine scaling parameter
% scale = derive_scale(Nrh,Lrh,s,alpha);

% Generate the stuff
[rh,wh] = gq(Nrh,s,alpha,scale);
[L,Dr] = PL_spherical_laplacian(rh,wh,s,alpha,scale);

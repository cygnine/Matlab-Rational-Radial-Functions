% script for debugging stuff

cd ..
N = 100;       % number of unknowns
scale = 4.0;   % scaling parameter
alpha = 3/4;
s = 1;
f = @(x) exp(-x.^2);
df = @(x) -2*x.*exp(-x.^2);
ddf = @(x) -2*exp(-x.^2) - 2*x.*-2.*x.*exp(-x.^2);

[x,laplacian,diff_mat,w] = Compute_2D_radial_Laplacian_Phi(N,scale,s,alpha);

lap_exact = ddf(x) + 1./x.*df(x);
lap_approx = laplacian*f(x);

fprintf('Unweighted Phi functions\n');
fprintf('Derivative error is %1.6e\n', norm(diff_mat*f(x)-df(x)));
fprintf('Laplacian error is %1.6e\n', norm(lap_approx-lap_exact));

[x,laplacian,diff_mat,w] = Compute_2D_radial_Laplacian_phi(N,scale,s,alpha);

lap_exact = ddf(x) + 1./x.*df(x);
lap_approx = laplacian*f(x);

fprintf('\n Weighted phi functions\n');
fprintf('Derivative error is %1.6e\n', norm(diff_mat*f(x)-df(x)));
fprintf('Laplacian error is %1.6e\n', norm(lap_approx-lap_exact));

cd debug

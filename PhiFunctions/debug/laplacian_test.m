% script for debugging stuff

cd ..
N = 100;       % number of unknowns
scale = 4.0;   % scaling parameter
alpha = -1/2;  
s = 1;
f = @(x) exp(-x.^2);
df = @(x) -2*x.*exp(-x.^2);
ddf = @(x) -2*exp(-x.^2) - 2*x.*-2.*x.*exp(-x.^2);

[x,w] = gq(N,s,alpha,scale);
[laplacian,diff_mat] = PL_spherical_laplacian(x,w,s,alpha,scale);

lap_exact = ddf(x) + 1./x.*df(x);
lap_approx = laplacian*f(x);

fprintf('Unweighted PL functions\n');
fprintf('Derivative error is %1.6e\n', norm(diff_mat*f(x)-df(x)));
fprintf('Laplacian error is %1.6e\n', norm(lap_approx-lap_exact));

[x,w] = pgq(N,s,alpha,scale);
[laplacian,diff_mat] = pl_spherical_laplacian(x,w,s,alpha,scale);

lap_exact = ddf(x) + 1./x.*df(x);
lap_approx = laplacian*f(x);

fprintf('\n Weighted pl functions\n');
fprintf('Derivative error is %1.6e\n', norm(diff_mat*f(x)-df(x)));
fprintf('Laplacian error is %1.6e\n', norm(lap_approx-lap_exact));

cd debug

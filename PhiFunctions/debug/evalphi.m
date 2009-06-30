% script for debugging stuff

cd ..
N = 100;
scale = 1.5;
alpha = -1/2;
s = 1;
f = @(x) exp(-(x).^2);
df = @(x) -2*x.*exp(-x.^2);
df_divided_x = @(x) -2.*exp(-(x).^2);

[x,w] = grq(N,s,alpha,scale);
pl = Phi_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('GRQ Mass matrix error is %1.6e\n', norm(mass(1:end-1,1:end-1)-eye(N-1)))

[x,w] = gq(N,s,alpha,scale);
pl = Phi_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('GQ Mass matrix error is %1.6e\n', norm(mass-eye(N)))

dpl = dPhi_eval(x,0:(N-1),s,alpha,scale);
dpl_x = dPhi_divided_x_eval(x,0:(N-1),s,alpha,scale);
modes = (spdiags(w,0,N,N)*pl)'*f(x);
d_reconstruction = dpl*modes;
fprintf('Derivative error at the nodes is %1.6e\n', norm(df(x) - ...
  d_reconstruction));
d_reconstruction_x = dpl_x*modes;
fprintf('Derivative divided by x error at the nodes is %1.6e\n', norm(df_divided_x(x) - ...
  d_reconstruction_x));

[x,w] = pgrq(N,s,alpha,scale);
pl = phi_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('PGRQ Mass matrix error is %1.6e\n', norm(mass(1:end-1,1:end-1)-eye(N-1)))

[x,w] = pgq(N,s,alpha,scale);
pl = phi_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('PGQ Mass matrix error is %1.6e\n', norm(mass-eye(N)))

dpl = dphi_eval(x,0:(N-1),s,alpha,scale);
dpl_x = dphi_divided_x_eval(x,0:(N-1),s,alpha,scale);
modes = (spdiags(w,0,N,N)*pl)'*f(x);
d_reconstruction = dpl*modes;
fprintf('Derivative error at the nodes is %1.6e\n', norm(df(x) - ...
  d_reconstruction));
d_reconstruction_x = dpl_x*modes;
fprintf('Derivative divided by x error at the nodes is %1.6e\n', norm(df_divided_x(x) - ...
  d_reconstruction_x));

cd debug

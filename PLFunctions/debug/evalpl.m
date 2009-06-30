% script for debugging stuff

cd ..
N = 100;
scale = 1.5;
alpha = -1/2;
s = 3/4;
f = @(x) exp(-(x-1).^2);
df = @(x) -2*(x-1).*exp(-(x-1).^2);

[x,w] = grq(N,s,alpha,scale);
pl = PL_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('GRQ Mass matrix error is %1.6e\n', norm(mass(1:end-1,1:end-1)-eye(N-1)))

[x,w] = gq(N,s,alpha,scale);
pl = PL_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('GQ Mass matrix error is %1.6e\n', norm(mass-eye(N)))

dpl = dPL_eval(x,0:(N-1),s,alpha,scale);
modes = (spdiags(w,0,N,N)*pl)'*f(x);
d_reconstruction = dpl*modes;
fprintf('Derivative error at the nodes is %1.6e\n', norm(df(x) - ...
  d_reconstruction));

[x,w] = pgrq(N,s,alpha,scale);
pl = pl_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('PGRQ Mass matrix error is %1.6e\n', norm(mass(1:end-1,1:end-1)-eye(N-1)))

[x,w] = pgq(N,s,alpha,scale);
pl = pl_eval(x,0:(N-1),s,alpha,scale);
mass = (spdiags(w,0,N,N)*pl).'*pl;
fprintf('PGQ Mass matrix error is %1.6e\n', norm(mass-eye(N)))

dpl = dpl_eval(x,0:(N-1),s,alpha,scale);
modes = (spdiags(w,0,N,N)*pl)'*f(x);
d_reconstruction = dpl*modes;
fprintf('Derivative error at the nodes is %1.6e\n', norm(df(x) - ...
  d_reconstruction));

cd debug

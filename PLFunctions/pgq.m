function[x,w] = pgq(N,s,alpha,scale)

pl_parameters;
beta = 2*s-2;

[x,w] = gq(N,s,alpha,scale); % Takes care of scale
r = x_to_r(x,scale);

weight = 1/2*(1+r).^(2*s);

x = r_to_x(r,scale);
w = w./weight;

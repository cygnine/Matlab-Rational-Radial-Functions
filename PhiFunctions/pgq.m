function[x,w] = pgq(N,s,alpha,scale)

phi_parameters;
beta = s-3/2;

[x,w] = gq(N,s,alpha,scale); % Takes care of scale
r = x_to_r(x,scale);

weight = (1+r).^(s);

x = r_to_x(r,scale);
w = w./weight;

function[weight] = weight_PL(x,alpha,s,scale)

pl_parameters;
r = x_to_r(x,scale);

weight = (1-r).^(alpha+1/2).*(1+r).^s;

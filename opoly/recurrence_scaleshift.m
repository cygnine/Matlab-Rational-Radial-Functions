% <<< jacobi_recurrence
% <<< hermite_recurrence
% <<< laguerre_recurrence
%
% Script for augmenting recurrence parameters to take into account affine
% mappings: shift and scale
% shift: new coordinate when x = 0
% scale: (new coordinate) - shift, when x = 1
%
% x --> w = c*x+d
% a_n <-- d + c*a_n
% b_n <-- c^2*b_n
%
% 20080623 -- acn

as = shift+scale*as;
bs = bs*scale^2;

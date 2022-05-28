function [trat] = tempratoblique(mn1,gamma);
%%%%%%%%
% this function computes the static temperature ratio across an oblique shock
%%%%%%%%%

trat = (1+2*gamma/(gamma+1)*(mn1^2-1))*(2+(gamma-1)*mn1^2)/((gamma+1)*mn1^2);
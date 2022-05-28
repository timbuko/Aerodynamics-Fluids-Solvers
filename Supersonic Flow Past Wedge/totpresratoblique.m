function [p_totalrat] = totpresratoblique(mn1,gamma);
%%%%%%%%
% this function computes the static pressure ratio across an oblique shock
%%%%%%%%%

p_totalrat = ((gamma+1)*mn1^2/((gamma-1)*mn1^2)+2)^(gamma/(gamma-1))*...
                    ((gamma+1)/(2*gamma*mn1^2-(gamma-1)))^(1/(gamma-1));
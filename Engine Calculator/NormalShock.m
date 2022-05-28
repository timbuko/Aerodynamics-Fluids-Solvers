function [mn2] = NormalShock(mn1,gamma);
% this function computes the normal component of the mach number
% down stream of an oblique shock
%%%%%%%%%%%%%%%%%%

mn2 = sqrt((1 + (gamma-1)/2*mn1^2)/(gamma*mn1^2-(gamma-1)/2));


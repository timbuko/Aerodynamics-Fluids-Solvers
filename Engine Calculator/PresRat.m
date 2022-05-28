function [prat] = PresRat(mn1,gamma);
%%%%%%%%
% this function computes the static pressure ratio across an oblique shock
%%%%%%%%%

prat = 1 + 2*gamma/(gamma + 1) *(mn1^2 -1);
function [prat] = presratoblique(mn1,gamma);
%%%%%%%%
% this function computes the static pressure ratio across an oblique shock
% p2/p2
%%%%%%%%%

prat = 1 + 2*gamma/(gamma + 1) *(mn1^2 -1);
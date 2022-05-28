function [prat] = presratpm(m1,m2,gamma);
%%%%%%%
%  this function computes the pressure ratio across an expansion wave
%%%%%%
    numer = 1 + (gamma-1)/2*m2^2;
    denom = 1  + (gamma-1)/2*m1^2;
    exper = gamma/(gamma-1);
    prat = (numer/denom)^exper;
    prat = 1.0/prat;
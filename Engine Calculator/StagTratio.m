function T0_T = StagTratio(M,p,gamma)
%Calculate the stagnation temperature to static temperature ratio

T0_T = 1+((gamma-1)*M^2)/2;

function p0_p = StagPratio(M,p,gamma)
%Calculate stagnation pressure to static pressure ratio
if nargin < 3
    gamma=1.4;
end

p0_p = (1+((gamma-1)*M^2)/2)^(gamma/(gamma-1));

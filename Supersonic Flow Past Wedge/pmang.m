function pmang = pmang(min,gamma);
%%%%
% this function computes the prandtl-meyer function nu(m) in radians
%%%%%
min2 = double(min);
val = sqrt((gamma-1)/(gamma+1)*(min2^2-1));
pmang = sqrt((gamma+1)/(gamma-1))*atan(val) - atan(sqrt(min2^2 -1));
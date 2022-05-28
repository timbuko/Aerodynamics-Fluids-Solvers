function [mafter] = findm2_new(nu,gamma);
%%%%%%%%
%  nu must be in radians!!!
% finds the mach number associated with a given prandtl meyer angle
%%%%%%%

syms  m   
	f = sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma -1)/(gamma+1)*(m^2-1)))-atan(sqrt(m^2-1)) == nu;
        ans1 =	solve(f,m);
	mafter = abs(ans1);
	

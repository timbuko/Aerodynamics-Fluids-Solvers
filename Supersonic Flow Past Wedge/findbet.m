function [beta] = findbet(m1,thet,gamma);
%%%%%%%%%%%%
%  calculate the shock angle beta given the upstream mach number and the 
% turning angle thet -- computes the weak shock solution
%angle in radians!!!!!
%%%%%%%%%%%%%

    
  lam = sqrt( ( m1^2 -1)^2 - 3*(1 + (gamma-1)/2*m1^2)*(1 + ...
       (gamma+1)/2*m1^2)*(tan(thet))^2);
    
  chi = (  (m1^2 - 1)^3 - 9*(1+ (gamma-1)/2*m1^2)*(1 + (gamma-1)/2*m1^2 + ...
       (gamma +1)/4*m1^4)*(tan(thet))^2)/lam^3;
  delta = 1;
  beta = atan( (m1^2 - 1 + 2 * lam*cos((4*pi*delta + acos(chi))/3))/...
    (3 * (1 + (gamma-1)/2*m1^2)*tan(thet)));

function [mafter] = findm2pm_Generic(nu,gamma);
%%%%%%%%
%  nu must be in radians!!!
% finds the mach number associated with a given prandtl meyer angle
%       ONLY VALID FOR 
%  Modified by Allan Bonhomme-Isaiah
%%%%%%%
      lmda=(gamma-1)/(gamma+1)
      TEMP1=strcat('(m^2-1)*',num2str(lmda));
      TEMP2=strcat('sqrt(',TEMP1);
      TEMP3=strcat(TEMP2,')');
      TEMP4=strcat('atan(',TEMP3);
      TEMP5=strcat(TEMP4,')*');
      TEMP6=strcat('sqrt(',num2str(1/lmda));
      TEMP7=strcat(TEMP6,')');
      TEMP8=strcat(TEMP5,TEMP7);
      TEMP9 = strcat(TEMP8,'-atan(sqrt(m^2-1))-');
      temp=strcat(TEMP9,num2str(nu));
      m= abs(solve(temp,'m'));
      mafter=m(1);

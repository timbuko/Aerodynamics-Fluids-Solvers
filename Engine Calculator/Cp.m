function [Cp,gamma]= Cp(Temperature,R)
%Use corelation fit to calc Cp (J/kgK) and gamma (Metric)
%Approximate T0~T 
if nargin < 2
    R=287;
end
Cp=1000*.9584*exp(0.000168*Temperature); 
gamma=Cp/(Cp-R);
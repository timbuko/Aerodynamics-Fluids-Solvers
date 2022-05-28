function [P03,T03,Cp_c,gamma_c,Tau_c]=compressor(P02,T02,pi_c,eta_c)
%Calculate flow properties in compressor. Iterate through gamma until
%consistent T03

convergence=0.001; %Degree K
T=T02;
error=1; %Difference in T03 from previous calc
T03=0;
while abs(error)>convergence
T03_i=T03;
[Cp_c,gamma_c] = Cp(T);
Tau_c= 1+(pi_c^((gamma_c-1)/gamma_c)-1)/eta_c;
T03= T02*Tau_c;
T= (T03+T02)/2;
error=T03-T03_i;
end

P03 = P02*pi_c;
function [P05,T05,Cp_t,gamma_t]=turbine(P04,T04,f,eta_m,eta_t,Cp_c,T03,T02);
%Calculate flow properties in turbine. Iterate Cp_t to find T05

convergence=0.001; %Degree K
T=T04;
error=1; %Difference in T03 from previous calc
T05=0;

while abs(error)>convergence
T05_i=T05;
[Cp_t,gamma_t] = Cp(T);
T05=T04-Cp_c*(T03-T02)/(eta_m*(1+f)*Cp_t);
T= (T05+T04)/2;
error=T05-T05_i;
end

Tau_t=T05/T04;
pi_t=(1-(1-Tau_t)/eta_t)^(gamma_t/(gamma_t-1));
P05 = P04*pi_t;
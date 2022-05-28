% Solve for all sections of a jet engine. Finds specific thrust and TSFC
% and max engine temp

clear all
close all
%% Effect of compression ratio and Mach number

NozzleType=2;      % 1 is design convdiv, 2 is choked conv only

P=[101325 53331 26500 12112];
T=[288.16 255.04 223 216.66];

for i=1:length(P)


M_a   =   0.8     ;%        Ambient Mach #
P_a   =   P(i)    ;%Pa      Ambient Pressure
T_a   =   T(i)    ;%K       Ambient Temperature
T04   =   1400    ;%K       Max Engine Temp
DH    =   4.14E7  ;%J/kg   Energy Density of Fuel

R     =   287     ;%J/kgK  Gas Constant
Cp_a  =   Cp(i)   ;%J/kgK  Specific heat of air
gamma_a=  Cp_a/(Cp_a-R)  ;%        Specific heat ratio of air

pi_d  =   0.92    ;%        Diffuser pressure ratio
pi_c  =   12      ;%        Compressor pressure ratio
pi_b  =   0.95    ;%        Burner pressure ratio
pi_n  =   1       ;%        Nozzle pressure ratio

eta_c =   0.90    ;%        Compressor efficiency factor
eta_b =   0.91    ;%        Burner efficiency factor
eta_m =   0.98    ;%        Mechanical efficiency factor (turbine to comp)
eta_t =   0.85    ;%        Turbine efficiency factor
eta_n =   0.96    ;%        Nozzle efficiency factor

%% Ambient
P0_a=P_a*StagPratio(M_a,P_a,gamma_a);
T0_a=T_a*StagTratio(M_a,P_a,gamma_a);

%% Diffuser
P02=inlet(P0_a,pi_d);
T02=T0_a;

%% Compressor
[P03,T03,Cp_c,gamma_c]=compressor(P02,T02,pi_c,eta_c);

%% Burner

[P04,T04,f,Cp_b]= burner(P03,T03,pi_b,eta_b,DH,gamma_a,T04);


%% Turbine

[P05,T05,Cp_t,gamma_t]=turbine(P04,T04,f,eta_m,eta_t,Cp_c,T03,T02);

%% Nozzle
T06=T05;P06=P05;
[P_e,T_e,v_e,Cp_n]=nozzle(P06,T06,P_a,eta_n,pi_n,R,NozzleType);
T08=T06;

%% Calc Specific Thrust and TSFC
v_a=M_a*sqrt(gamma_a*R*T_a);

TO4(i)=T04;
F_ma(i)=(1+f)*v_e-v_a+(P_e-P_a)*(1+f)*R*T_e/(P_e*v_e);
TSFC(i)=f./F_ma(i);
end


F_ma,TSFC

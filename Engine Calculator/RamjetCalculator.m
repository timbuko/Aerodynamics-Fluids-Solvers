% Solve for all sections of a jet engine. Finds specific thrust and TSFC

%% Input
NozzleType=2;      % 1 is design convdiv, 2 is choked conv only

M_a   =   2.0     ;%        Ambient Mach #
P_a   =   101325  ;%Pa      Ambient Pressure
T_a   =   288.2   ;%K       Ambient Temperature
T04   =   1666.5  ;%K       Max Engine Temp
DH    =   4.3E7   ;%kJ/kg   Energy Density of Fuel

Cp_a  =   1004    ;%kJ/kgK  Specific heat of air
gamma_a=  1.4     ;%        Specific heat ratio of air
R     =   287     ;%J/kgK  Gas Constant

pi_d  =   0.8     ;%        Diffuser pressure ratio
pi_b  =   0.95    ;%        Burner pressure ratio
pi_n  =   1       ;%        Nozzle pressure ratio
eta_b =   0.91    ;%        Burner efficiency factor
eta_n =   0.92    ;%        Nozzle efficiency factor

%% Calculate

%Ambient
P0_a=P_a*StagPratio(M_a,P_a,gamma_a);
T0_a=T_a*StagTratio(M_a,P_a,gamma_a);

%Diffuser
P02=inlet(P0_a,pi_d);
T02=T0_a;

%Burner
P03=P02;T03=T02;
[P04,temp]= burner(P03,T03,pi_b,eta_b,Cp_a,DH,gamma_a,T04);
if temp<1
    f=temp;
else
    T04=temp;
end

%Nozzle
T06=T04;P06=P04;
[P_e,T_e,v_e]=nozzle(P06,T06,P_a,eta_n,pi_n,Cp_a,gamma_a,R,NozzleType);
T08=T06;

%Calc Specific Thrust and TSFC
v_a=M_a*sqrt(gamma_a*R*T_a);

F_ma=(1+f)*v_e-v_a+(P_e-P_a)*(1+f)*R*T_e/(P_e*v_e);
TSFC=f/F_ma;


%Output the data
location={'Ambient','2','3','4','6','8'};
t=[T0_a,T02,T03,T04,T06,T08];
p=[P0_a,P02,P03,P04,P06,P06];
headers={'T0','P0'};
Ramjet=table(t',p');
Ramjet.Properties.VariableNames=headers;
Ramjet.Properties.RowNames=location

disp('Converging only (choked)')

f,F_ma,TSFC,P_e,v_e

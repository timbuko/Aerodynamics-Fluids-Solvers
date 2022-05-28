% Solve for all sections of a jet engine. Finds specific thrust and TSFC

%% Input
NozzleType=2;      % 1 is design convdiv, 2 is choked conv only

M_a   =   .6     ;%        Ambient Mach #
P_a   =   35600   ;%Pa      Ambient Pressure
T_a   =   216.7   ;%K       Ambient Temperature
T04   =   1700    ;%K       Max Engine Temp
DH    =   4.5E7   ;%kJ/kg   Energy Density of Fuel

Cp_a  =   1004    ;%kJ/kgK  Specific heat of air
gamma_a=  1.4     ;%        Specific heat ratio of air
R     =   287     ;%J/kgK  Gas Constant

pi_d  =   0.97     ;%        Diffuser pressure ratio
pi_c  =   30      ;%        Compressor pressure ratio
pi_b  =   0.96    ;%        Burner pressure ratio
pi_n  =   1       ;%        Nozzle pressure ratio

eta_c =   0.85    ;%        Compressor efficiency factor
eta_b =   0.96    ;%        Burner efficiency factor
eta_m =   0.99    ;%        Mechanical efficiency factor (turbine to comp)
eta_t =   0.90    ;%        Turbine efficiency factor
eta_n =   0.98    ;%        Nozzle efficiency factor

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

F_ma=(1+f)*v_e-v_a+(P_e-P_a)*(1+f)*R*T_e/(P_e*v_e);
TSFC=f/F_ma;


%% Output the data
% location={'Ambient','2','3','4','5','6','8'};
% t=[T0_a,T02,T03,T04,T05,T06,T08];
% p=[P0_a,P02,P03,P04,P05,P06,P06];
% headers={'T0','P0'};
% Ramjet=table(t',p');
% Ramjet.Properties.VariableNames=headers;
% Ramjet.Properties.RowNames=location
% 
% disp('Converging only (choked)')
% 
% f,F_ma,TSFC,P_e,v_e

fprintf('T_02 = %4.4f\n',T02)
fprintf('P_02 = %4.4f\n',P02)
fprintf('gamma_a = %4.4f\n',gamma_a)
fprintf('T_03 = %4.4f\n',T03)
fprintf('P_03 = %4.4f\n',P03)
fprintf('Cp_c = %4.4f\n',Cp_c)
fprintf('gamma_c = %4.4f\n',gamma_c)
fprintf('f = %4.4f\n',f)
fprintf('T_04 = %4.4f\n',T04)
fprintf('P_04 = %4.4f\n',P04)
fprintf('Cp_b = %4.4f\n',Cp_b)
fprintf('T_05 = %4.4f\n',T05)
fprintf('P_05 = %4.4f\n',P05)
fprintf('Cp_t = %4.4f\n',Cp_t)
fprintf('gamma_t = %4.4f\n',gamma_t)
fprintf('Cp_n = %4.4f\n',Cp_n)
fprintf('gamma_n = %4.4f\n',Cp_n/(Cp_n-R))
fprintf('P* = %4.4f\n',P_e)
fprintf('u8 = %4.4f\n',v_e)
fprintf('Area/m_dot = %4.4f\n',(1+f)/(P_e*v_e/T_e/R))
fprintf('ua = %4.4f\n',v_a)
fprintf('F/m_dot = %4.4f\n',F_ma)
fprintf('TSFC = %4.4e\n',TSFC)
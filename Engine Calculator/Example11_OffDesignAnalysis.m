% Solve for all sections of a jet engine. Finds specific thrust and TSFC
clear all
close all

%% Input

%Assumed mc2 and Nc2
% mc2=59.9;
% Nc2=9200;
mc20=60.381; 
Nc20=9193;

mc2=mc20;
Nc2=Nc20;

% %iterate mc2 to get pi_t_map=pi_t_calc
% error_pi_t=1;
% cnt=0;
% while(abs(error_pi_t)>.01)
% cnt=cnt+1
% if cnt>300
%     break
% end
% mc2=mc20+1-.01*cnt;

%%iterate Nc2 to get mc5_calc=mc5_map 
% error_mc5=1;
% cnt=0;
% while(abs(error_mc5)>.01)
% cnt=cnt+1;
% if cnt>300
%     break
% end
% Nc2=Nc20+2000-10*cnt;


NozzleType=2;      % 1 is design convdiv, 2 is choked conv only

M_a   =   .75     ;%        Ambient Mach #
P_a   =   101325   ;%Pa      Ambient Pressure
T_a   =   288.17  ;%K       Ambient Temperature
f     =   0.021    ;%        Fuel to air ratio
DH    =   1.78E7*2.326   ;%J/kg   Energy Density of Fuel

Cp_a  =   1004    ;%kJ/kgK  Specific heat of air
gamma_a=  1.4     ;%        Specific heat ratio of air
R     =   287     ;%J/kgK  Gas Constant

pi_d  =   0.92    ;%        Diffuser pressure ratio
pi_c  =   0       ;%        Compressor pressure ratio
pi_b  =   0       ;%        Burner pressure ratio
pi_t =    0       ;%        Turbine pressure ratio
pi_n  =   1       ;%        Nozzle pressure ratio

eta_c =   0    ;%        Compressor efficiency factor
eta_b =   0.91    ;%        Burner efficiency factor
eta_m =   0.995    ;%        Mechanical efficiency factor (turbine to comp)
eta_t =   0    ;%        Turbine efficiency factor
eta_n =   0.96    ;%        Nozzle efficiency factor

%zero means we have to solve for it later

%% Ambient
P0_a=146928; %Pa
T0_a=320.56; %Pa

%% Diffuser
P02=135173; %Pa
T02=320.56; %Pa

%% Compressor
[pi_c,eta_c]=comp_mapkg_11_1(mc2,Nc2);
while isnan(pi_c)|isnan(eta_c)
    i=i+1;
    mc2=mc2-.005;[pi_c,eta_c]=comp_mapkg_11_1(mc2,Nc2);
    if i>10000
        break
    end
end
[P03,T03,Cp_c,gamma_c]=compressor(P02,T02,pi_c,eta_c);

%% Burner
mc3=mc2*sqrt(T03/T02)/pi_c;
f_O=f*T_a/T03;
pi_b=burner_mapkg_11_1(mc3,f_O);

while isnan(pi_b)|isnan(eta_b)
    i=i+1;
    mc3=mc3-.005;[pi_b]=burner_mapkg_11_1(mc3,f_O);
    if i>10000
        break
    end
end

[P04,T04,f,Cp_b]= burner(P03,T03,pi_b,eta_b,DH,gamma_a,f);

%% Turbine
mc4=(1+f)*mc2*sqrt(T04/T03*T03/T02)/(pi_c*pi_b);
Nc4=Nc2/sqrt(T04/T03*T03/T02);
[pi_t,eta_t]=turb_mapkg_11_1(mc4,Nc4);
mc40=mc4;
i=0;
while isnan(pi_t)|isnan(eta_t)
    i=i+1
    if i<=10000
        mc4=mc40+.005*i;[pi_t,eta_t]=turb_mapkg_11_1(mc4,Nc4);
    elseif i>10000
        mc4=mc40-.005*(i-10000);[pi_t,eta_t]=turb_mapkg_11_1(mc4,Nc4);
    else 
       break
    end
    
end
    

[P05,T05,Cp_t,gamma_t]=turbine(P04,T04,f,eta_m,eta_t,Cp_c,T03,T02);

pi_t_map=pi_t
pi_t_calc=P05/P04

%Loop to find mc2
error_pi_t=pi_t_map-pi_t_calc;

%% Nozzle


T06=T05;P06=P05;
[P_e,T_e,v_e,Cp_n]=nozzle(P06,T06,P_a,eta_n,pi_n,R,NozzleType);
T08=T06;

%Iterate Nc2 so mc5 is equal
mc5_calc=8.8*sqrt(T05/T04)/pi_t
mc5_map=32.5;%Assume P06/P_a >1.7 --> mc5=32.5

error_mc5=mc5_calc-mc5_map;
%% while end
% end
%% Calc Specific Thrust and TSFC
v_a=M_a*sqrt(gamma_a*R*T_a);

F_ma=(1+f)*v_e-v_a+(P_e-P_a)*(1+f)*R*T_e/(P_e*v_e);
TSFC=f/F_ma;

Thrust=F_ma*32.5;
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

fprintf('Mach = %4.4f\n',0.75)
fprintf('f = %4.4f\n',0.021)
fprintf('mc2 = %4.4f\n',mc2)
fprintf('Nc2 = %4.4f\n',Nc2)
fprintf('eta_c = %4.4f\n',eta_c)
fprintf('gamma_c = %4.4f\n',gamma_c)
fprintf('Cp_c = %4.4f\n',Cp_c)
fprintf('mc3 = %4.4f\n',mc3)
fprintf('f_c = %4.4f\n',f_O)
fprintf('pi_b = %4.4f\n',pi_b)
fprintf('Cp_b = %4.4f\n',Cp_b)
fprintf('Tau_b = %4.4f\n',T04/T03)
fprintf('Cp_t = %4.4f\n',Cp_t)
fprintf('T04 = %4.4f\n',T04)
fprintf('mc4 = %4.4f\n',mc4)
fprintf('Nc4 = %4.4f\n',Nc4)
fprintf('eta_t = %4.4f\n',eta_t)
fprintf('Cp_t = %4.4f\n',Cp_t)
fprintf('pi_t map = %4.4f\n',pi_t_map)
fprintf('pi_t_calc = %4.4f\n',pi_t_calc)
fprintf('pi_t_calc = %4.4f\n',pi_t_calc)
fprintf('F/m_dot = %4.4f\n',F_ma)
fprintf('TSFC = %4.4e\n',TSFC)
fprintf('Thrust = %4.4f\n',Thrust)
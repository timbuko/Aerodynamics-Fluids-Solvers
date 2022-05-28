function [P_e,T_e,v_e,Cp_n] = nozzle(P06,T06,P_a,eta_n,pi_n,R,NozzleType)
%Calculate flow properties in nozzle. NozzleType=1=Design ConvDiv
%NozzleType=2=Choked Conv only

[Cp_n,gamma]=Cp(T06);

if NozzleType==1
    syms x
    P_e=P_a; 
    v_e=sqrt(2*Cp_n*eta_n*T06*(1-(P_e/P06)^((gamma-1)/gamma)));
    T_e=eval(solve((T06/x)==1+((gamma-1)*(v_e*sqrt(gamma*R*x))^2)/2,x));
    T_e=abs(T_e(1));   
elseif NozzleType==2
    
    P_star=P06*(1+(1-gamma)/(eta_n*(1+gamma)))^(gamma/(gamma-1));
    T_star=T06*2/(gamma+1);
   
    P_e=P_star; %Flow is choked
    T_e=T_star;
    v_e=sqrt(gamma*R*T_star);
end
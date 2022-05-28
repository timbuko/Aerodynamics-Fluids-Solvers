function [P04,T04,f,Cp_b] = burner(P03,T03,pi_b,eta_b,DH,gamma_a,T04orf)
%Calculate flow properties in burner. Can input T04 or f. Assume constant Cp
if T04orf<1
    f=T04orf;
    convergence=0.001;
    error=1;
    T04=0;
    T=T03;
    syms x
    while abs(error)>convergence
    T04_i=T04;
    [Cp_b,gamma_b]=Cp(T);
    T04=(eta_b*f*DH+Cp_b*T03)/((1+f)*Cp_b);
    T=(T04+T03)/2;
    error=T04_i-T04;
    end
    T04orf=T04;
    
else
    T04=T04orf;
    [Cp_b,gamma_b]=Cp((T03+T04)/2);
    
    f= Cp_b*(T04-T03)/(eta_b*DH-Cp_b*T04);
    
    T04orf=f;
end

P04=pi_b*P03;
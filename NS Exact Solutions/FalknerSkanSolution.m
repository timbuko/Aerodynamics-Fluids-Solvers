%Calculate Boundary Layer Profile for Falkner-Skan Flows of various Beta
%using shooting method and 4th order runga kutta

% clear all; close all; clc;
a = 0;
b = 10; 
h = .05;%Step size
N = ((b-a)/h); %Number of steps?
n(1) = 0;         %Declares initial n (eta) at 0
y(1) = 0;         %Boundary Condition G(0) = 0
z(1) = 0;         %Boundary Condition dG/dn(0) = 0


% %guess for f''(0) Boundary Condition ; Compare to f'(inf) = 1 for accuracy
% aBC = [4.5475e-16 4.5475e-15 0.02, 0.1286362263, 0.3384520231, 0.4696,...
%     0.774754584147709, 1.23258759711676, 1.47722388694746, 1.52151375980262];  
% aBCsupper=[4.5476e-15 4.5474e-15 1 .129 .3385 .4697 .775 1.233 1.4775, 1.522];
% aBCslower=[0 0 0 .128 .3384 .469 .773 1.232 1.477, 1.52] ;
% 
% B = [-1.9 -0.5 -0.19884 -0.18 -0.089 0.000 0.3 1.000 1.5, 1.6];   %Defines arrary for freestream pressure gradient (Beta)
j = 1;
% 

aBC = [ 4.5475e-15 0.4696 1 1.52151375980262];  
aBCsupper=[4.5475e-10 .4697 1.5 1.522];
aBCslower=[0 .469  .5 1.52] ;
B = [ -0.2 0.000 .5 1.6];   %Defines arrary for freestream pressure gradient (Beta)

% B=-1.9;
% aBC=4.5475e-16;
% aBCsupper=4.5476e-15;
% aBCslower=0;

while j <= length(B)  
 
     
    
 cnt=0;
    while true 
         
        cnt=cnt+1
        for ii=1:3
            if ii == 1
                aBCs(j)=aBCsupper(j);
            elseif ii==2
                aBCs(j)=aBCslower(j);
            elseif ii==3
                aBCs(j)=aBC(j);
            end
            
            
        BGen = B(j);   %Sets generic B values for current array index and overwrites for each iteration    
        y = y(1);   %Sets original initial conditions for y, z and a    
        z = z(1);   
        a = aBCs(j); %Sets a to guesses in f'' boundary conditions arrary 
    
        Fxyza = @(n,y,z,a) z;                  %dy/dn = z   A.K.A   f'    
        Gxyza = @(n,y,z,a) a;                  %dz/dn = a   A.K.A   f''   
        Hxyza = @(n,y,z,a) BGen*(z*z - 1) - y*a;  %da/dn = [function(n,y,z,a)]  A.K.A  f''' 
        

            
            %Solve f,f',f''        
            i = 1;            %Loop control variable;?    
            while i <= N;     %Runs for the interval from a to b??     



            %Runge-Kutta 4th order algorithm based off?        %k0, k1, k2 and k3 to estimate next y, z and a values??        
                k0(i) = h*Fxyza(n(i),y(i),z(i),a(i));     %Declares arrary for k0 values?        
                L0(i) = h*Gxyza(n(i),y(i),z(i),a(i));     %Declares arrary for L0 values?        
                M0(i) = h*Hxyza(n(i),y(i),z(i),a(i));     %Declares arrary for M0 values??        

                k1(i) = h*Fxyza(n(i) + .5*h, y(i) + .5*k0(i), z(i) + .5*L0(i), a(i) + .5*M0(i));  %Determines arrary for k1 values?        
                L1(i) = h*Gxyza(n(i) + .5*h, y(i) + .5*k0(i), z(i) + .5*L0(i), a(i) + .5*M0(i));  %Determines arrary for L1 values?        
                M1(i) = h*Hxyza(n(i) + .5*h, y(i) + .5*k0(i), z(i) + .5*L0(i), a(i) + .5*M0(i));  %Determines arrary for M1 values??        

                k2(i) = h*Fxyza(n(i) + .5*h, y(i) + .5*k1(i), z(i) + .5*L1(i), a(i) + .5*M1(i));  %Determines arrary for k2 values?        
                L2(i) = h*Gxyza(n(i) + .5*h, y(i) + .5*k1(i), z(i) + .5*L1(i), a(i) + .5*M1(i));  %Determines arrary for L2 values?        
                M2(i) = h*Hxyza(n(i) + .5*h, y(i) + .5*k1(i), z(i) + .5*L1(i), a(i) + .5*M1(i));  %Determines arrary for M2 values??        

                k3(i) = h*Fxyza(n(i) + h, y(i) + k2(i), z(i) + L2(i), a(i) + M2(i));  %Determines arrary for k3 values?        
                L3(i) = h*Gxyza(n(i) + h, y(i) + k2(i), z(i) + L2(i), a(i) + M2(i));  %Determines arrary for L3 values?        
                M3(i) = h*Hxyza(n(i) + h, y(i) + k2(i), z(i) + L2(i), a(i) + M2(i));  %Determines arrary for M3 values??        

                y(i+1) = y(i) + (1/6)*(k0(i) + 2*k1(i) + 2*k2(i) + k3(i));      %Using Rk4 Algorithm to estimate y at next step?        
                z(i+1) = z(i) + (1/6)*(L0(i) + 2*L1(i) + 2*L2(i) + L3(i));      %Using Rk4 Algorithm to estimate z at next step?        
                a(i+1) = a(i) + (1/6)*(M0(i) + 2*M1(i) + 2*M2(i) + M3(i));      %Using Rk4 Algorithm to estimate a at next step??        

                n(i+1) = n(i) + h;  %Increments to the next step size??        
                %End of Runge-Kutta algorithm??        

                i = i+1;            %Adds one to the loop control variable??   

                zend(ii)=z(end);%put f'(inf) in a vector

            end
        end
        
  %Put in half interval search here%%%%%%%%%%%    
        tol = 1e-11;
        if abs(zend(3)-1)<tol || cnt>10000 ||zend(1)==zend(2)
            zendfinal(j)=zend(3);
            break
%         elseif abs(zend(1)-1)<abs(zend(3)-1)
        elseif zend(3)>1
            aBCsupper(j)=aBC(j);
%         elseif abs(zend(2)-1)<abs(zend(3)-1)
        elseif zend(3)<1
            aBCslower(j)=aBC(j);    
        end
        aBC(j)=(aBCsupper(j)+aBCslower(j))/2;
   end

    %Graph of the boundary-layer velocity profiles?    
    %for the different values of B (Beta)?    
    figure(1)   
     hold on;    
     grid on;   
     plot(n,z,'LineWidth',1.5);    
     Legend{j}=sprintf('B = %g',B(j));

     axis([0 6 0 1]);  
     title('Non Dimensionalized Boundary Layer Velocity Profile');   
     ylabel('f\prime (u/U_{\infty})');   
     xlabel('\eta'); 
     
     eta99(j)=n(find(z>0.99,1));
     
     figure(2)
     hold on;plot(n/eta99(j),z,'LineWidth',1.5);axis([0 1 0 1])
     title('Non Dimensionalized Boundary Layer Velocity Profile');   
     ylabel('f\prime (u/U_{\infty})');   
     xlabel('\eta'); 
    
     j = j+1;  %Increments Beta(j) to next array value and runs till end of Beta array?
end
legend(Legend,'Location','Best')
%Plot f'' as func of Beta

figure
plot(B,aBCs,'LineWidth',1.5)
title('f"(0) vs \beta'),xlabel('\beta'),ylabel('f"(0)')

vpa(aBCs,15) %Final f''(0) that was found
vpa(zendfinal,15) %Final f'(inf) value


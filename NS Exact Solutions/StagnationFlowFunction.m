%This script creates the function  f(eta) for 2d Stagnation flow
%Solves for phi''' IC using half interval search

clear all

BC = 1.232534222946897;  %guess for f''(0) Boundary Condition ; Compare to f'(inf) = 1 for accuracy
BCupper=1.24;
BClower=1.23;
        
 cnt=0;
    while true 
         
        cnt=cnt+1;
        for ii=1:3
            if ii == 1
                BCs=BCupper;
            elseif ii==2
                BCs=BClower;
            elseif ii==3
                BCs=BC;
            end
            
           [eta,phi]=ode45(@Stagnation2dODE,[0,10],[0,0,BCs]); %% Difeq to solve
           
            phi_prime(ii)=phi(end,2);%put f'(inf) in a vector
        end
        
  %Put in half interval search here%%%%%%%%%%%    
        tol = 1e-15;
        if abs(phi_prime(3)-1)<tol || cnt>10000 ||phi_prime(1)==phi_prime(2)
            f_primeFinal=phi_prime(3);
            break
%         elseif abs(zend(1)-1)<abs(zend(3)-1)
        elseif phi_prime(3)>1
            BCupper=BC;
%         elseif abs(zend(2)-1)<abs(zend(3)-1)
        elseif phi_prime(3)<1
            BClower=BC;    
        end
        BC=(BCupper+BClower)/2;
    end
   
    cnt=cnt
    fprintf('phi''(inf)= %.15f\n',f_primeFinal)
    fprintf('phi"(0)= %.15f\n',BCs)
    
    plot(phi(:,2),eta)
    ylabel('\eta=sqrt(a/\nu)')
    xlabel('\phi''=u/U')
    title('Steady Stagnation Flow')
    hold on

    
    %% 3D Stagnation Flow
    
    
    
BC3 = 1.5;  %guess for f''(0) Boundary Condition ; Compare to f'(inf) = 1 for accuracy
BC3upper=2;
BC3lower=1;
        
 cnt=0;
    while true 
         
        cnt=cnt+1;
        for ii=1:3
            if ii == 1
                BC3s=BC3upper;
            elseif ii==2
                BC3s=BC3lower;
            elseif ii==3
                BC3s=BC3;
            end
            
           [eta,phi]=ode45(@Stagnation3dODE,[0,10],[0,0,BC3s]); %% Difeq to solve
           
            phi_prime(ii)=phi(end,2);%put f'(inf) in a vector
        end
        
  %Put in half interval search here%%%%%%%%%%%    
        tol = 1e-5;
        if abs(phi_prime(3)-1)<tol || cnt>10000 ||phi_prime(1)==phi_prime(2)
            phi3_primeFinal=phi_prime(3);
            break
%         elseif abs(zend(1)-1)<abs(zend(3)-1)
        elseif phi_prime(3)>1
            BC3upper=BC3;
%         elseif abs(zend(2)-1)<abs(zend(3)-1)
        elseif phi_prime(3)<1
            BC3lower=BC3;    
        end
        BC3=(BC3upper+BC3lower)/2;
    end
   
    cnt=cnt
    
    fprintf('2D phi''(inf)= %.15f\n',f_primeFinal)
    fprintf('2D phi"(0)= %.15f\n',BCs)
    fprintf('3D phi''(inf)= %.15f\n',phi3_primeFinal)
    fprintf('3D phi"(0)= %.15f\n',BC3s)
    
    plot(phi(:,2),eta)
    xlim([0,1.1])
    ylim([0,3.2])
    plot([.99,.99],[0,3.2])
    legend('2D Flow','Axisymmetric Flow')
    
    
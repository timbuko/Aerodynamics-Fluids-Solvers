%Solve Boundary Layer Using Pohlhausen Method to solve integral equation
%Using 4th order polynomial to approx BL profile

%by Timothy Bukowski for AME 70731 Viscous Flow

clear all; close all
tic
Uinf=23; %m/s
R=2.54/100;
dx=R/100;
%Air
nu=1.81e-5; %Pa*s
rho=1.275; %kg/m^3
%Water
% nu=1e-3; %Pa*s
% rho=1000; %kg/m^3
mu=rho*nu;

%Define Potential Flow Solution

ue=@(xx) 2*Uinf*sin(xx/R);
% ue=taylor(ue,xx,'Order',9);

syms xx
due=diff(ue(xx));
ddue=diff(due);

%% Initial Conditions

x(1)=R/1000;
K(1)=0.0770;
eta=[0:0.01:1.2];
Lambda(1)=7.052;
z(1)=K(1)/subs(due,xx,x);
% % z(1)=d2(1)/nu; %doesnt' work
dz(1)=-0.0652*subs(ddue,xx,x)/subs(due,xx,x)^2;
i=1;
while Lambda>=-12
    
    d2(i)=sqrt(nu*K(i)/eval(subs(due,xx,x(i))));
    d(i)=d2(i)/(1/63*(37/5-Lambda(i)/15-Lambda(i)^2/144));
    d1(i)=d(i)*(3/10-Lambda(i)/120);
    T(i)=mu*ue(x(i))/d(i)*(2+Lambda(i)/6);
%     eta(:,i)=y/d(i);
    u_ue(:,i)=(2*eta-2*eta.^3+eta.^4)...
        +Lambda(i)/6*(eta-3*eta.^2+3*eta.^3-eta.^4);
    F(i)=2*(37/315-1/945*Lambda(i)-1/9072*Lambda(i)^2)*...
        (2-116/315*Lambda(i)+(2/945+1/120)*Lambda(i)^2+2/9072*Lambda(i)^3);
    
    if i==1;F(1)=dz/dx*ue(x);end
    dz(i+1)=F(i)/ue(x(i))*dx;
    
    i=i+1;
%     if x(i-1)*180/pi>106.5;dx=0.001;end
    x(i)=x(i-1)+dx; 
    z(i)=z(i-1)+dz(i-1);
    K(i)=z(i)*subs(due,xx,x(i));
    
    %Bisection Method to solve for Lambda
        if i>4
            BCupper=Lambda(i-1);
        else
            BCupper=7.052;
        end
        BClower=-13;
        BC = (BCupper+BClower)/2;  %guess for Lambda ;
        cnt=0;
            while true 

                cnt=cnt+1;
                for ii=1:3
                    if ii == 1
                        L=BCupper;
                    elseif ii==2
                        L=BClower;
                    elseif ii==3
                        L=BC;
                    end

                  K_try(ii)=(37/315-1/945*L-1/9072*L^2)^2*L; %%%EQ to solve
                  
                end
  
                tol = 1e-6;
                if abs(K_try(3)-K(i))<tol || cnt>10000 ||K_try(1)==K_try(2)
                    Lambda(i)=L;
                    break
        %         elseif abs(zend(1)-1)<abs(zend(3)-1)
                elseif K_try(3)>K(i)
                    BCupper=BC;
        %         elseif abs(zend(2)-1)<abs(zend(3)-1)
                elseif K_try(3)<K(i)
                    BClower=BC;    
                end
                BC=(BCupper+BClower)/2;
            end
    if i>5/dx*0.01
        if Lambda(i)==Lambda(i-floor(3/dx*0.01))&Lambda(i-floor(3/dx*0.01))==Lambda(i-floor(5/dx*0.01))
            error('Lambda Calc is messed up')
        end
    end
    
%     figure(1)
%     plot(x(i-1)/R*180/pi,T(i-1)/(.5*rho*Uinf^2)*sqrt(Uinf*R/nu),'*');hold on
%     xlim([0 150])
% %     ylim([0 2.5])
if i>10000; break;end
    
end
toc
fprintf('\n\nSeparation at theta = %.1f\n',x(end-1)/R*180/pi)

%% Plots

%Plot Tau
close all
figure(1)
plot(x(1:length(T))./R*180/pi,T./(.5*rho*Uinf^2)*sqrt(Uinf*R/nu));hold on
xlim([0 120])
xlabel('Surface Location $\theta$ (Deg)',...
    'Interpreter','latex','fontsize',15,'fontweight','bold')
ylabel('$\frac{\tau_w}{\frac{1}{2} \rho U^2_\infty} \sqrt{\frac{U_\infty R}{\nu}}$',...
    'Interpreter','latex','fontsize',15,'fontweight','bold')

    % Add grid lines-----------------------------------------------------------
    arraynums = (0:20:120);
    for i = 1:length(arraynums)
        xline(arraynums(i),'color',[0,0,0]+0)
    end

    arraynums = (0:1:5);
    for i = 1:length(arraynums)
        yline(arraynums(i),'color',[0,0,0]+0)
    end

    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
         'YGrid', 'on',   'LineWidth', 1)  

%Plot Flow Profile
figure(2)
theta=[0:35:105,x(end-1)/R*180/pi];
xPlot=theta*R*pi/180;

y=eta*d(1);
plot(u_ue(:,1),y/R*sqrt(Uinf*R/nu));hold on;Legend{1}=sprintf('%stheta = %g','\',theta(1));
for i=2:length(xPlot)
   idx=find(x<=xPlot(i),1,'last');
   y=eta*d(i);
   plot(u_ue(:,idx),y/R*sqrt(Uinf*R/nu))
   Legend{i}=sprintf('%stheta = %g','\',theta(i));
end
xlim([0,1])
% ylim([0 3])
xlabel('$\frac{u}{U_\infty}$',...
    'Interpreter','latex','fontsize',15,'fontweight','bold')
ylabel('$\frac{y}{R} \sqrt{\frac{U_\infty R}{\nu}}$',...
    'Interpreter','latex','fontsize',15,'fontweight','bold')
grid on
legend(Legend,'Location','Northwest')

% %Plot F(K) vs K
% figure(3)
% plot(K(1:length(F)),F)
% xlim([-0.12, .08])
% ylim([0 1.3])
% grid on
% xlabel('K')
% ylabel('F(K)')

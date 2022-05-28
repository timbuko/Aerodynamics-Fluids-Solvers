%Taylor Maccoll Solve

%Assume
% Steady Shock, Axisymmetric Flow, Infinitely Long Sharp Cone, adiabatic,
% cal perf gas

clear all; close all
addpath('H:\My Drive\60630 Compressible Flow\Shock Functions')

n=3000;

% beta = linspace(0,pi/2,n); %Shock Angle
M = [1.25, 2, 6, 10];    %Mach numbers to solve
% M=[1.002:.001:1.005,1.01,1.025:.005:1.055,1.06:.05:1.5,1.6:.1:2,2.2:.2:4,4.5,5,6,8,10,20,9999]; %Full Range of M

conv=1e-5; %accuracy of theta wall

%Constants
gamma=1.4;


T=300; %K
P=1e5; %Pa
R=287; %J/kgK
rho=P/(R*T); %kg/m^3
a=sqrt(gamma*R*T);
Cp=gamma*R/(gamma-1);

%Calc Displacement angle right after shock
for i=1:length(M)
    for j=1:length(beta)
        
        %Calc deflection angle for 2D oblique shock
        delta(i,j) = atan(2*cot(beta(j))*(M(i)^2*sin(beta(j))^2-1)/(M(i)^2*(gamma+cos(2*beta(j)))+2));
        
        %Calc ur (ur doesn't change across a shock)
        ur(i,j)= M(i)*a*cos(beta(j));

        %Calc Mach after shock
        mn(i,j)=M(i)*sin(beta(j));
        mn2(i,j)=oblique(mn(i,j),gamma);
        M2(i,j)=mn2(i,j)/sin(beta(j)-delta(i,j));
        
        %nondimensional velocity u_tilde
        u_tild(i,j)=1/sqrt((2/((gamma-1)*M2(i,j)^2))+1);
        ur_tild(i,j)=u_tild(i,j)*cos(beta(j)-delta(i,j));
        ut_tild(i,j)=-u_tild(i,j)*sin(beta(j)-delta(i,j));
       
    end   
end


ut_tild(isinf(ut_tild))=0; %Set inf values to zero so they skip in ode45


%% Solve T-M

options=odeset('Events',@events,'RelTol',1e-11);
endtheta=1e-10; %cant actually go to zero

%Start solve at theta=beta the step theta towards zero
%Stop solve when u_theta=0 (means you are at wall)-->this theta is the cone angle
for i=1:length(M)
    for j=1:length(beta)
        [thetaWall,UrSol] = ode15s(@taylorMaccoll,[beta(j),endtheta],[ur_tild(i,j),ut_tild(i,j)],options,gamma);
        if abs(UrSol(end,2))<conv
            cone_angle(i,j)=thetaWall(end); %Cone half Angle
            
            for k=1:length(thetaWall)
                ur_tild(i,j,k)=UrSol(k,1);
                ut_tild(i,j,k)=UrSol(k,2);
                u_tild(i,j,k)=sqrt(ur_tild(i,j,k)^2+ut_tild(i,j,k)^2);
                M_post(k,j,i)=sqrt(2*u_tild(i,j,k)^2/((gamma-1)*(1-u_tild(i,j,k)^2))); %Mach behind shock
%                 if M_post(k,j,i)==0;M_post(k,j,i)=nan;end
                theta(k,j,i)=thetaWall(k); %need to flip theta order
                
            end
            M_post(M_post==0)=nan;theta(M_post==0)=nan;
            
            Mc(i,j)=sqrt(2*u_tild(i,j,k)^2/((gamma-1)*(1-u_tild(i,j,k)^2)));%Mach at cone surface
            Mc(Mc==0)=nan;
        else
            cone_angle(i,j)=nan;
            M_post(i,j,:)=nan;
            theta(i,j,:)=nan;
        end
    end
end

    
    
%% Plot Cone Theta Beta M

load TaylorMaccol_NACASolution

%Find Theta max (limit of weak shock)
% maxThetaLine=[0,90]; %XX
for i=1:length(M)
    [~,maxThetaIndex(i)]=max(cone_angle(i,:));
%      maxThetaLine(i+1,:)=[180/pi*cone_angle(i,maxThetaIndex(i)),180/pi*beta(maxThetaIndex(i))];%XX
    M_label(i,:)=[180/pi*cone_angle(i,maxThetaIndex(i)),180/pi*beta(maxThetaIndex(i))];
end
maxThetaLine=TaylorMaccol_NACASolution.maxThetaLine; %use max theta from large data set

%Smooth out max theta line
xinterp=linspace(min(maxThetaLine(:,1)),max(maxThetaLine(:,1)),100);
yinterp=interp1(maxThetaLine(:,1),maxThetaLine(:,2),xinterp);

%        Plot Conical Theta-Beta-M       %
figure()
c=hsv(length(M));%Get Colors for plot
for i=1:length(M)
    temp=~isnan(cone_angle(i,2:end));%get first point (fp) that is solved to extend to 0
    for j=1:length(temp)-1
        temp2(j)=double(temp(j))+double(temp(j+1));
    end
    fp(i)=find(temp2==2,1)+1;
    
    ca=[0,180/pi*cone_angle(i,fp(i):maxThetaIndex(i));180/pi*beta(fp(i)),180/pi*beta(fp(i):maxThetaIndex(i))]';
    plot(ca(:,1),ca(:,2),'LineWidth',2,'Color',c(i,:));hold on;
    text(M_label(i,1),M_label(i,2),num2str(M(i)),'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',9)
end
% plot(maxThetaLine(:,1),maxThetaLine(:,2)) %non smooth max theta
plot(xinterp,yinterp) %smoothed max theta
legend('M=1.25', 'M=2.0','M=6.0','M=10.0','\theta_{max}');
title('Conical flow shock angle vs cone angle');xlabel('Cone Angle, \theta_c (deg)');ylabel('Shock Angle, \theta_s (deg)')

%%           Plot Concical and 2D Theta-Beta-M      %
figure()
hh=plot(xinterp,yinterp,'k','LineWidth',2);hold on %smoothed max theta
hhh=plot(TaylorMaccol_NACASolution.maxThetaLine2D(:,1),TaylorMaccol_NACASolution.maxThetaLine2D(:,2),'--k');
ylim([0,90]);xlim([0,60]);title('2D \theta-\beta-M');xlabel('\theta');ylabel('\beta')
for i=1:length(M)
    plot([0,180/pi*cone_angle(i,fp(i):end)],[180/pi*beta(fp(i)),180/pi*beta(fp(i):end)],'LineWidth',2,'Color',c(i,:))
    hold on;
    plot(delta(i,:)*180/pi,beta*180/pi,'--','Color',c(i,:))
    text(M_label(i,1),M_label(i,2),num2str(M(i)),'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',9)

end
legend([hh,hhh],'Conical','Wedge')


%    Plot Mach at wall vs cone angle     %
figure()
for i=1:length(M)
    plot(cone_angle(i,2:maxThetaIndex(i))*180/pi,Mc(i,1:maxThetaIndex(i)-1))
    hold on
end
title('Mach at cone surface M_c vs cone half angle'),xlabel('\theta_c (deg)'),ylabel('M_c')
legend('M=1.25', 'M=2.0','M=6.0','M=10.0')


%%   Calc and plot p/p_inf and M vs theta for cone half angle of 7, M_inf=6
angle=7;
mach=6;

[~,ind]=min(abs(cone_angle(find(M==mach),:)*180/pi-angle));%find cone angle closest to 7 for m=6
[rr,cc]=find(cone_angle==cone_angle(find(M==mach),ind));
i=rr(1);
j=cc(1);



% Plot Mach vs Theta
figure()
plot(theta(:,j,i)*180/pi,M_post(:,j,i))
title(['Mach vs Theta for Cone Angle of ',num2str(cone_angle(i,j)*180/pi),' and M_{\infty}=6']),xlabel('Theta (deg)'),ylabel('Mach')
set (gca,'Xdir','reverse')

% %Calc and plot P/P_inf vs Theta
p2_pinf=presratoblique(mn(i,j),gamma);
[~,~,p2_po,~,~] = flowisentropic(gamma, M2(i,j), 'mach');
for k=1:length(M_post(:,1,1))
[~,~,p_po(k),~,~] = flowisentropic(gamma, M_post(k,j,i), 'mach');
p_pinf(k)=p2_pinf/p2_po*p_po(k);
end

figure()
plot(theta(:,j,i)*180/pi,p_pinf)
title(['P/P_{\infty} vs Theta for Cone Angle of ',num2str(cone_angle(i,j)*180/pi),' and M_{\infty}=6']),xlabel('Theta (deg)'),ylabel('P/P_{\infty}')
set (gca,'Xdir','reverse')


%% 3D plot of Theta Beta M

% bb=beta;%matlab doesnt like beta name
% figure()
% subplot(2,2,1)%3D view
% for i=1:length(M)
%     plot3(M(i)*ones(length(bb)),bb*180/pi,cone_angle(i,:)*180/pi,'Color',c(i,:)) %Conical
%     hold on
% end
% for i=1:length(M)
%     plot3(M(i)*ones(length(bb)),bb*180/pi,delta(i,:)*180/pi,'--','Color',c(i,:)) %Wedge
% end
% plot3([0,TaylorMaccol_NACASolution.M(1:end-1)],maxThetaLine(1:end-1,2),maxThetaLine(1:end-1,1),'k')
% plot3([0,TaylorMaccol_NACASolution.M(1:end-1)],TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,2),TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,1),'--k');
% 
% xlabel('Mach'),ylabel('Beta'),zlabel('Theta'),ylim([0,90]),zlim([0,60])
% grid
% 
% subplot(2,2,2)%theta vs beta
% for i=1:length(M)
%     plot(bb*180/pi,delta(i,:)*180/pi,'--','Color',c(i,:))
%     hold on
% end
% for i=1:length(M)
%     plot(bb*180/pi,cone_angle(i,:)*180/pi,'Color',c(i,:))
% end
% plot(maxThetaLine(1:end-1,2),maxThetaLine(1:end-1,1),'k')
% plot(TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,2),TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,1),'--k');
% xlabel('Beta'),ylabel('Theta'),xlim([0,90]),ylim([0,60])
% 
% subplot(2,2,3)%mach vs beta
% for i=1:length(M)
%     plot(M(i)*ones(length(bb)),bb*180/pi,'--','Color',c(i,:))
%     hold on
% end
% for i=1:length(M)
%     plot(M(i)*ones(length(bb)),bb*180/pi,'Color',c(i,:))
% end
% plot([0,TaylorMaccol_NACASolution.M(1:end-1)],maxThetaLine(1:end-1,2),'k')
% plot([0,TaylorMaccol_NACASolution.M(1:end-1)],TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,2),'--k');
% xlabel('Mach'),ylabel('Beta'),ylim([0,90])
% 
% subplot(2,2,4)%theta vs mach
% for i=1:length(M)
%     plot(M(i)*ones(length(bb)),delta(i,:)*180/pi,'--','Color',c(i,:))
%     hold on
% end
% for i=1:length(M)
%     plot(M(i)*ones(length(bb)),cone_angle(i,:)*180/pi,'Color',c(i,:))
% end
% plot([0,TaylorMaccol_NACASolution.M(1:end-1)],maxThetaLine(1:end-1,1),'k')
% plot([0,TaylorMaccol_NACASolution.M(1:end-1)],TaylorMaccol_NACASolution.maxThetaLine2D(1:end-1,1),'--k');
% xlabel('Mach'),ylabel('Theta'),ylim([0,60])

%% Misc other codes I did

% %%%%%%%%%% Plot a bunch of Mach
% figure()
% b=TaylorMaccol_NACASolution.beta;m=TaylorMaccol_NACASolution.M;mti=TaylorMaccol_NACASolution.maxThetaIndex;mtl=TaylorMaccol_NACASolution.maxThetaLine;caa=TaylorMaccol_NACASolution.cone_angle;
% c=hsv(length(m));
% for i=1:length(m)
%     plot(caa(i,1:mti(i)),b(1:mti(i)),'LineWidth',2,'Color',c(i,:));hold on;
%     temp=~isnan(caa(i,2:end));%get first point (fp) that is solved to extend to 0
%         for j=1:length(temp)-1
%             temp2(j)=double(temp(j))+double(temp(j+1));
%         end
%         fp(i)=find(temp2==2,1)+1;
%     ca=[0,caa(i,fp(i):mti(i));b(fp(i)),b(fp(i):mti(i))]';
%     plot(ca(:,1),ca(:,2),'LineWidth',2,'Color',c(i,:));hold on;
% end
% % plot(maxThetaLine(:,1),maxThetaLine(:,2)) %non smooth max theta
% plot(xinterp,yinterp) %smoothed max theta
% figure()
% for i=1:length(m)-1
%     plot3(m(i)*ones(length(b),1),b,caa(i,:)) %Plot 3D plot
%     hold on
% end


% % %Save file with a bunch of Mach
% % TaylorMaccol_NACASolution.beta=beta*180/pi;
% % TaylorMaccol_NACASolution.M=M;
% % TaylorMaccol_NACASolution.maxThetaIndex=maxThetaIndex;
% % TaylorMaccol_NACASolution.maxThetaLine=maxThetaLine;
% % TaylorMaccol_NACASolution.cone_angle=cone_angle*180/pi;
% % 
% % save TaylorMaccol_NACASolution TaylorMaccol_NACASolution


% % % Calc and save 2D max theta line
% M=TaylorMaccol_NACASolution.M;
% for i=1:length(M)
%     for j=1:length(beta)
%         
%         %Calc deflection angle for 2D oblique shock
%         delta(i,j) = atan(2*cot(beta(j))*(M(i)^2*sin(beta(j))^2-1)/(M(i)^2*(gamma+cos(2*beta(j)))+2));
%     end
% end
% 
%  maxThetaLine2D=[0,90];
% for i=1:length(M)
%     [~,maxThetaIndex2D(i)]=max(delta(i,:));
%      maxThetaLine2D(i+1,:)=[180/pi*delta(i,maxThetaIndex2D(i)),180/pi*beta(maxThetaIndex2D(i))];
% end
% xinterp=linspace(min(maxThetaLine2D(:,1)),max(maxThetaLine2D(:,1)),100);
% yinterp=interp1(maxThetaLine2D(:,1),maxThetaLine2D(:,2),xinterp);
% 
% TaylorMaccol_NACASolution.maxThetaLine2D=maxThetaLine2D;%[xinterp',yinterp'];



% % % % delete bad points
% % % ind=[find(M==1.11),find(M==1.46),find(M==1.002),find(M==1.035),find(M==1.003)];
% % % cone_angle(ind,:)=[];
% % % maxThetaLine(ind+1,:)=[];
% % % maxThetaIndex(ind)=[];
% % % M(ind)=[];
%%
% Other Plotter
% title({'Taylor and Maccoll Solution';'Right Circular Cone at \alpha = 0'},'fontsize',18)
% 
% xlabel('Cone Angle $\theta_c$ [degrees]',...
%     'Interpreter','latex','fontsize',20,'fontweight','bold')
% ylabel('Shock Angle $\theta_s$ [degrees]',...
%     'Interpreter','latex','fontsize',20,'fontweight','bold')
% xlim([0 60])
% ylim([0 90])
% hold on
% 
% Add grid lines-----------------------------------------------------------
% arraynums = (0:.5:60);
% for i = 1:length(arraynums)
%     xline(arraynums(i),'color',[0,0,0]+0)
% end
% 
% arraynums = (0:1:90);
% for i = 1:length(arraynums)
%     yline(arraynums(i),'color',[0,0,0]+0)
% end
% 
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'LineWidth', 1)  
% caxis([0, .5]);
% legend('M_{\infty}')

%%

function [value,isterminal,direction] = events(theta,ur,gamma)
% ODE Solver option to Locate the time when u_theta=0 and stop integration.

value=zeros(2,1);
isterminal=zeros(2,1);
direction=zeros(2,1);

%Stop when u_theta goes positive at the wall (u_theta=0 at wall)
value(2)=1;
if (ur(2)>0.0)
    value(2)=0.0;
end
isterminal(2)=1; % terminal on ur(2)=0 (stops when u_theta=dur_dtheta=0)
direction(2)=0; % detects from rising or falling slope

end
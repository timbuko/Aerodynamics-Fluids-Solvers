%Solve Shock Tube through time (Only solves till shock hits wall)

%rho,p,u,M,up, and speeds of waves are all checked against data
%But h0 is still increasing with time

%Check assumptions that go into equations to see if theres something that
%would allow h0 to increase. Where would the energy be coming from

%%%
%values from https://altairuniversity.com/wp-content/uploads/2012/08/Example_13.pdf
% p1= 0.2e5 ;%Pa
% p4= 5e5 ; %Pa
% T1= 303 ; %K
% T4= 303 ; %K
%%%

%%%
%values from https://www.grc.nasa.gov/www/wind/valid/stube/stube.html
% p1= 6894.76 ;%Pa (1psia)
% p4= 68947.6 ; %Pa (10 psia)
% T1= 231.111 ; %K (416R)
% T4= 288.889 ; %K (520R)
% %%

% %values from Anderson Compressible Flow problem 7.10
% p1= 1e5 ;%Pa (1atm)
% p4= 5*p1 ; %Pa (10 atm)
% T1= 300; %K (416R)
% T4= 300; %K (520R)

clear all, close all

%ASSUMPTIONS
%calorically perfect gas - constant gamma,R
%1-D, adiabatic, no body force,
%assumption that u3=u2 (made in solving expansion section)
%constant flow properties across u2 and across u3
%no viscous or thermal diffusion (inviscid)
%constant properties in region 1 and 2 ->steady state req for Hugioniot eq

%% Inputs
%gas 1
gamma1=1.4;
R1=287; %J/kgK
Cp1=gamma1*R1/(gamma1-1);

%gas 4
gamma4=gamma1;
R4=R1;
Cp4=gamma4*R4/(gamma4-1);

%%Initial Condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p1= 1e5 ;%Pa (1atm)
p4= 5*p1 ; %Pa (10 atm)
T1= 300; %K (416R)
T4= 300; %K (520R)

%SET UP TUBE
L = 12; %m Length of Shock Tube
diaphragm = 3/12; %Diaphragm halfway down tube 
xorigin=L*diaphragm;

n_t= 1e-4; %Time Step Size
n_x=1000;  %discretization of tube (# of elements)

%% Make x-t Diagram
%Find Speed of the Shock, contact surface, and expansion waves
syms ppp mr
a1=sqrt(gamma1*R1*T1);
a4=sqrt(gamma4*R4*T4);
temp=ppp/p1*(1-((gamma4-1)*a1/a4)*(ppp/p1-1)/(sqrt(2*gamma1*(2*gamma1+(gamma1+1)*(ppp/p1-1)))))^(-2*gamma4/(gamma4-1));
p2=vpa(solve(p4/p1==temp,ppp));
p2=p2(real(p2)>0&imag(p2)==0);




W=a1*sqrt((gamma1+1)/(2*gamma1)*(p2/p1-1)+1); %Speed of the incident shock
up=a1/gamma1*(p2/p1-1)*sqrt(2*gamma1/(gamma1+1)/(p2/p1+(gamma1-1)/(gamma1+1))); %Speed of Contact Surface

Ms=W/a1;
Mr=vpa(solve(mr/(mr^2-1)==Ms/(Ms^2-1)*sqrt(1+2*(gamma1-1)/(gamma1+1)^2*(Ms^2-1)*(gamma1+1/Ms^2)),mr)); 
Mr=Mr(real(Mr)>0&imag(Mr)==0);
a2=(W-up)/sqrt(1+(gamma1-1)/2*(W/a1)^2/(gamma1*(W/a1)^2-(gamma1-1)/2));
Wr=Mr*a2+up; %Speed of Reflected Shock
% p5=p2*(1+2*gamma1/(gamma1-1)*(Mr^2-1));
% a5=Wr/sqrt(1+(gamma1-1)/2*Mr^2/(gamma1*Mr^2-(gamma1-1)/2));


u_expHead=a4;
u_expTail=a4-up*(gamma4-1)/2;


%Solve until shock or expansion hits wall (whichever is first)
t_sReflect = (L-xorigin)/W; %Shock hit wall
t_eReflect = xorigin/a4; %Expansion hits wall
t_sMeetCS = t_sReflect+solve(up*(t_sReflect+mr)+xorigin==L-Wr*(mr),mr); %Reflected Shock hits contact surface
tmax = t_sMeetCS; %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Set which time I'm plotting to  %%%%%%%%%%%%%%
tmaxsolve= t_eReflect; %%%%%%%% Set which time I'm solving to
    
%Plot x-t Diagram
figure(1);hold on
plot([xorigin,xorigin-u_expHead*tmax],[0,tmax],'g');xlabel('Location'),ylabel('Time'),xlim([0,L])
plot([xorigin,xorigin+(up-u_expTail)*tmax],[0,tmax],'g');
plot([xorigin,up*tmax+xorigin],[0,tmax],'b');
plot([xorigin,W*t_sReflect+xorigin],[0,t_sReflect],'m');
if tmax==t_sMeetCS
    plot([L,L-Wr*(tmax-t_sReflect)],[t_sReflect,tmax],'m');
end



%% Solve Properties across whole tube for each time step
% Region 1 %
    T1=T1; %K
    rho1=p1/(R1*T1); %kg/m^3
    u1=0; 
    a1=sqrt(gamma1*R1.*T1); %m/s
    M1=u1./a1;
    T01=findT0_T(M1,gamma1).*T1; %K
    h01=Cp1.*T01;

    % Region 4 %
    p4=p4;
    T4=T4; %K
    rho4=p4/(R4*T4); %kg/m^3
    u4=0;
    a4=sqrt(gamma4*R4.*T4); %m/s
    M4=u4./a4; 
    T04=findT0_T(M4,gamma4).*T4; %K
    h04=Cp4*T04;
    
    
for t=[1:floor(tmaxsolve/n_t)]    
    fprintf('Step:%d/%d\n',t,floor(tmaxsolve/n_t))
    
    x=[1:n_x];
    x4=x(1:floor((xorigin-u_expHead*t*n_t)/L*n_x)); %region 4 (driver) (0:15)
    xx=x(floor((xorigin-u_expHead*t*n_t)*n_x/L)+1:floor((xorigin+(up-u_expTail)*t*n_t)*n_x/L));%expansion region (15:35)
    if floor((xorigin-u_expHead*t*n_t)*n_x/L)==floor((xorigin+(up-u_expTail)*t*n_t)*n_x/L);xx=x(floor((xorigin-u_expTail*t*n_t)*n_x/L));end
    x3=x(floor((xorigin+(up-u_expTail)*t*n_t)*n_x/L)+1:floor((xorigin+up*t*n_t)*n_x/L));%region 3 (35:60)
    x2=x(floor((xorigin+up*t*n_t)*n_x/L)+1:floor((xorigin+W*t*n_t)*n_x/L));%region 2 (60:85)
    if floor((xorigin+up*t*n_t)*n_x/L)==floor((xorigin+W*t*n_t)*n_x/L);x2=x(floor((xorigin+W*t*n_t)*n_x/L));end
    x1=x(floor((xorigin+W*t*n_t)*n_x/L)+1:end);% region 1 (driven) (85:100
    
    

    % Region 1 fill out vector for plotting
    p(x1,t)=p1;
    T(x1,t)=T1; %K
    rho(x1,t)=rho1; %kg/m^3
    u(x1,t)=0; 
    a(x1,t)=sqrt(gamma1*R1.*T(x1,t)); %m/s
    M(x1,t)=u(x1,t)./a(x1,t);
    T0(x1,t)=findT0_T(M(x1,t),gamma1).*T(x1,t); %K
    h0(x1,t)=Cp1.*T0(x1,t);
    % Region 4 fill out vector for plotting
    p(x4,t)=p4;
    T(x4,t)=T4; %K
    rho(x4,t)=rho4; %kg/m^3
    u(x4,t)=0;
    a(x4,t)=sqrt(gamma4*R4.*T(x4,t)); %m/s
    M(x4,t)=u(x4,t)./a(x4,t); 
    T0(x4,t)=findT0_T(M(x4,t),gamma4).*T(x4,t); %K
    h0(x4,t)=Cp4*T0(x4,t);

    % Region 2 after shock %
    %Anderson eq for p2, u2, T2, rho2
    p1=p(x1(1),t);
    p(x2(end),t)=p2;
    u(x2(end),t)=a(x1(1),t)./gamma1*(p(x2(end),t)./p1-1).*sqrt((2*gamma1/(gamma1+1))/(p(x2(end),t)./p(x1(1),t)+(gamma1-1)/(gamma1+1)));
    T(x2(end),t)=T1*(p(x2(end),t)./p1*((gamma1+1)/(gamma1-1)+p(x2(end),t)/p1)/(1+p(x2(end),t)/p1*(gamma1+1)/(gamma1-1)));
    rho(x2(end),t)=rho1*(1+(gamma1+1)/(gamma1-1)*p(x2(end),t)/p1)/((gamma1+1)/(gamma1-1)+p(x2(end),t)/p1);
%     rho(x2(end),t)=p(x2(end),t)./(R1*T(x2(end),t));
    a(x2(end),t)=sqrt(gamma1*R1.*T(x2(end),t)); %m/s
    M(x2(end),t)=u(x2(end),t)/a(x2(end),t);
    T0(x2(end),t)=findT0_T(M(x2(end),t),gamma1).*T(x2(end),t); %K
    h0(x2(end),t)=Cp1*T0(x2(end),t);

    % Rest of Region 2 %
    p(x2,t)=p(x2(end),t);
    u(x2,t)=u(x2(end),t);
    T(x2,t)=T(x2(end),t);
    rho(x2,t)=rho(x2(end),t);
    a(x2,t)=a(x2(end),t);
    M(x2,t)=M(x2(end),t);
    T0(x2,t)=T0(x2(end),t);
    h0(x2,t)=h0(x2(end),t);
    
    % Region of expansion %
    u(x3,t)=u(x2(1),t);
    u(xx,t)=linspace(u(x4(end),t),u(x3(1),t),length(xx));
    %u(xx)=2/(gamma4-1)*(a4+xx./t);

    %equations from notes
    a(xx,t)=a4*(1-(gamma4-1)/2*(u(xx,t)/a4));
    T(xx,t)=T4*(1-(gamma4-1)/2*(u(xx,t)/a4)).^2;
    p(xx,t)=p4*(1-(gamma4-1)/2*(u(xx,t)/a4)).^(2*gamma4/(gamma4-1));
    rho(xx,t)=rho4*(1-(gamma4-1)/2*(u(xx,t)/a4)).^(2/(gamma4-1));

    M(xx,t)=u(xx,t)./a(xx,t);
    T0(xx,t)=findT0_T(M(xx,t),gamma4).*T(xx,t); %K
    h0(xx,t)=Cp4*T0(xx,t);

    % Region 3 (assume constant properties in 3)
    T(x3,t)=T(xx(end),t);
    p(x3,t)=p(xx(end),t);
    a(x3,t)=sqrt(gamma4*R4.*T(x3,t)); %m/s
    M(x3,t)=u(x3,t)./a(x3,t); 
    rho(x3,t)=rho(xx(end),t);
    T0(x3,t)=findT0_T(M(x3,t),gamma4).*T(x3,t); %K
    h0(x3,t)=Cp4*T0(x3,t);

    T0Avgxx(t)=sum(T0(xx,t))/range(xx);%Average T0 over the expansion region
    T0Avg23(t)=sum(T0([x2,x3],t))/range([x2,x3]); %Average T0 over region 2&3
    rhoAvgxx(t)=sum(rho(xx,t))/range(xx); %Average density of expansion
    rhoAvg23(t)=sum(rho([x2,x3],t))/range([x2,x3]);
    rhoAvgnotx(t)=sum(rho([x1,x4,x2,x3],t))/range([x1,x4,x2,x3]);

end
p0= p.*findP0_P(M,gamma1);


%% Plot

param_name='Speed of Sound';
param = a;  %%%%% Choose what you want to plot (M,p,T,T0,rho,u,a,h0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make Gif
[~,tmstp]=size(param);
figure(2)
for i=1:tmstp
    pp=plot(x/n_x,param(:,i));
    title('Shock Tube until Shock hits wall'),ylabel(param_name),xlabel('Location')
%     ylim([0,max(max(param))])
    Mov(i)=getframe;

end

%Play Gif
gif_length = 3; %seconds
rep = 2;        %how many times movie reps
movie(Mov,rep,tmstp/gif_length)

%Plot Integral of T0 vs time
figure(4)
plot([n_t:n_t:tmaxsolve],sum(T0))
xlabel('Time')
ylabel('Integral of Total Temperature over tube')


%Plot a  couple things
t=40;

figure
plot(x,M(:,t),'color',[0 0.4470 0.7410],'LineWidth',2)
yyaxis right
plot(x,p0(:,t),'k','LineWidth',2)
hold on
plot(x,p(:,t),'k','LineWidth',2)
legend('M','P','P_0')

figure
plot(x,M(:,t),'color',[0 0.4470 0.7410],'LineWidth',2)
yyaxis right
plot(x,T0(:,t),'color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(x,T(:,t),'color',[0.8500 0.3250 0.0980],'LineWidth',2)
legend('Mach','T_0','T')

%%%%%%%% Plot Characteristic Lines %%%%%%%%%%%
%Plot x-t Diagram
figure(10);hold on
plot([xorigin,xorigin-u_expHead*tmax],[0,tmax/n_t],'g');xlabel('Location'),ylabel('Time'),xlim([0,L])
plot([xorigin,xorigin+(up-u_expTail)*tmax],[0,tmax/n_t],'g');
plot([xorigin,up*tmax+xorigin],[0,tmax/n_t],'b');
plot([xorigin,W*t_sReflect+xorigin],[0,t_sReflect/n_t],'m');
if tmax==t_sMeetCS
    plot([L,L-Wr*(tmax-t_sReflect)],[t_sReflect/n_t,tmax/n_t],'m');
end

%Plot Characteristic lines
num1=10; %decrease for more lines in x
num2=4;  % for lines in y
xx=x([1:num1:end])/n_x*L;
tt=[1:num2:tmaxsolve/n_t];
quiver(xx,tt,...
    (u([1:num1:end],[1:num2:end])+a([1:num1:end],[1:num2:end]))',...
    ones(size(a([1:num1:end],[1:num2:end])'))/n_t,'ShowArrowHead','off')
title('Characteristic line slopes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Comments
% % T0 is increasing with time. This doesn't make sense because there is no heat
% % input. I also expected avg T0 under the expansion to be constant, but thats
% % not the case either. So notsure if this is an algorithm thing or if its physical. 
% % 
% % The region between contact surface and shock is growing faster than the 
% % region from expansion head to contact surface, so I guess I would expect T0 
% % to increase. 
% % 
% % Maybe T0 will decrease again with shock and exp. reflection.
% 
% 
% % The easier to follow idea is that h0 should be constant 
% % h0=e+p/rho+1/2u^2, so we can see the breakup of all the energy forms
% % however I found three different ways of calulating 'e' and they are all
% % equal
% % (the three h0 are much closer for NASA example than for Altair example)
% internalEnergy1=h0-p./rho-.5*u.^2; %from def of h0
% internalEnergy2=5/2*p./rho; %internal energy of ideal gas
% internalEnergy3=(Cp1-R1)*T; %e=cvT
% 
% figure;hold on;plot(internalEnergy1(:,t),'LineWidth',4);plot(internalEnergy2(:,t),'LineWidth',2);plot(internalEnergy3(:,t),'--','LineWidth',2)
% legend('e=h_0-P\nu-1/2u^2','e=5/2p\nu','e=c_vT'),title('Internal Energy')
% 
% h02=internalEnergy2+p./rho+0.5*u.^2;h03=internalEnergy3+p./rho+0.5*u.^2;
% figure;hold on;plot(h0(:,t),'LineWidth',4);plot(h02(:,t),'LineWidth',2);plot(h03(:,t),'--','LineWidth',2),title('h0')
% legend('e=h_0-P\nu-1/2u^2','e=5/2p\nu','e=c_vT')
% 
% figure;subplot(1,3,1);plot(sum(h0));subplot(1,3,2);plot(sum(h02));title('Integral of h_0 vs Time');subplot(1,3,3);plot(sum(h03))
% 
% figure
% plot(x,internalEnergy1(:,t),'LineWidth',2)
% hold on
% plot(x,p(:,t)./rho(:,t),'LineWidth',2)
% plot(x,.5*u(:,t).^2,'LineWidth',2)
% plot(x,h0(:,t),'LineWidth',2)
% legend('Internal Energy','P\nu','Specific KE','h_0')

%three different methods of calc internal energy all the same, two of which
%doesn't include finding T0, so my calculations are coherent. and The
%isentropic relation used to find T0 is not faulty. Also I know that T is
%calculated right cuz P and rho are calc without T and i get same internal
%energy 5/2p/rho=CvT

% So I'm pretty convinced h0 is not constant, but where would this energy be coming from??
% There is potential energy of the material holding back the higher pressure gas that is released when
% the diaphragm breaks, so that makes sense why you might have some extra energy as first
% but how do you get a linear increase with time?


%% Functions

function T0_T = findT0_T(M,gamma)
    T0_T=1+(gamma-1)/2*M.^2;
end

function P0_p = findP0_P(M,gamma)
    P0_p=(1+(gamma-1)/2*M.^2).^(gamma/(gamma-1));
end

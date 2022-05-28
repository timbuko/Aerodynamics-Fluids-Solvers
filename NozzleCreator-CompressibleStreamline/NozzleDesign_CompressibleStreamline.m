%Calculate nozzle shape using compressible streamline theory
%
% This code used compressible streamline theory derived in Dr. Jumper's
% viscous flow course. 
clear; close all;
tic
%So we are neglecting heat transfer for the flat plate 
%If you make dy too small it wont work (like 1/20 is too small)

%% Inputs
%Start nozzle with a boundary layer growing along a flat plate

HeightThroat=0.04/2 /100; %m
Twall=450;  %K Wall Temperature
Uinf=340; %m/s inlet velocity
Pinf=101325;  %Pa inlet static pressure
Tinf=2000;    %K inlet static temperature

dx=0.000004; 
dy=HeightThroat/22;
y_centerInlet=0.045/100;%m centerline height at the inlet 
    %will be defining nozzle by the throat size so this y is arbitrary. May
    %need to adjust inlet profile if inlet is smaller than the BL, since
    %flow will be fully developed
    
%Constants
R=287; %J/kgK
Pr=0.74; %Prandtl #  μCp/k

%Flags
PressureProfile =1; % 1 for Notes; 2 for AIAA

% k=110 /1000; %W/mK
% gamma=1.4; %Heat capacity ratio
% cp=1004;  %J/kgK

%% Input Nozzle Pressure Profile
%Calc values
rhoInf=Pinf/(R*Tinf); %kg/m^3 inlet density 
cp=@(T) 1000*.9584*exp(0.000168*T); %J/kgK
gamma=@(T) 1000*.9584*exp(0.000168*T)./( 1000*.9584*exp(0.000168*T)-R);%cp/(cp-R);
mu=@(T) 18.27e-6*((291.15+120)./(T+120).*(T/291.15).^1.5);
k=@(T) 6.4779e-05*T+0.0062; %W/mK
nu=mu(Tinf)/rhoInf;

%Grid
x=-0.4/100:dx:0.4/100;
if PressureProfile==2;x=-0.4/100:dx:0.3/100;end
y=[0:dy:y_centerInlet]'; %DEFINE GRID IN Y

%Get points for pressure profile
if PressureProfile == 1 %Using profile from notes
    P_theory=@(x) 1-0.979*exp(-(9.7e-3)./(x*100/2.54-0.018));
    P_x=[0 0.005 0.00725 0.01 0.0125 0.015 0.0175 0.02 0.0225 0.023]/100*2.54;
    P_y=[1 0.9951154804 0.9896187008 0.9798615245 0.9677467384...
        0.9519943399 0.9317019646 0.905188353 0.8686136719 0.8593138331];

    cf=polyfit(P_x*1000,P_y,4);
    top=find(x<=P_x(end),1,'Last');bot=find(x>P_x(end),1);
    Pfit=polyval(cf,x(1:top)*1000);
    Pfit(1:find(Pfit>1,1,'last'))=1;
    P_P0=([Pfit(1:top),P_theory(x(bot:end))]);
elseif PressureProfile == 2  %Using profile from AIAA paper
    P_P0=0.55-0.45*tanh(8800*(x-0.00055))-30*x;
    P_P0(P_P0>1)=1;
end

for i=1:50
   P_P0=smooth(P_P0); 
end
throat=x(find(P_P0<0.528,1));
figure(11);subplot(2,1,1);
plot(x*100,P_P0,'k','Linewidth',2);hold on
ylabel('P/P_0');xlabel('x [cm]')
title('Pressure Ratio Distribution in Nozzle')
grid on

%% Preallocate stuff
U=zeros(length(y),length(x));
T=ones(length(y),length(x))*Twall;
M=zeros(length(y),length(x));
dudPsi=zeros(length(y),length(x)-1);
d2ud2Psi=zeros(length(y),length(x)-1);
dMuRhoUdPsi=zeros(length(y),length(x)-1);
dkRhoUdPsi=zeros(length(y),length(x)-1);
dTdPsi=zeros(length(y),length(x)-1);
d2Td2Psi=zeros(length(y),length(x)-1);
rho=zeros(length(y),length(x));
y_streamlines=repmat(y,1,length(x)-1);

%% Inlet flat plate calculation
%Calc boundary layer development using Blasius 

xin=1:find(x>=0,1);
for i=xin(2:end)
%Inlet velocity profile
tspan=y*sqrt(Uinf/(nu*(x(i)+0.004)));
if max(tspan)<5;error('Span of y values is too short and doesn''t include full boundary layer');end
[eta,f]=ode45(@blasiusode,tspan,[0,0,0.33206]);
U(:,i)=f(:,2)*Uinf;
    
%Inlet temperature profile
tspan=y*Pr^(-1/3)*sqrt(Uinf/(nu*(x(i)+0.004)));
[eta,f]=ode45(@blasiusode,tspan,[0,0,0.33206]);
T(:,i)=f(:,2)*(Tinf-Twall)+Twall;
end
U(2:end,1)=ones(length(y)-1,1)*Uinf;T(2:end,1)=ones(length(y)-1,1)*Tinf;

%Calc Mach and rho
M(:,xin)=U(:,xin)./sqrt(gamma(T(:,xin))*R.*T(:,xin));
rho(:,xin)=Pinf./(R*T(:,xin));

%Convert y to Psi
ruAvg=(rho(1:end-1,xin(end)).*U(1:end-1,xin(end))+rho(2:end,xin(end)).*U(2:end,xin(end)))/2;
dPsi=ruAvg.*dy;
Psi=cumsum([0;dPsi]);

% figure
% plot(y*1000,U(:,xin(end)));hold on;ylabel('Velocity [m/s]')
% yyaxis right;
% plot(y*1000,T(:,xin(end)));ylabel('Temperature [K]')
% xlabel('Distance from wall [mm]')
% title('Profiles at nozzle inlet, x=0')

%% Nozzle calculation
P=P_P0*Pinf';  %Assuming Pinf=P0 since our P_P0 dist starts at 1.  
              % based off Minf our P_P0=0.9055
dPdx=diff(P)/dx;

for j=length(xin):length(x)-1
       %Calc next U profile
       dudPsi(:,j)=Order1(U(:,j),dPsi);
       d2ud2Psi(:,j)=Order2(U(:,j),dPsi);
       dMuRhoUdPsi(:,j)=Order1(mu(T(:,j)).*rho(:,j).*U(:,j),dPsi);
       U(2:end,j+1)=U(2:end,j)+dx*(...
           -1./(rho(2:end,j).*U(2:end,j))*dPdx(j)...
           +dMuRhoUdPsi(2:end,j).*dudPsi(2:end,j)...
           +mu(T(2:end,j)).*rho(2:end,j).*U(2:end,j).*d2ud2Psi(2:end,j)); U(1,j+1)=U(1,j);    
       
       %Calc next T profile
       dkRhoUdPsi(:,j)=Order1(k(T(:,j)).*rho(:,j).*U(:,j),dPsi);
       dTdPsi(:,j)=Order1(T(:,j),dPsi);
       d2Td2Psi(:,j)=Order2(T(:,j),dPsi);
       T(2:end,j+1)=T(2:end,j)+dx*(...
           1./(rho(2:end,j).*cp(T(2:end,j)))*dPdx(j)...
           +mu(T(2:end,j)).*rho(2:end,j).*U(2:end,j)./cp(T(2:end,j)).*(dudPsi(2:end,j)).^2 ...
           +1./cp(T(2:end,j)).*(dkRhoUdPsi(2:end,j).*dTdPsi(2:end,j)+k(T(2:end,j)).*rho(2:end,j).*U(2:end,j).*d2Td2Psi(2:end,j))); T(1,j+1)=T(1,j);    
       
       %Calc M and rho
       M(:,j+1)=U(:,j+1)./sqrt(gamma(T(:,j+1))*R.*T(:,j+1));
       rho(:,j+1)=P(j+1)./(R*T(:,j+1));
    
       %Convert Psi back to y
       ruAvg(:,j)=(rho(1:end-1,j).*U(1:end-1,j)+rho(2:end,j).*U(2:end,j))/2;
        for i=2:size(rho,1)
         y_streamlines(i,j)=sum(dPsi(1:i-1)./ruAvg(1:i-1,j));
        end
end
toc

%Get Nozzle Wall streamline
wallidx=find(min(y_streamlines,[],2)<HeightThroat,1,'last');

%% Plot
spacing=1; %Adjust spacing to plot faster by not plotting all x locations
xx=x(1:spacing:size(y_streamlines,2))*100;

%Shift up profiles so nozzle wall is curved and centerline is straight
% YY=flipud(y_streamlines(1:wallidx,1:spacing:end));
YY=(max(y_streamlines(1:wallidx,1:spacing:end))-y_streamlines(1:wallidx,1:spacing:end))*100;
UU=U(1:wallidx,1:spacing:end-1);
TT=T(1:wallidx,1:spacing:end-1);
MM=M(1:wallidx,1:spacing:end-1);

figure(11) %Plot Nozzle Wall
subplot(2,1,2)
p(1)=plot(xx,YY(1,:),'k','Linewidth',2);hold on
% for i=2:wallidx
%     plot(xx,YY(i,:),'k')
%     xlabel('x location [cm]');ylabel('y location [cm]');title('Nozzle streamlines (not spaced out')
% end


figure %Velocity Contour
X=repmat(xx,wallidx,1);
contourf([X;X],[-YY;flipud(YY)],[UU;flipud(UU)],'ShowText','on','LabelSpacing',220);hold on
colormap parula
h2 = colorbar;
set(get(h2,'title'),'string','U [m/s]');
plot(xx,YY(1,:),'k');plot(xx,-YY(1,:),'k');title('Velocity')

figure %Temperature Contour
contourf([X;X],[-YY;flipud(YY)],[TT;flipud(TT)],'ShowText','on','LabelSpacing',300);hold on
colormap parula
h2 = colorbar;
set(get(h2,'title'),'string','T [K]');
plot(xx,YY(1,:),'k');plot(xx,-YY(1,:),'k');title('Temperature')

figure %Mach Contour
Machs=[0.1 0.3 0.5 0.6 0.7 0.8 0.9 1 1.2 1.4 1.5 1.7 2.0 max(max(MM))];
contourf([X;X],[-YY;flipud(YY)],[MM;flipud(MM)],Machs,'ShowText','on','LabelSpacing',220);hold on
colormap parula
h2 = colorbar('Ticks',Machs);
set(get(h2,'title'),'string','Mach');
plot(xx,YY(1,:),'k');plot(xx,-YY(1,:),'k')
title('Mach Number')

% InviscidNozzle=load('InviscidNozzle');
% figure(11) %overlay other nozzles
% if PressureProfile==1
% p(2)=plot(InviscidNozzle.x,InviscidNozzle.Height(InviscidNozzle.M),'r','Linewidth',2);hold on
% set(gca,'Children',[p(1) p(2)])
% legend(p,'Viscous','Inviscid')
% xlabel('x [cm]');ylabel('y [cm]');grid on;
% set(gca,'yMinorTick','on')
% end
% ylim([0 0.06])

% 
% %Compare Inviscid to viscous mach number
% figure
% plot(x*100,M(end,:));hold on
% plot(InviscidNozzle.x,InviscidNozzle.M);
% ylabel('Mach');xlabel('x [cm]');xlim([-0.1 .4])
% legend('Viscous nozzle centerline','Inviscid 1D Nozzle')


%Velocity profile at end
c=hsv(5);
figure
[~,throatidx]=min(YY(1,:));
[x1idx]=find(x>.001,1);
[x2idx]=find(x>.002,1);
plot(UU(:,xin(end)),YY(:,xin(end)),'color',c(1,:),'linewidth',1.5);hold on
plot(UU(:,throatidx),YY(:,throatidx),'color',c(2,:),'linewidth',1.5);
plot(UU(:,x1idx),YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(UU(:,x2idx),YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(UU(:,end),YY(:,end),'color',c(5,:),'linewidth',1.5)
plot(UU(:,xin(end)),-YY(:,xin(end)),'color',c(1,:),'linewidth',1.5)
plot(UU(:,throatidx),-YY(:,throatidx),'color',c(2,:),'linewidth',1.5)
plot(UU(:,x1idx),-YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(UU(:,x2idx),-YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(UU(:,end),-YY(:,end),'color',c(5,:),'linewidth',1.5)
xlabel('Velocity [m/s]');ylabel('y position [cm]')
legend('Inlet x=0','Throat','x=0.1 cm','x=0.2 cm','Exit x=0.4 cm');xlim([0 1.1*max(max(UU))])

%Temperature profile at end
figure
plot(TT(:,xin(end)),YY(:,xin(end)),'color',c(1,:),'linewidth',1.5);hold on
plot(TT(:,throatidx),YY(:,throatidx),'color',c(2,:),'linewidth',1.5);
plot(TT(:,x1idx),YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(TT(:,x2idx),YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(TT(:,end),YY(:,end),'color',c(5,:),'linewidth',1.5)
plot(TT(:,xin(end)),-YY(:,xin(end)),'color',c(1,:),'linewidth',1.5)
plot(TT(:,throatidx),-YY(:,throatidx),'color',c(2,:),'linewidth',1.5)
plot(TT(:,x1idx),-YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(TT(:,x2idx),-YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(TT(:,end),-YY(:,end),'color',c(5,:),'linewidth',1.5)
xlabel('Temperature[K]');ylabel('y position [cm]')
legend('Inlet x=0','Throat','x=0.1 cm','x=0.2 cm','Exit x=0.4 cm');xlim([Twall 1.01*Tinf])

%Mach Profiles
figure
plot(MM(:,xin(end)),YY(:,xin(end)),'color',c(1,:),'linewidth',1.5);hold on
plot(MM(:,throatidx),YY(:,throatidx),'color',c(2,:),'linewidth',1.5);
plot(MM(:,x1idx),YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(MM(:,x2idx),YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(MM(:,end),YY(:,end),'color',c(5,:),'linewidth',1.5)
plot(MM(:,xin(end)),-YY(:,xin(end)),'color',c(1,:),'linewidth',1.5)
plot(MM(:,throatidx),-YY(:,throatidx),'color',c(2,:),'linewidth',1.5)
plot(MM(:,x1idx),-YY(:,x1idx),'color',c(3,:),'linewidth',1.5);
plot(MM(:,x2idx),-YY(:,x2idx),'color',c(4,:),'linewidth',1.5);
plot(MM(:,end),-YY(:,end),'color',c(5,:),'linewidth',1.5)
xlabel('Mach Number');ylabel('y position [cm]');grid on
legend('Inlet x=0','Throat','x=0.1 cm','x=0.2 cm','Exit x=0.4 cm');xlim([0 1.1*max(max(MM))])

% %Plot Streamlines
% figure
% PSI=repmat(Psi(1:wallidx)',size(X,2),1)';
% contourf(X,YY,PSI,0:7e-4:.023,'ShowText','on');hold on
% contourf(X,-YY,PSI,0:7e-4:.023,'ShowText','on');
% colormap gray
% h2 = colorbar;,caxis([0,1e-10])
% set(get(h2,'title'),'string','\psi');

%%
function dadb=Order1(a,db)
    dadb=zeros(length(a),1);
    dadb(1)=(a(2)-a(1))/db(1);
    dadb(end)=(a(end)-a(end-1))/(-db(end));%Notes has it divide by -db but that doesnt make sense
    for i=2:length(a)-1
    dadb(i)=(db(i-1)*a(i+1)-db(i)*a(i-1)-(db(i-1)-db(i))*a(i))/(2*db(i-1)*db(i));
    end
end

function d2ad2b=Order2(a,db)
    d2ad2b=zeros(length(a),1);
    d2ad2b(1)=2*(((db(1)+db(1+1))*a(1+1)-db(1)*a(1+2)-db(1+1)*a(1))/...
        (db(1)^2*(db(1)+db(1+1))-db(1)*(db(1)+db(1+1))^2));
    d2ad2b(end)=2*(((db(end-1)+db(end))*a(end-1)-db(end-1)*a(end-2)-db(end)*a(end))/...
        ((db(end-1)+db(end))*db(end-1)^2-(db(end-1)+db(end))^2*db(end-1)));
    for i=2:length(a)-1
    d2ad2b(i)=2*((db(i)*a(i-1)+db(i-1)*a(i+1)-(db(i-1)+db(i))*a(i)) /...
        (db(i)*db(i-1)^2+db(i)^2*db(i-1)));
    end
end
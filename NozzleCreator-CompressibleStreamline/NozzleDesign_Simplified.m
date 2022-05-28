% Calculate nozzle shape using compressible streamline theory. A boundary layer grows over a flat plate before reaching the nozzle inlet.
clear all; close all;
%If you make dy too small it wont work (like 1/20 is too small)
%% Inputs
HeightThroat=0.04/2 /100; %m
Twall=450;  %K Wall Temperature
Uinf=340; %m/s inlet velocity
Pinf=101325;  %Pa inlet static pressure
Tinf=2000;    %K inlet static temperature

dx=0.000004; 
dy=HeightThroat/10;
y_centerInlet=0.1/100;%m centerline height at the inlet 
    %will be defining nozzle by the throat size so this y is arbitrary. May
    %need to adjust inlet profile if inlet is smaller than the BL, since
    %flowwill be fully developed
 
%Constants
R=287; %J/kgK
Pr=0.74; %Prandtl#   μCp/k
k=0.11; %W/mK
cp=1004; %J/kgK
gamma=1.4;
mu=6.3793e-05;

%% Input Nozzle Pressure Profile
%Calc values
rhoInf=Pinf/(R*Tinf); %kg/m^3 inlet density 
nu=mu/rhoInf;      

%Define grid in (X,Y)
x=[0:dx:0.004];
y=[0:dy:y_centerInlet]';

P_theory=@(x) 1-0.979*exp(-(9.7e-3)./(x*100-0.018));
P_x=[0 0.005 0.00725 0.01 0.0125 0.015 0.0175 0.02 0.0225 0.023]/100;
P_y=[1 0.9951154804 0.9896187008 0.9798615245 0.9677467384...
    0.9519943399 0.9317019646 0.905188353 0.8686136719 0.8593138331];
cf=polyfit(P_x*1000,P_y,4);
top=find(x<=P_x(end),1,'Last');bot=find(x>P_x(end),1);
Pfit=polyval(cf,x([1:top])*1000);idx=find(Pfit>1,1,'last');Pfit(1:idx)=1;
P_P0=([Pfit(1:top),P_theory(x(bot:end))]);

figure;plot(x,P_P0,'k','Linewidth',2);grid on
%% Inlet flat plate calculation
%Calc boundary layer development using Blasius

%Inlet velocity profile after a 0.4cm flat plate
tspan=y*sqrt(Uinf/(nu*0.004));
[eta,f]=ode45(@blasiusode,tspan,[0,0,0.33206]);
U(:,1)=f(:,2)*Uinf;
    
%Inlet temperature profile
tspan=y*Pr^(-1/3)*sqrt(Uinf/(nu*0.004));
[eta,f]=ode45(@blasiusode,tspan,[0,0,0.33206]);
T(:,1)=f(:,2)*(Tinf-Twall)+Twall;

M(:,1)=U(:,1)./sqrt(gamma*R.*T(:,1));
rho(:,1)=Pinf./(R*T(:,1));

%Convert y to Psi
ruAvg=(rho(1:end-1,1).*U(1:end-1,1)+rho(2:end,1).*U(2:end,1))/2;
dPsi=ruAvg.*dy;
Psi=cumsum([0;dPsi]);

%% Nozzle calculation
P=P_P0*Pinf';  %Assuming Pinf=P0 since our P_P0 dist starts at 1.  
              % based off Minf our P_P0=0.9055
dPdx=diff(P)/dx;

for j=1:length(x)-1
     %Calc next U profile
       dudPsi(:,j)=Order1(U(:,j),dPsi);
       d2ud2Psi(:,j)=Order2(U(:,j),dPsi);
       dMuRhoUdPsi(:,j)=Order1(mu*rho(:,j).*U(:,j),dPsi);
       U(2:end,j+1)=U(2:end,j)+dx*(...
           -1./(rho(2:end,j).*U(2:end,j))*dPdx(j)...
           +dMuRhoUdPsi(2:end,j).*dudPsi(2:end,j)...
           +mu*rho(2:end,j).*U(2:end,j).*d2ud2Psi(2:end,j)); U(1,j+1)=U(1,j);      
     %Calc next T profile
       dkRhoUdPsi(:,j)=Order1(k*rho(:,j).*U(:,j),dPsi);
       dTdPsi(:,j)=Order1(T(:,j),dPsi);
       d2Td2Psi(:,j)=Order2(T(:,j),dPsi);
       T(2:end,j+1)=T(2:end,j)+dx*(...
           1./(rho(2:end,j).*cp)*dPdx(j)...
           +mu*rho(2:end,j).*U(2:end,j)./cp.*(dudPsi(2:end,j)).^2 ...
           +1./cp*(dkRhoUdPsi(2:end,j).*dTdPsi(2:end,j)+k*rho(2:end,j).*U(2:end,j).*d2Td2Psi(2:end,j))); T(1,j+1)=T(1,j);   
     %Calc next M and rho
       M(:,j+1)=U(:,j+1)./sqrt(gamma*R.*T(:,j+1));
       rho(:,j+1)=P(j+1)./(R*T(:,j+1));  
     %Convert Psi back to y
       ruAvg(:,j)=(rho(1:end-1,j).*U(1:end-1,j)+rho(2:end,j).*U(2:end,j))/2;
       for i=2:size(rho,1)
         y_streamlines(i,j)=sum(dPsi(1:i-1)./ruAvg(1:i-1,j));
       end
end
%%
function dAdeta = blasiusode(eta,A)
    %Blasius equation 2*A4+A1*A3=0
dAdeta=zeros(3,1);
dAdeta(1)=A(2);
dAdeta(2)=A(3);
dAdeta(3)=-A(1)*A(3)/2;
end
function dadb=Order1(a,db)
    dadb=zeros(length(a),1);
    dadb(1)=(a(2)-a(1))/db(1);
    dadb(end)=(a(end)-a(end-1))/(-db(end));
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
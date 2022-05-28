%Calculate properties in nonsimple region of Expansion Reflection

%update 12/12/21 changed slope inputs to findslope for right running wave
%to be -a becasue it should be finding slope of 1/u+a (now wave are
%diverging after reflection as they should)

%Assume: Calorically Perfect
clear all; close all
syms x
%% Inputs
L = 50;     %m      expansion tube length
T4 = 425;   %K      Driver section static temperature
p4 = 1e6;   %Pa     Driver section static pressure
u3 = 175;    %m/s    Mass motion velocity

n = 10;      %--     Discretation number

%Constants 
R = 287.15;    %J/kgK  Gas Constant
gamma=1.4;  %--     Ratio of heat capacity

%allocate
a=zeros(n); %m/s    speed of sound
u=zeros(n);u(1,[1:n])=linspace(0,u3,n); %m/s velocity
J_plus=zeros(n); %right running wave Reiman Invariant
J_neg=zeros(n); %left running wave Reiman Invariant
xcoord=zeros(n,n+1);
tcoord=zeros(n,n+1);

%% Solve 

%Solve node 1 (first wave to reflection)
a(1,1)=sqrt(gamma*R*T4); %speed of sound in region 4
J_plus(1,1)=u(1,1)+2*a(1,1)/(gamma-1);
J_neg(1,1)=-J_plus(1,1);
xcoord(1,1)=0; %x coordinate of 0 at wall
tcoord(1,1)=-L/(u(1,1)-a(1,1));

%Solve nodes at entrance of wave into nonsimple region
for j=2:n
    a(1,j)=a(1,1)-(gamma-1)/2*u(1,j);
    J_plus(1,j)=u(1,j)+2*a(1,j)/(gamma-1);
    J_neg(1,j)=u(1,j)-2*a(1,j)/(gamma-1);
    
    m02=1/(u(1,j)-a(1,j)); %left
    m12= findSlope(u(1,j-1),-a(1,j-1),u(1,j),-a(1,j)); %right
    xcoord(1,j)=solve(m12*(x-xcoord(1,j-1))+tcoord(1,j-1)==m02*(x-L),x);
    tcoord(1,j)=(xcoord(1,j)-L)*m02;
    if j==n;slopef(1)=m12;end
end

for i=2:n
    %Solve wall points
    u(i,i)=0;
    J_neg(i,i)=J_neg(i-1,i);
    J_plus(i,i)=-J_neg(i,i);
    a(i,i)=(gamma-1)/4*(J_plus(i,i)-J_neg(i,i));
    
    m12=findSlope(u(i-1,i),a(i-1,i),u(i,i),a(i,i));
    xcoord(i,i)=0;
    tcoord(i,i)=m12*(xcoord(i,i)-xcoord(i-1,i))+tcoord(i-1,i);
    
    if i==n;slopef(i)=-m12;end
    
    %Solve other points
    for j=i+1:n
        J_neg(i,j)=J_neg(i-1,j);
        J_plus(i,j)=J_plus(i,j-1);
        u(i,j)=(J_plus(i,j)+J_neg(i,j))/2;
        a(i,j)=(gamma-1)/4*(J_plus(i,j)-J_neg(i,j));
        
        m12=findSlope(u(i,j-1),-a(i,j-1),u(i,j),-a(i,j));  %right
        m02=findSlope(u(i-1,j),a(i-1,j),u(i,j),a(i,j)); %left
        xcoord(i,j)=solve(m02*(x-xcoord(i-1,j))+tcoord(i-1,j)... 
            ==m12*(x-xcoord(i,j-1))+tcoord(i,j-1),x);
        tcoord(i,j)=m02*(xcoord(i,j)-xcoord(i-1,j))+tcoord(i-1,j);
        
        if j==n;slopef(i)=m12;end
       
    end
end

%% Plot Lines

ts=0;xs=L; %Starting point
plot(xcoord,tcoord,'ko');xlabel('x(m)');ylabel('t(sec)')

%initial slopes
for i=1:n
    slope0(i)=1/(u(1,i)-a(1,i));
end

%Front of wave gets to diaphram at...
tcoord(:,n+1)=.01+slopef(1)*(L-xcoord(1,n))+tcoord(1,n); %%%%%%%%%%%%%%%%%% added .01 a bit to make it go past diaphram
xcoord(:,n+1)=(tcoord(:,n+1)-tcoord(:,n))./slopef'+xcoord(:,n);

%Make vectors of points following each wave
for i=1:n
    for j=i:n
        line(j,1,i)=xcoord(i,j); %line(index,x or t,line#)
        line(j,2,i)=tcoord(i,j);
        line(i,1,j)=xcoord(i,j);
        line(i,2,j)=tcoord(i,j);
    end
end

%add starting and end point and plot
hold on
for i=1:n%1:round(n/4):n
    line(n+1,:,i)=[xcoord(i,n+1),tcoord(i,n+1)];
    Line(:,:,i)=[[L,0];line(:,:,i)];
    plot(Line(:,1,i),Line(:,2,i),'linewidth',2)
end


%% Plot Pressure at wall vs Time

% P=p4*(1-(gamma-1)/2*(diag(u)./diag(a))).^(2*gamma/(gamma-1));
c=p4/a(1,1)^(2*gamma/(gamma-1));
P=c*diag(a).^(2*gamma/(gamma-1));
figure()
plot([0,diag(tcoord)'],[1,P'*10^-6])
title('Pressure at Wall vs Time')
xlabel('Time (s)')
ylabel('Pressure (MPa)')


function dt_dx = findSlope(u1,a1,u2,a2)
    %Finds left running slope
    dt_dx=tan((atan(1/(u1-a1))+atan(1/(u2-a2)))/2);
end






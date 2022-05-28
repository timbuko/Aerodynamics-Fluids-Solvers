%Solve for properties in a variable area duct with heating

%% Set up solve
dx=0.00001; %step size m

%Given
gamma =1.4;     %        heat capacity ratio
R     =287;     %J/kgK   gas constant
Cp    =1004.5;   %J/kgK   Specific Heat (constant P)
M(1)  =6;       %        Mach at inlet
T0(1)=600;      %K       Total Temperature at inlet
P0(1)=200000;   %Pa      Total Pressure at inlet

%Convert P0 to P
P(1)=P0(1)*(1+(gamma-1)*M(1)^2/2)^(gamma/(gamma-1))^-1;

%Convert T0 to T
T(1)=T0(1)*(1+(gamma-1)/2*M(1)^2)^-1;

%Inlet Conditions
x(1)    =   0;                      %m     Axial distance 
r(1)    =   -0.01*x(1)+1.1;         %m     Radius of tube
u(1)    =   M(1)*sqrt(gamma*R*T(1));%m/s   Flow velocity
A(1)    =   pi*(r(1))^2;            %m^2   Area
rho(1)  =   P(1)/(R*T(1));          %kg/m^3 Density
mdot    =   rho(1)*u(1)*A(1)        %kg/s   Mass flow rate

QQ(1)   =   10000*sin(x(1));        %W/m^2 Heat flux into the duct
Q(1)    =   2*pi*r(1)*QQ(1);        %W/m Heat per m in axial direction
                           % So to get total heat transfer in the pipe
                           % integrate over the length of pipe (=sum(Q*dx))
%% Solve where current value = value plus current slope
i=1;
while M(i)>1
    
   if i>10E8
    disp('Did not solve')
    break
   end
       
    x(i+1)=x(i)+dx;
    r(i+1)=-0.01*x(i+1)+1.1;
    QQ(i+1)=10000*sin(x(i+1));
    
    A(i+1) = pi*(r(i+1))^2;
    dA(i) = pi*(r(i+1)^2-r(i)^2);
    
    dQ(i)=(QQ(i+1)*pi*r(i+1)+QQ(i)*pi*r(i))*dx;
    
    dT0(i) = dQ(i)/(mdot*Cp);
    T0(i+1)=T0(i)+dT0(i);
    
    dM(i) = M(i)*(...
        (1+(gamma-1)*M(i)^2/2)/(M(i)^2-1)*dA(i)/A(i)...
        -(gamma*M(i)^2+1)/(2*(M(i)^2-1))*dQ(i)/(mdot*Cp*T(i))...
        );  
    M(i+1)=M(i)+dM(i);
    
    T(i+1)=T0(i+1)*(1+(gamma-1)/2*M(i+1)^2)^-1; 
    dT(i)=T(i+1)-T(i);
    
    du(i)= u(i)*(dM(i)/M(i)+0.5*dT(i)/T(i));
    u(i+1)=u(i)+du(i);
    
    drho(i)=rho(i)*(-dA(i)/A(i)-du(i)/u(i));
    rho(i+1)=rho(i)+drho(i);
    
    dp(i)= P(i)*(drho(i)/rho(i)+dT(i)/T(i));
    P(i+1) = P(i)+dp(i);
    
    
     i=i+1;
end
%Chop off last point cuz i+1 point might be wack
M(end)=[];T(end)=[];P(end)=[];rho(end)=[];x(end)=[];A(end)=[];r(end)=[];QQ(end)=[];T0(end)=[];u(end)=[];


Mach = M(end)
L=x(end) 

%get P0 and T0
P0=P.*(1+(gamma-1)./2.*M.^2).^(gamma/(gamma-1));

%Convert T0 to T
T0=T.*(1+(gamma-1)./2.*M.^2);


%% Find time in duct
t=L/mean(u)

%% Static Pressure at choke
StaticP=P(end)

%% Static Temp at choke
StaticT=T(end)

%% Plot

figure()
plot(x,M),title('M vs position'),xlabel('Axial Position (m)'),ylabel('Mach')

figure()
plot(x,P./1000),title('P-static vs position'),xlabel('Axial Position (m)'),ylabel('Static Pressure (kPa)')

figure()
plot(x,P0./1000),title('P-total vs position'),xlabel('Axial Position (m)'),ylabel('Total Pressure (kPa)')

figure()
plot(x,T),title('T-Static vs position'),xlabel('Axial Position (m)'),ylabel('Static Temp (K)')

figure()
plot(x,T0),title('T-Total vs position'),xlabel('Axial Position (m)'),ylabel('Total Temp (T)'),ylim([0,T0(1)*1.1])


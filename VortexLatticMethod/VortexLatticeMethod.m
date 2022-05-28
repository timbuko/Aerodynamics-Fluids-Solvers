%Vortex Lattice Method uses wing geometry to calculate the coefficient of
%lift. It does this by calculating gamma for a bunch of pannels along the
%span of the wing.

%This method improves upon Prandtl Lifting Line and allows for taper,
%sweep

%% Inputs
syms y

N      = 20;         %Number of Panels on one half of wing
b      = 6;         %Wing Span
cr     = 1;          %Wing root chordlength
lambda = 1;          %taper ratio lambda=ct/cr
sweep  = 20*pi/180; 
vinf   = 10;         %Bulk wind velocity
alfa   = 2*pi/180;   %angle of attack
yaw    = 0;          %yaw angle
rho    = 1.225;      %Density of Air 
s=cr*(1+lambda)*b/2; %Planform Area

vmat=[vinf*cos(alfa),vinf*sin(yaw),vinf*sin(alfa)];


%% Define Equations and Mesh

A=[b]; %initialize with sym so you can add syms later
B=[b];
K=[b];

%Position equations
xle=abs(y)*tan(sweep); %Leading edge x position
c=cr+2/b*cr*(lambda-1)*abs(y); %chord length
xqrt=xle+c/4; %Quarter chord x position
x3qrt=xle+3*c/4; %Three quarter chord x posistion

for n=1:2*N %Get coordinates for all horseshoe filament vertices and BCs
   A(n,2)=b*(n-1)/(2*N)-b/2;
   A(n,1)=subs(xqrt,y,A(n,2));
   B(n,2)=b*(n)/(2*N)-b/2;
   B(n,1)=subs(xqrt,y,B(n,2));
   K(n,2)=(A(n,2)+B(n,2))/2;
   K(n,1)=subs(x3qrt,y,K(n,2));
   A(n,3)=0;B(n,3)=0;K(n,3)=0;
   
end

%Plot filaments vertices and boundary points
figure(1)
scatter(A(:,1),A(:,2))
hold on
scatter(B(:,1),B(:,2)),scatter(K(:,1),K(:,2))
fplot(finverse(xle),sort([0,b/2*tan(sweep)]),'k')
fplot(-finverse(xle),sort([0,b/2*tan(sweep)]),'k')
fplot(finverse(xle+c),sort([cr,b/2*tan(sweep)+cr*lambda]),'k')
fplot(-finverse(xle+c),sort([cr,b/2*tan(sweep)+cr*lambda]),'k')
legend('A','B','K','Wing Edges')
title('Wing Geometry'),xlabel('x-position'),ylabel('y-position')

%% Solve for Ciruclation Gamma, and Lift Force

%Create Matrix to Solve for Gamma
hmat=ones(N,N);
for i=1:N
    for j=1:N
        hmat(i,j)=dot(horseshoe(A(j,:),B(j,:),K(i,:))...
            +horseshoe(A(2*N+1-j,:),B(2*N+1-j,:),K(i,:)),[0,0,1]);
    end
end                                    

%Solve For Gamma

G= hmat^-1*(dot(-vmat,[0,0,1]).*ones(N,1));

%Calculate Cl
CL=4*sum(G)/(N*vinf*cr*(lambda+1))
L=CL/2*rho*vinf^2*s

%Plot Gamma vs span(y)
ybodypoints=transpose(K(:,2));
gamgam= [transpose(G),fliplr(transpose(G))];

figure(2)
plot(ybodypoints,gamgam),title('Gamma vs span'),xlabel('Span'),ylabel('Gamma')


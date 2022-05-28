%Solve and plot for 1D pressure driven flow between two plates

close all
clear all
%% Define Variables
rho =  1261.0;  %kg/m^3
nu  =  0.001;   %m^2/s
dp  =  -126.1;  %N/m
h   =  0.12;    %m

N   =  25;       %mesh density
dy=h/(N-1); 

tt = 1;         %time range
dt = .2;       %time step
u= zeros(N,1); %initial condition

solveConvergence = 1; %Solve for Convergence or given time length?

    conv = 1e-10;  %Required difference between to points to be considered 
                   %converged of |u(n+1)-u(n)|

%% Create Matrix of Equations

%Constants
A=nu/(2*dy^2);
B=1/dt+nu/dy^2;
C=1/dt-nu/dy^2; 


% lmat*u^(n+1)= rmat*u^n + pmat
lmat=zeros(N);
rmat=zeros(N);
pmat=zeros(N,1);


lmat(1,1)=1;
lmat(N,N)=1;

for j=2:N-1
    pmat(j)= -dp/rho;
    lmat(j,j-1)=-A;
    lmat(j,j)= B;
    lmat(j,j+1)=-A;
    rmat(j,j-1)=A;
    rmat(j,j)= C;
    rmat(j,j+1)=A;
end


%% Solve Equations by looping through time
%u^(n+1)=rmat^-1*lmat*u^n+pmat

if solveConvergence~=1
%Solve for given number of time steps
convergence=1;
for i=1:tt/dt
    u(:,i+1)=lmat^-1*(rmat*u(:,i)+pmat);
end
end


if solveConvergence==1
%Solve to convergence
convergence=1;
i=1;
u(:,i+1)=lmat^-1*(rmat*u(:,i)+pmat);

while any(abs(u(:,i+1)-u(:,i))>conv)
    i=i+1;
    iteration=i
    u(:,i+1)=lmat^-1*(rmat*u(:,i)+pmat);
    if i>5000 & all(abs(u(:,1000)-u(:,900))<abs(u(:,i)-u(:,i-100)))
        convergence=0;
        disp('Does not converge')
        break
    end
end
end


%% Plot Velocity Profile
y=linspace(0,h,N)';

%Plot Flow Movie
if convergence==1
disp('Plotting.....')
[~,tmstp]=size(u);

c = parula(30); %For colorgradient with time

for i=1:30
    p=plot(u(:,i),y,'Color',c(1+end-i,:));
    text(0,h/2,['Time=',int2str(i*dt),'sec'])
    title('Velocity Profile'),ylabel('y position (m)'),xlabel('Velocity (m/s)')
    ylim([0,h]),xlim([0,max(max(u))])
    M(i)=getframe;
    hold on

end

%Plot analytical solution
syms q w
pp=fplot(solve(1/(rho*nu)*dp*(q^2/2-h*q/2)==w,q),'k*-');
legend([p(end),pp(end)],'Numerical','Analytical')

% %Play Movie
% figure(2)
% gif_length = 3; %seconds
% rep = 2;        %how many times movie reps
% movie(M,rep,tmstp/gif_length)

end
%Potential Flow Solver for flow field in a duct with a step

%Assume: Steady, incompressible, inviscid

clear all; close all
%% Input

d=.01; %Discretization step size
conv=10^-6; %Convergence factor

rho=1000; %kg/m^2 fluid density 

Q=6; %m^3/s Flowrate into duct

use_sor=1; %set to 1 to use successive over relaxation (SOR)
beta=1.99;  %SOR condition

%Discretize Field
N=round(2/d);
M=round(6/d);

%initialize
psi=zeros(N,M);
psi_old=psi+1;
cnt=0;

%% Calc Psi
while (sum(sum(abs(psi_old-psi)))>conv)
cnt=cnt+1;

if mod(cnt,100)==1
fprintf('Count: %d\tConv:%f\n',cnt,sum(sum(abs(psi_old-psi))))
end

% if cnt>10^5; break; end

psi_old=psi;

%First
for j=1:floor(M/2)
    %Wall BC
    psi(1,j)=0;
    psi(N,j)=Q;
    
    %Solve Middle Region
    if j>1 && j<floor(M/2)
        for i=2:N-1
            if use_sor==1
                psi(i,j)=(1-beta)*psi(i,j)+beta/4*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1));
            else
                psi(i,j)=1/4*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1));
            end
        end
    end
    
    %Entrance BC
    if j<4
        for i=2:N-1
            psi(i,1)=1/3*(4*psi(i,2)-psi(i,3));
        end
    end
end

%Wall
for i=1:floor(N/2)
    psi(i,floor(M/2))=0;
end

%Second Half
for j=floor(M/2):M
    %Wall BC
    psi(floor(N/2),j)=0;
    psi(N,j)=Q;
    
    
    %Solve Middle Region
    if j<M
        for i=floor(N/2)+1:N-1
            if use_sor==1
                psi(i,j)=(1-beta)*psi(i,j)+beta/4*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1));
            else
                psi(i,j)=1/4*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1));
            end
        end
    end
    
    %Exit BC
    if j>M-4
        for i=floor(N/2)+1:N-1
            psi(i,M)=1/3*(4*psi(i,M-1)-psi(i,M-2));
        end
    end
end


            
end
fprintf('Count: %d\tConv:%f\n',cnt,sum(sum(abs(psi_old-psi))))

%% Velocity Field

%First
for j=1:floor(M/2)
    v(1,j)=0; %velocities on top and bottom wall
    v(N,j)=0;
    u(1,j)=0;
    u(N,j)=0;
    
        for i=2:N-1
            u(i,j)=(psi(i+1,j)-psi(i-1,j))/(2*d);
            if j>1 && j<floor(M/2)
                v(i,j)=-(psi(i,j+1)-psi(i,j-1))/(2*d);   
            end
            
        end

    
end

%Wall
for i=1:floor(N/2)
    u(i,floor(M/2))=0;
    v(i,floor(M/2))=0;
end

%Second Half
for j=floor(M/2):M
    u(floor(N/2),j)=0;%velocities on top and bottom wall
    u(N,j)=0;
    v(floor(N/2),j)=0;
    v(N,j)=0;
    

        for i=floor(N/2)+1:N-1
            u(i,j)=(psi(i+1,j)-psi(i-1,j))/(2*d);
            if j<M
                v(i,j)=-(psi(i,j+1)-psi(i,j-1))/(2*d);
            end
                
        end
    
    
end

%% Calc Pressure 

P=.5*rho*(u.^2+v.^2);

%% Plot

%Streamlines
figure; hold on
contourf(psi,20); colorbar;title('\Psi - Streamlines');
contour(psi,'r','LevelList',0);
contour(psi,'r','LevelList',Q);
axis equal

%Velocity Field
figure()
uu=downsample(u,N/10);uu=downsample(uu',round(N/10))';
vv=downsample(v,round(N/10));vv=downsample(vv',round(N/10))';
quiver(uu,vv);title('Velocity Field')
axis equal

%Pressure
figure()
contourf(P,100); colorbar;title('Pressure Field'); set(gca,'ColorScale','log')
axis equal

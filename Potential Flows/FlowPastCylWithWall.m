%This script plots the velocity field of flow around a cylinder in a wind
%tunnel

clear
syms  Uinf  Vinf   x  y gamma lamda y0 x0 K;

%Graph Bounds
xbound=12;
ybound=6;

%% Define Elementary Stream Functions
psi_free = Uinf*y - Vinf*x;                             %freestream
psi_vort = gamma/2/pi*log(sqrt((x-x0)^2 + (y-y0)^2 ));  %vortex
psi_sourcesink = lamda*atan2((y-y0),(x-x0))/(2*pi());   %source/sink
psi_doublet = -K*(y-y0)/(2*pi()*((x-x0)^2+(y-y0)^2));   %doublet


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forming a wall with two vortices/sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% lam=20;
% gam=20;
% 
% psi_source=subs(psi_sourcesink,[lamda,y0,x0],[lam,0,-2]);
% psi_source2= subs(psi_sourcesink, [lamda,y0,x0],[lam,0,2]);
% psi_vort1=subs(psi_vort,[gamma,x0,y0],[gam,0,2]);
% psi_vort2=subs(psi_vort,[gamma,x0,y0],[-gam,0,-2]);
% 
% psi_total1= psi_source+psi_source2;
% psi_total2= psi_vort1+psi_vort2;
% 
% figure(1)
% ff=fcontour(psi_total1,[-xbound  xbound  -ybound  ybound])
% ff.LevelList = [-20:.5:20]
% figure(2)
% ffv=fcontour(psi_total2,[-xbound  xbound  -ybound  ybound])
% ffv.LevelList = [-10:.2:10]
% 
% %Check if meeting in sources plot is impermeable
% u=diff(psi_total1,y);
% v=-diff(psi_total1,x);
% u2=diff(psi_total2,y);
% v2=-diff(psi_total2,x);
% 
% u_origin=subs(u,[x,y],[0,1]);
% v_origin=subs(v,[x,y],[0,1]);
% 
% u_check=subs(u2,[x,y],[0,0]);
% v_check=subs(v2,[x,y],[0,0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Flow around a cylinder in a tunnel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=3;
h=6;
uinf=1;
vinf=0;
k=2*pi*uinf*R^2;

psi_doub=subs(psi_doublet,[K,x0,y0], [2*pi*uinf*R^2,0,0]);
psi_free1=subs(psi_free,[Uinf,Vinf],[uinf,vinf]);
for i=1:5 %1/2 of the number of vortices
    psi_doubWall(1,i)=subs(psi_doublet,[K,x0,y0],[k,0,i*2*h]);
    psi_doubWall(2,i)=subs(psi_doublet,[K,x0,y0],[k,0,i*2*-h]);
end
psi_walls=sum(sum(psi_doubWall));
psi_total3=psi_doub+psi_free1+psi_walls;

F=fcontour(psi_total3,[-xbound  xbound  -ybound  ybound]);
F.LevelList = [-100:.5:100];
axis equal

 
%Find velocity
u=diff(psi_total3,y);
v=-diff(psi_total3,x);

%Check for wall around cylinder
Ur= u*x/R+v*y/R;
theta= pi;
Ur_check=vpa(subs(Ur,[x,y],[R*cos(theta),R*sin(theta)]),4)
% Ur_check1=vpa(subs(Ur,[x,y],[R*cos(pi/2+theta),R*sin(pi/2+theta)]),4)
% u_cyl=vpa(subs(u,[x,y],[R*cos(theta),R*sin(theta)]),4)
% u_cyl1=vpa(subs(u,[x,y],[R*cos(pi/2+theta),R*sin(pi/2+theta)]),4)
% v_cyl=vpa(subs(v,[x,y],[R*cos(theta),R*sin(theta)]),4)
% v_cyl1=vpa(subs(v,[x,y],[R*cos(pi/2+theta),R*sin(pi/2+theta)]),4)

%Check for tunnel wall below
xx=5;
vtop_check1=eval(subs(v,[x,y],[xx,h]))
vbottom_check1=eval(subs(v,[x,y],[xx,-h]))
% vbottom_check2=eval(subs(v,[x,y],[-xx,-h]))
% vbottom_check3=eval(subs(v,[x,y],[xx^2,-h]))
% ubottom_check1=eval(subs(u,[x,y],[xx^2,-h]))
% ubottom_check2=eval(subs(u,[x,y],[xx^2,-h]))
% ubottom_check3=eval(subs(u,[x,y],[xx^2,-h]))

%the number of vortices after like 5 doesn't change the v velocity on the
%wall and it remains not quite zero. The velocity across the surface of the
%cylinder eventually reads NaN


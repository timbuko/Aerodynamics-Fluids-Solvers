
syms  Uinf  Vinf   x  y gamma lamda y0 x0 k;

% Elementary Stream Functions
psi_free = Uinf*y - Vinf *x;                            %freestream
psi_vort = gamma/2/pi*log(sqrt((x-x0)^2 + (y-y0)^2 ));  %vortex
psi_sourcesink = lamda*atan2((y-y0),(x-x0))/(2*pi());   %source/sink
psi_doublet = -k*(y-y0)/(2*pi()*((x-x0)^2+(y-y0)^2));   %doublet
psi_corner = Uinf*sqrt(x^2+y^2)^n*sin(n*atan2((y),(x))); %corner n=pi/theta


% if you want to substitute in values for Uinf and Vinf
psi_free2 = subs(psi_free,[Uinf,Vinf],[10,1]);

% create contour plot of the streamfunction to visualize streamlines
% make all the lines the same color - black
% figure(1)
% fcontour(psi_free2, [-6  6  -4  4],'k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stream function for a vortex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms gamma

% psi_vort = gamma/2/pi*log(sqrt(x^2 + y^2 ));
% psi_vort2 = subs(psi_vort,gamma, 1);

% figure(2)
% % fcontour(psi_vort2,[-5 5 -5 5]);
% axis equal;

%%%%%%%%%%%%%%%%%%
% Doublet
%%%%%%%%%%%%%%%%%%%%%

psi_doublet = -k*(y-y0)/(2*pi()*((x-x0)^2+(y-y0)^2));
% psi_doublet1 = subs(psi_doublet,[k,x0,y0], [20,0,0]);
% psi_free4 = subs(psi_free,[Uinf,Vinf],[5,0]);
% fcontour([psi_doublet1+psi_free4], [-6  6  -4  4])
% axis equal;
% 
% u=diff(psi_doublet1+psi_free4,x)
% v=diff(psi_doublet1+psi_free4,y)
% 
% solve([u==0,y==0],[stag(1),stag(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Flow past Rankine Oval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms lamda y0 x0 k
% 
% tit= 'Rankine Oval';
% xbound=6;
% ybound=4;
% num=10; %how many surface points 
% R=2.7; %based on stagnation points (doesnt actually make radius R)
% 
% %Create Psi equation and contour plot
% psi_sourcesink = lamda*atan2((y-y0),(x-x0))/(2*pi());
% psi_source = subs(psi_sourcesink,[lamda,x0,y0], [10,-2,0]);
% psi_sink = subs(psi_sourcesink,[lamda,x0,y0], [-10,2,0]);
% psi_free3 = subs(psi_free,[Uinf,Vinf],[2,0]);
% 
% psi_total=(psi_source+psi_sink+psi_free3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Flow past a cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tit = 'Flow past cylinder';
% xbound=6;
% ybound=5;
% num=10; %how many surface points 
% 
% U__inf= 1;
% R= 2;
% 
% Create Psi equation
% psi_free5= subs(psi_free,[Uinf,Vinf],[U__inf,0]);
% psi_doublet2 = subs(psi_doublet,[k,x0,y0], [2*pi*U__inf*R^2,0,0]);
% psi_total= psi_free5+psi_doublet2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Flow around  rotating cylinder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tit = 'Flow past rotating cylinder';
xbound=10;
ybound=10;
num=10; %how many surface points 

U__inf= 1;
R= 2;
gam = 4*U__inf*pi*R;

%surface values for R=2 U=1
ysurf=[ 0,-1.2571,-1.6630,-1.8856,-1.9876,-1.9876,-1.8856,-1.6630,-1.2571,0,0,1.2571,...
 1.6630,1.8856,1.9876,1.9876,1.8856,1.6630,1.2571,0];
xvals=[-2.0000,-1.5556,-1.1111,-0.6667,-0.2222,0.2222,0.6667,1.1111,1.5556,2.0000,-2.0000,-1.5556...
   -1.1111,-0.6667   -0.2222    0.2222    0.6667    1.1111    1.5556    2.0000];

%Define Psi
psi_free6= subs(psi_free,[Uinf,Vinf],[U__inf,0]);
psi_doublet3 = subs(psi_doublet,[k,x0,y0], [2*pi*U__inf*R^2,0,0]);
psi_vort3 = subs(psi_vort,[gamma,x0,y0], [gam,0,0]);
psi_total= psi_free6+psi_doublet3+psi_vort3;





%Plot Psi graph   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
fcontour(psi_total, [-xbound  xbound  -ybound  ybound],'LevelStep',1),title(tit) %LevelStep Changes number of lines on plot

hold on
axis equal
%Find Stagnation points and streamline
u_vel=diff(psi_total,y);
v_vel=-diff(psi_total,x);
[xstag,ystag] = vpasolve([u_vel==0,v_vel==0],[x,y]);
%Plot stagnation streamline
stag_psi = subs(psi_total,[x,y],[min(xstag),min(ystag)]);
fcontour(psi_total,[-xbound  xbound  -ybound  ybound],'r','LevelList',[stag_psi],'LineWidth',2)
hold off
% %Find Surface Points   %%%%%%%%%
% xvals= [linspace(-R,R,num),linspace(-R,R,num)];
% ysurf = zeros(1,2*num);
% for i= 1:num
% 
%     psi = subs(psi_total,x,xvals(i));
%         temp = vpasolve(psi-stag_psi,y,3);
%         j=0;
%         temp2=temp(1);
%         for k=1:length(temp)  %Sort out solutions and choose only the real solution
%             if isreal(temp(k))
%                 j=j+1;
%                 temp2(j)=temp(k);
%             end
%         end
%         
%         temp3=temp2(1);  %sometimes there is a pos, neg, and zero solution, so this just pics first one
%         
%         ysurf(i) = temp3;
%         ysurf(i+num)=-temp3;
%     
%         
% end

%Find Velocity vectors at surface points
Vel = zeros(length(xvals),2);
for i=1:length(xvals)
    Vel(i,1)= subs(u_vel,[x,y],[xvals(i),ysurf(i)]);
    Vel(i,2)= subs(v_vel,[x,y],[xvals(i),ysurf(i)]);
end
%Calc Coefficient of Pressure   (bottom then top)
uinf=subs(u_vel,[x,y],[-10,0]);
vinf=subs(v_vel,[x,y],[-10,0]);
V_inf= sqrt(uinf^2+vinf^2);
for i=1:length(xvals)
    V_mag(i)=sqrt(Vel(i,1)^2+Vel(i,2)^2);
end
for i=1:length(xvals)
    Cp(i) = 1 - V_mag(i)^2/V_inf^2;
end
subplot(1,2,2)
plot(xvals(1:num),-Cp(1:num)), title('Cp vs x'), xlabel('x position'), ylabel('-Cp')
hold on
plot(xvals(num+1:2*num),-Cp(num+1:2*num)), legend('Bottom','Top','Location','South')
hold off

% %Calc Cl and Cd  for c=1
% for i=1:2*num
%     Fy(i)= -Cp(i)*2*pi*R/length(Cp)*sin(atan2(ysurf(i),xvals(i)));
%     Fx(i)=-Cp(i)*2*pi*R/length(Cp)*cos(atan2(ysurf(i),xvals(i)));
% end
% Cl=vpa(sum(Fy),5)
% Cd=vpa(sum(Fx),5)

%Method 2 for calc Cl
for i=2:num-1
    Fy2(i)=-Cp(i)/(2*R)*abs(xvals(i-1)-xvals(i));
    Fy2(2*i)=-Cp(2*i)/(2*R)*abs(xvals(i-1)-xvals(i));
end
Cl=vpa(sum(Fy2),5)

%%%%%%%%%%%%%%%%%%%%
% problem 3
%%%%%%%%%%%%%%%%%%%

% psi_3= 0.5*x^2*y^2;
% fcontour([psi_3], [-10  10  -10  10])

%%%%%%%%%%%%%%5
% problem 7
%%%%%%%%%%%%%%%5

% psi_7= atan(y,x)-atan(x,y);
% fcontour([psi_7], [-10  10  -10  10])
% Calculate lift and Drag for supersonic flow past a wedge

%% Inputs
M_inf   =    2.0; %Mach#
p_inf   =    100000 ; %Pa  Upstream static pressure
alpha   =    5  ; %deg  Angle of attack
theta_f =    4  ; %deg  Forward wedge half angle
theta_a =    4  ; %deg  Aft wedge half angle
chord   =     1 ; %  Chord length
gamma   =    1.4; %specific heat ratio

%% Plot wing section

x_p = chord*tand(theta_a)/(tand(theta_f)+tand(theta_a));
y_p = x_p*tand(theta_f);
x= [0,0;x_p,y_p;chord,0;x_p,-y_p;0,0];

plot(x(:,1),x(:,2))
text(x_p/2,y_p/2,'1'),text(x_p/2,-y_p/2,'2'),text(3*x_p/2,y_p/2,'3'),text(3*x_p/2,-y_p/2,'4') %label surfaces
axis equal

%% Calculate Surface Pressures 

%convert angles to radians
theta_f=theta_f*pi/180;
theta_a=theta_a*pi/180;
alpha=alpha*pi/180;


%pass alpha,theta,M,P,gamma
%Front top surface (1)
if alpha<theta_f
    if alpha>0 
        theta=alpha-theta_f;
    else
        theta=-alpha+theta_f; 
    end
    beta1=findbet(M_inf,theta,gamma);
    mn_inftop=M_inf*sin(beta1);
    p1_pinf=presratoblique(mn_inftop,gamma)
    p_1=p_inf*p1_pinf;
    
    mn1=oblique(mn_inftop,gamma);
    M1=mn1/sin(beta1-theta)
elseif alpha>theta_f
    theta=alpha-theta_f;
    pmanginf = pmang(M_inf,gamma);
    pmang1=pmanginf+theta;
    M1=findm2pmgen(pmang1,gamma);
    p1_pinf=presratpm(M_inf,M1,gamma);
    p_1=p_inf*p1_pinf;
else
    M1=M_inf;
    p_1=p_inf;
end

%Front bottom surface (2)
if alpha>-theta_f
    if alpha<0 
        theta2=-alpha-theta_f;
    else
        theta2=alpha+theta_f; 
    end
    beta2=findbet(M_inf,theta2,gamma);
    mn_infbot=M_inf*sin(beta2);
    p2_pinf=presratoblique(mn_infbot,gamma);
    p_2=p_inf*p2_pinf;
    
    mn2=oblique(mn_infbot,gamma);
    M2=mn2/sin(beta2-theta2);
elseif alpha<-theta_f
    theta2=-alpha-theta_f;
    pmanginf = pmang(M_inf,gamma);
    pmang2=pmanginf+theta2;
    M2=findm2pmgen(pmang2,gamma);
    p2_pinf=presratpm(M_inf,M2,gamma);
    p_2=p_inf*p2_pinf;
else
    M2=M_inf;
    p_2=p_inf;
end

%Back top surface (3)

theta3=abs(theta_f-(-theta_a));
pmang1 = pmang(M1,gamma);
pmang3=pmang1+theta3;
M3=findm2pmgen(pmang3,gamma);
p3_p1=presratpm(M1,M3,gamma);
p_3=p_1*p3_p1;

%Back bottom surface (4)

theta4=abs(theta_f-(-theta_a));
pmang2 = pmang(M2,gamma);
pmang4=pmang2+theta4;
M4=findm2pmgen(pmang4,gamma);
p4_p2=presratpm(M2,M4,gamma);
p_4=p_2*p4_p2;


%% Find Forces

Fx=p_1*y_p-p_3*y_p+p_2*y_p-p_4*y_p
Fy=p_2*x_p-p_1*x_p+p_4*(chord-x_p)-p_3*(chord-x_p)

L=Fy*cos(alpha)-Fx*sin(alpha)
D=Fx*cos(alpha)+Fy*sin(alpha)

Cl=L/(0.5*gamma*p_inf*M_inf^2*chord)
Cd=D/(0.5*gamma*p_inf*M_inf^2*chord)
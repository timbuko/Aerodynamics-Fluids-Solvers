%Aka Priliminary Method: This method uses the Prandtl Lifting line eqn to solve for the Gamma dist.
%over the span of a wing

%this method is best for high Aspect ration, rectangular flat wings
syms A1 A3 A5 A7 alpha o y


%% Inputs

%Set number of panels
theta=[pi/4,pi/2,pi/8,3*pi/8];
n=[1,2,3,4];

%Parameters of wing as a function of theta
%  calculated as linear gradients between root and tip airfoils


cr=1;                                       %Root Chord length
lambda=0.6;                                 %Taper ratio
AR=7;                                       %Aspect ratio
Vinf=  300000/3600;                         %Velocity 
alfa= 1.8    *pi/180;                       %Angle of attack

root.alfL = -1.2    *pi/180;                 %Zero lift AoA at root (converted to rad)
root.a0 =   0.1     *180/pi;                 %Airfoil Lift curve slope (converted to 1/rad)
root.alfa = 0       *pi/180;                 %Twist at root (converted to rad)

tip.alfL = -2.0    *pi/180;              
tip.a0 =   0.098   *180/pi;                 
tip.alfa = 3       *pi/180;                 
b=AR*cr*(1+lambda)/2;                       %Span

%get quantities as function of theta
alfL= (tip.alfL-root.alfL)*cos(o)+root.alfL;       
a0=(tip.a0-root.a0)*cos(o)+root.a0;     
alf=alpha-((tip.alfa-root.alfa)*cos(o)+root.alfa); %Angle of attack seen by airfoil
c=(lambda-1)*cr*abs(cos(o))+cr;                    %Chord length



%% Solve For A's

%Set matrix of A's and n theta's to solve for A's
for i= 1:4 %theta
    for j=1:4 %n
        a0i=subs(a0,o,theta(i));
        ci=subs(c,o,theta(i));
        eq(i,j)=sin((2*n(j)-1)*theta(i))*(8*(b/2)/(a0i*ci)+(2*n(j)-1)/sin(theta(i)));
        
    end
    %Set matrix of alpha-alpha_L0's to solve for A's
    alfL1(i)=subs(alfL,o,theta(i));
    alf1(i)=subs(alf,o,theta(i));
    
end
        
eqinv=eq^-1;
alff=(alf1-alfL1);
%Solve for A's
A=vpa(eqinv*transpose(alff),4)

%% Solve for Gamma distribution over span of wing

%Substitute alpha into A's 
a=subs(A,alpha,alfa);

%Summation piece of gamma equation
stuff=0;
alstuff=0;
for i=1:4
    stuff=stuff+a(i)*sin((2*i-1)*o);     %Summation with value alpha
    alstuff=alstuff+A(i)*sin((2*i-1)*o); %Summation with symbolic alpha
end

%Gamma distribution
gamma=vpa(4*b/2*Vinf*stuff,4);          %function of theta
gamVsSpan(y)=subs(gamma,o,acos(2*y/b)); %function of span y
subplot(1,2,1)
fplot(gamVsSpan(y),[-b/2,b/2]),title('Gamma vs span'),xlabel('Span'),ylabel('Gamma')

%Wing Lift Coefficient vs alpha
CL=pi*AR*A(1);  
figure(1)
subplot(1,2,2)
fplot(CL,[-7*pi/180,7*pi/180]),title('C_L vs. \alpha'),xlabel('\alpha'),ylabel('C_L')
hold on
%Wing tip airfoil lift coefficient vs alpha
Cltip=4*b*subs(alstuff,o,pi)/subs(c,o,pi); 
fplot(Cltip,[-7*pi/180,7*pi/180]),title('C_l/L vs. \alpha'),xlabel('\alpha'),ylabel('C_l/L')
%Wing root airfoil lift coefficient vs alpha
Clroot=4*b*subs(alstuff,o,pi/2)/subs(c,o,pi/2); 
fplot(y,[-7*pi/180,7*pi/180]), legend('C_L','C_ltips','C_lroot')


%% Prandtl Glauert Correction
CL_0=subs(CL,alpha,alfa);
M=linspace(.2,.8,100);
CL_comp=zeros(100,1);
for i=1:100
    CL_comp(i)=CL_0/sqrt(1-M(i)^2);
end

figure(2)
plot(M,CL_comp),title('C_L vs Mach number'),xlabel('Mach number'),ylabel('C_L')


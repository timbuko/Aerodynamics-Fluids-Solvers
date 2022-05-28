%Solve a wave Reflection then the reflection and refraction with contact
%surface

syms p6_pp p6_p5 Wrr Wref u6

%input
W=
a1=

%Solve region 5
Ms=W/a1;




%Define
p2_p1a=3.507999999985390120001957675104;
gamma=1.4;
u2=429.2152454863963714236; %same as u3 for shock relecting off contact surface
a2=339.4;%a3 for shock relecting off contact surface
a1b

eq1=p5_p2*p2_p1a-p5_p1;
eq2=-p5_p2+1+2*gamma/(gamma+1)*((Wr+u2)^2/a2^2-1);
eq3=-p5_p1+1+2*gamma/(gamma+1)*((Wref)^2/a1b^2-1);
eq4=-Wref+a1b*sqrt((gamma+1)/(2*gamma)*(p5_p1-1)+1);
eq5=-u5+a1/gamma*(p5_p1-1)*sqrt(2*gamma/(gamma+1)/(p5_p1+(gamma-1)/(gamma+1)));

%Solve system of eq (set all eq=0)

S=solve(eq1==0,eq2==0,eq3==0,eq4==0,eq5==0)
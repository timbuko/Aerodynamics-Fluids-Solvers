%Find air properties inside RamJet

ramjet.mach(1)=4;
ramjet.pressure(1)=30; %kPa
ramjet.temperature(1)=-45; %C
ramjet.temperature(1)=ramjet.temperature(1)+273.15; %Convert to Kelvin
ramjet.total_pressure(1)=ramjet.pressure(1)*(1+(gamma-1)*ramjet.mach(1)^2/2)^(gamma/(gamma-1))
gamma=1.4;

theta(1)=15; %deg
theta(2)=15;

theta=theta*pi/180; %Convert to rad

%First shock
beta=findbet(ramjet.mach(1),theta(1),gamma);
mn_inf=ramjet.mach(1)*sin(beta);
p1_pinf=presratoblique(mn_inf,gamma);
ramjet.pressure(2)=ramjet.pressure(1)*p1_pinf;
t1_tinf=tempratoblique(mn_inf,gamma);
ramjet.temperature(2)=ramjet.temperature(1)*t1_tinf;

mn1=oblique(mn_inf,gamma);
ramjet.mach(2)=mn1/sin(beta-theta(1));
ramjet.total_pressure(2)=ramjet.pressure(2)*(1+(gamma-1)*ramjet.mach(2)^2/2)^(gamma/(gamma-1));

%Second shock
beta2=findbet(ramjet.mach(2),theta(2),gamma);
mn2=ramjet.mach(2)*sin(beta2);
p3_p2=presratoblique(mn2,gamma);
ramjet.pressure(3)=ramjet.pressure(2)*p3_p2;
t3_t2=tempratoblique(mn2,gamma);
ramjet.temperature(3)=ramjet.temperature(2)*t3_t2;

mn3=oblique(mn2,gamma);
ramjet.mach(3)=mn3/sin(beta2-theta(2));
ramjet.total_pressure(3)=ramjet.pressure(3)*(1+(gamma-1)*ramjet.mach(3)^2/2)^(gamma/(gamma-1));

%Third normal shock
p4_p3=presratoblique(ramjet.mach(3),gamma);
ramjet.pressure(4)=ramjet.pressure(3)*p4_p3;
t4_t3=tempratoblique(ramjet.mach(3),gamma);
ramjet.temperature(4)=ramjet.temperature(3)*t4_t3

ramjet.mach(4)=oblique(ramjet.mach(3),gamma)
ramjet.total_pressure(4)=ramjet.pressure(4)*(1+(gamma-1)*ramjet.mach(4)^2/2)^(gamma/(gamma-1));

%If only a normal shock
pp_pinf=presratoblique(M_inf,gamma);
p_p=ramjet.pressure(1)*pp_pinf
tp_tinf=tempratoblique(M_inf,gamma);
t_p=ramjet.temperature(1)*tp_tinf

Mp=oblique(ramjet.mach(1),gamma)
pt_p=p_p*(1+(gamma-1)*Mp^2/2)^(gamma/(gamma-1))

ramjet.temperature=ramjet.temperature-273.15
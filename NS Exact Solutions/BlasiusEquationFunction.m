%This script creates the Blasius function   f(eta)

[eta,f]=ode45(@blasiusode,[0,10],[0,0,0.33206]) 
        %Third initial condition was found by others and gives f(inf)=1
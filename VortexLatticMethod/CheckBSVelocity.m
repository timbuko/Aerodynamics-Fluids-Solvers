%%%I CHANGED BSvelocity so it no longer inputs gamma%%%%%%

%Part 1  One side of Square Filament%%%%%%%%%%

a=2*sqrt(2)/pi;  %Side length
gamma=1; 
A=[0,0,0]; 
B=[0,a,0];
C=[a/2,a/2,0]; %point in middle of square

v=BSvelocity(A,B,C,gamma)
vmag=sqrt(dot(v,v))

%Part 2 Semi infinite filament 

A=[0,0,0];
B=[10000000,0,0];
C=[0,1,0];
gamma=1;

v=BSvelocity(A,B,C,gamma)
vmag=sqrt(dot(v,v))

% %Part 3 Semi infinite filament in opposite direction
A=[10000000,0,0];
B=[0,0,0];
C=[0,1,0];
gamma=1;

v=BSvelocity(A,B,C,gamma)
vmag=sqrt(dot(v,v))
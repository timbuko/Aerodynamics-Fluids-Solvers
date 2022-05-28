function V = horseshoe(A,B,K)
%This function takes in vertices of horseshoe filament A,B
% point of influence K, and circulation G

v1=BSvelocity([99999999,A(2),A(3)],A,K);
v2=BSvelocity(A,B,K);
v3=BSvelocity(B,[9999999,B(2),B(3)],K);

V=v1+v2+v3;

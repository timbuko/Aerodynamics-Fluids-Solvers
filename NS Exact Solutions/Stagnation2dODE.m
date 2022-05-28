function dAdeta = Stagnation2dODE(eta,A)
    %2D Steady Stagnation Flow Eq A4+A1*A3-A2^2+1=0
dAdeta=zeros(3,1);
dAdeta(1)=A(2);
dAdeta(2)=A(3);
dAdeta(3)=-A(1)*A(3)+A(2)^2-1;


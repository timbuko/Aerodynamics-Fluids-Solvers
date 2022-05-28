function dAdeta = blasiusode(eta,A)
    %Blasius equation 2*A4+A1*A3=0
dAdeta=zeros(3,1);
dAdeta(1)=A(2);
dAdeta(2)=A(3);
dAdeta(3)=-A(1)*A(3)/2;


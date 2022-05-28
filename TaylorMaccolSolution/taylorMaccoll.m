function dur_dtheta = taylorMaccoll(theta,ur,gamma)
    %Taylor Maccoll Eqn
    %ur=[ur ur']=[ur utheta]
    
    dur_dtheta=zeros(2,1);
    
    dur_dtheta(1)=ur(2);
    dur_dtheta(2)=(ur(2)^2*ur(1)-(gamma-1)/2*(1-ur(1)^2-ur(2)^2)*(2*ur(1)+ur(2)*cot(theta)))...
/((gamma-1)/2*(1-ur(1)^2-ur(2)^2)-ur(2)^2);
    
    
end

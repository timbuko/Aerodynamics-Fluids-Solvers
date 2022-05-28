function [mafter] = findm2pmgen(nu,gamma);
%%%%%%%%
%  nu must be in radians!!!
% finds the mach number associated with a given prandtl meyer angle
% Recursive Method, and dependent on pmang function
%%%%%%%
    
M_start = 1;    %recursion starts at the min possible mach
step = 0.001;   %mach increment step
counter = 0;    %a counter is a counter
while counter > -1  %infinite while loop until we find what we want
    M_test = M_start + counter*step;    %calculate the current test mach number
    nu_test = pmang(M_test,gamma);      %calculate the output nu using the pmang function
    if nu_test > nu                     %since we're starting with the mach number of smallest possible nu, eventually our test nu will exceed the nu we're looking for
        M_prev = (M_test - step);       %calculate the mach and nu for one step before we exceeded the nu we want
        nu_prev = pmang(M_prev,gamma);
        mafter = M_prev + (((nu - nu_prev)/(nu_test - nu_prev))*step);  %linear interoplation between the two mach numbers that has our nu in between
        counter = -1;   %break the loop since we found our mach
    else
        counter = counter + 1;  %increase the counter if our test nu didn't exceed our wanted nu yet
    end
end
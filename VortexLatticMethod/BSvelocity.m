function v =BSvelocity(A,B,C)
%this function inputs 2 points on a filament, a point, and gamma to find velocity
%through that point

r1=C-A;
r2=C-B;
r0=B-A;
r1mag=sqrt(dot(r1,r1));
r2mag=sqrt(dot(r2,r2));
r12=cross(r1,r2);

v=1/(4*pi)*r12/dot(r12,r12)*dot(r0,((r1/r1mag)-(r2/r2mag)));
clc;
clear all;
close all;

% define initial tangent
p.tx0 = 1;
p.ty0 = 0;
p.tz0 = 0;

% define initial normal
p.nx0 = 0;
p.ny0 = 1; 
p.nz0 = 0; 

% define initial binormal
p.bx0 = 0;
p.by0 = 0; 
p.bz0 = 1; 

% define initial point
p.x0  = 0;
p.y0  = 0;
p.z0  = 0;

% define clothoid parameters
p.dk  = 0.1;
p.k0  = 0;
p.dt  = 0.1;
p.t0  = 0;

% define lenths/n.points
L = 30;
h = 0.1;
n = floor(L/h)+1;
s = 0:h:L;


% compute canonical clothoid
[pos,~,~,~] = clothoid3d(s,p);

% compute general clothoid commutator-free
solMG4 = MG4(L,n,p);

% compute general clothoid commutator-free
solCF4GL = CF4GL(L,n,p);





% compare plots
plot3(pos(1,:),pos(2,:),pos(3,:)); hold on;
plot3(solCF4GL(10,:),solCF4GL(11,:),solCF4GL(12,:)); 
plot3(solMG4(10,:),solMG4(11,:),solMG4(12,:));
legend("Canonical closed form", "CF4GL", "MG4");

axis equal

% measure error on endpoint
Error_final_point_MG4   = norm(pos(:,end)-solMG4(10:12,n+1),2)
Error_final_point_CF4GL = norm(pos(:,end)-solCF4GL(10:12,n+1),2)



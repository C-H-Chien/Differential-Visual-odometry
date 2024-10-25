%% Verification of closed-form solution of 3D clothoid model by synthetic data
clc; clear all; close all;

% -- params to create synthetic curve --
s_c = 0:0.05:3;
k0 = 0;
tau0 = 0;
k_prime = 5;
tau_prime = 3.5;
init_T = [1;0;0];
init_N = [0;1;0];
init_B = [0;0;1];
init_C = [0;0;0];

tangent = init_T;
normal = init_N;
binormal = init_B;
curveStartPt = init_C;
synthCurve = [];
synthCurve = [synthCurve, init_C];
for i = 1:size(s_c,2)
    k = k0 + k_prime*s_c(1,i);
    tau = tau0 + tau_prime*s_c(1,i);
    
    [synthCurvePt, propogated_FF] = self_generateHelixFromModel(k, tau, tangent, normal, binormal, s_c(1,2), curveStartPt);
    synthCurve = [synthCurve, synthCurvePt];
    
    curveStartPt = synthCurvePt;
    tangent = propogated_FF.T;
    normal = propogated_FF.N;
    binormal = propogated_FF.B;
end


[GeneralClothoid, ~, sol_check] = self_generateGeneralClothoid(k0, tau0, k_prime, tau_prime, init_T, init_N, init_B, s_c, init_C);

[CanonicalClothoid, ~] = self_generateCanonicalClothoid(k_prime, tau_prime, init_T, init_N, init_B, s_c, init_C);

[GeneralCommutatorClothoid, ~, FregoCurve] = self_generateGeneralCommutatorClothoid(k0, tau0, k_prime, tau_prime, init_T, init_N, init_B, s_c, init_C);

% % from Dr. Frego
% % define lenths/n.points
% % replicate the parameters from my side
% % define initial tangent
% p.tx0 = 1;
% p.ty0 = 0;
% p.tz0 = 0;
% 
% % define initial normal
% p.nx0 = 0;
% p.ny0 = 1; 
% p.nz0 = 0; 
% 
% % define initial binormal
% p.bx0 = 0;
% p.by0 = 0; 
% p.bz0 = 1; 
% 
% % define initial point
% p.x0  = 0;
% p.y0  = 0;
% p.z0  = 0;
% 
% % define clothoid parameters
% p.dk  = 0.1;
% p.k0  = 0;
% p.dt  = 0.1;
% p.t0  = 0;
% 
% 
% L = s_c(1, end);
% h = 0.5;
% n = floor(L/h)+1;
% s = 0:h:L;
% solCF4GL = CF4GL(L,n,p);

% criteria = [];
% for i = 1:size(s,2)
%     val = (k0^2+tau0^2)*s(1,i) + (k0*k_prime + tau0*tau_prime)*s(1,i)^2 + (1/3)*(k_prime^2+tau_prime^2)*s(1,i)^3;
%     criteria = [criteria, val];
% end

figure;
%plot3(CanonicalClothoid(1,:), CanonicalClothoid(2,:), CanonicalClothoid(3,:), 'ro-', 'DisplayName', 'Canonical Clothoid');
%hold on;
plot3(FregoCurve(1,:), FregoCurve(2,:), FregoCurve(3,:), 'ro-', 'DisplayName', 'General Clothoid w Commutator');
hold on;
plot3(GeneralClothoid(1,:), GeneralClothoid(2,:), GeneralClothoid(3,:), 'go-', 'DisplayName', 'General Clothoid');
hold on;
%plot3(synthCurve(1,:), synthCurve(2,:), synthCurve(3,:), 'bo-', 'DisplayName', 'approximate Clothoid');
xlabel('x');
ylabel('y');
zlabel('z');
legend;
axis equal;
set(gcf,'color','w');
clc; clear all; close all;

% -- params to create synthetic curve --
s_c = 0:1:30;
k0 = 0.0001;
tau0 = 0.0001;
k_prime = 0.01;
tau_prime = 0.01;
init_T = [1;0;0];
init_N = [0;1;0];
init_B = [0;0;1];
init_C = [0;0;0];

% General 3d Clothoid
[GeneralClothoid, ~] = self_generateGeneralClothoid(0, 0, k_prime, tau_prime, init_T, init_N, init_B, s_c, init_C);

% Canonical 3d Clothoid
[CanonicalClothoid, ~] = self_generateCanonicalClothoid(k_prime, tau_prime, init_T, init_N, init_B, s_c, init_C);

figure;
plot3(CanonicalClothoid(1,:), CanonicalClothoid(2,:), CanonicalClothoid(3,:), 'ro-', 'DisplayName', 'Canonical Clothoid');
hold on;
plot3(GeneralClothoid(1,:), GeneralClothoid(2,:), GeneralClothoid(3,:), 'go-', 'DisplayName', 'General Clothoid');
xlabel('x');
ylabel('y');
zlabel('z');
legend;
axis equal;
set(gcf,'color','w');

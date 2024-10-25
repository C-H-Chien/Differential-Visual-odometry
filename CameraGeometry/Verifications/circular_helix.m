
% -- plot the helix using samples and true velocity, acceleration,
% curvature, and torsion --
clc; clear all; close all;

% -- define parameters --
a = 10;
b = 7;
start_pt = 20;
end_pt = 60;

% -- circular helix in terms of arc length s --
s = 0:100;
sqrt_ab = sqrt(a^2+b^2);
t = s / sqrt_ab;
init_r = zeros(size(s,2), 3);
r = zeros(size(s,2), 3);
vel = zeros(size(s,2), 3);
sped = zeros(size(s,2), 1);
acc = zeros(size(s,2), 3);
a0 = zeros(size(s,2), 3);
mag_acc = zeros(size(s,2), 1);
curv = zeros(size(s,2), 1);
tors = zeros(size(s,2), 1);

% -- positions --
init_r(:,1) = a*cos(s(1,:) / sqrt_ab);
init_r(:,2) = a*sin(s(1,:) / sqrt_ab);
init_r(:,3) = b*s(1,:)/sqrt_ab;

% -- add perturbations to the positions --
bnd_upper = 0.05;
bnd_lower = -0.05;
rng('default');

%noise_sigma = 0.05;
noise_sigma = 0;
r(:,1) = init_r(:,1) + noise_sigma*randn(size(init_r(:,1), 1),1) + 0;
r(:,2) = init_r(:,2) + noise_sigma*randn(size(init_r(:,2), 1),1) + 0;
r(:,3) = init_r(:,3) + noise_sigma*randn(size(init_r(:,3), 1),1) + 0;

% -- velocity --
vel(:,1) = -a*sin(t(1,:));
vel(:,2) = a*cos(t(1,:));
vel(:,3) = b;
for i = 1:size(vel, 1)
    sped(i, 1) = norm(vel(i,:));    % -- CORRECT! --
end

% -- acceleration --
% -- acceleration is the change of speed over time --
for i = 1:size(vel, 1)
    if i == 1
        a0(i,1) = 0;
    else
        a0(i,1) = (sped(i,1) - sped(i-1,1))/(t(1,i)-t(1,i-1));
    end
end

% -- acceleration in each direction --
acc(:,1) = -a*cos(t(1,:));
acc(:,2) = -a*sin(t(1,:));
acc(:,3) = 0;
for i = 1:size(acc, 1)
    mag_acc(i, 1) = norm(acc(i,:));
end

% -- derivative of acceleration --
acc_prime(:,1) = a*sin(t(1,:));
acc_prime(:,2) = -a*cos(t(1,:));
acc_prime(:,3) = 0;
for i = 1:size(acc, 1)
    mag_acc_prime(i, 1) = norm(acc_prime(i,:));
end

% -- curvature --
curv(:,1) = mag_acc(:,1) ./ (a^2 + b^2);

% -- torsion --
tors(:,1) = b / (a^2 + b^2);


% ====================================================================================
% -- construt the model for C(s) --
% -- redefining notations --
k_0 = curv;
tau_0 = tors;

omega_0 = zeros(size(curv));
for i = 1:size(curv, 1)
    omega_0(i,1) = sqrt(k_0(i,1)^2 + tau_0(i,1)^2);
end

estimated_s = zeros(size(s,2), 1);
for i = 1:size(s,2)
    estimated_s(i,1) = s(1,1) + sped(i,1)*t(1,i) + 0.5*acc(i, 1)*t(1,i)^2;
end

% -- first-order approximation --
% -- initial T, N, B when s = 0 --
T_0 = [0, a/sqrt_ab, b/sqrt_ab];
N_0 = [-1, 0, 0];
B_0 = cross(T_0, N_0);

C = zeros(3, size(s, 2));
for i = 1:size(s,2)
    ds = s(1,1) + sped(i,1)*t(1,i) + a0(i,1)*0.5*(t(1,i)^2);% + cg_mag_trip*(s(1,i)^3)/6;
    %ds = s(1,i);
    C_1 = [(tau_0(i,1)/omega_0(i,1))^2*ds + (k_0(i,1)^2/omega_0(i,1)^3)*sin(omega_0(i,1)*ds); ...
           (k_0(i,1)/(omega_0(i,1)^2))*(1-cos(omega_0(i,1)*ds)); ...
           (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^2)*ds - (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^3)*sin(omega_0(i,1)*ds)];
    
    C(:, i) = r(1,:)' + C_1(1,1)*T_0' + C_1(2,1)*N_0' + C_1(3,1)*B_0';
end
figure;
plot3(C(1,:), C(2,:), C(3,:), 'ro-');
xlabel('x');
ylabel('y');
zlabel('z');
set(gcf,'color','w');
% ====================================================================================
% -- plot the circular helix trajectory --
figure;
plot3(r(:,1), r(:,2), r(:,3), 'bo-', 'DisplayName', 's \in [0, 100]');
legend;
xlabel('x');
ylabel('y');
zlabel('z');
set(gcf,'color','w');

% ======================================================================================


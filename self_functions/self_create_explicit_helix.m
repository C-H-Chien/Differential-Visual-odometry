function [r, dynamicParams, geometryParams] = self_create_explicit_helix(s, a, b, add_noise, shift2origin)


init_r = zeros(size(s,2), 3);
r = zeros(size(s,2), 3);
vel = zeros(size(s,2), 3);
sped = zeros(size(s,2), 1);
acc = zeros(size(s,2), 3);
a0 = zeros(size(s,2), 3);
mag_acc = zeros(size(s,2), 1);
curv = zeros(size(s,2), 1);
tors = zeros(size(s,2), 1);

sqrt_ab = sqrt(a^2+b^2);
t = s / sqrt_ab;

if shift2origin
    initial_pt = [a*cos(s(1,1) / sqrt_ab), a*sin(s(1,1) / sqrt_ab), b*s(1,1)/sqrt_ab];
else
    initial_pt = [0, 0, 0];
end

% -- positions --
init_r(:,1) = a*cos(s(1,:) / sqrt_ab) - initial_pt(1);
init_r(:,2) = a*sin(s(1,:) / sqrt_ab) - initial_pt(2);
init_r(:,3) = b*s(1,:)/sqrt_ab - initial_pt(3);

rng('default');
% -- add perturbations to the positions --
if add_noise
    noise_sigma = 0.05;
else
    noise_sigma = 0;
end
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
    % -- MAY NOT BE USED --
    mag_acc_prime(i, 1) = norm(acc_prime(i,:));
end

% -- curvature --
curv(:,1) = mag_acc(:,1) ./ (a^2 + b^2);

% -- torsion --
tors(:,1) = b / (a^2 + b^2);

% -- compute Frenet frame at the first point --


% -- return the geometry parameters --
dynamicParams.vel_vec = vel;
dynamicParams.acc = a0;
dynamicParams.acc_vec = acc;
geometryParams.curvature = curv;
geometryParams.torsion = tors;

end
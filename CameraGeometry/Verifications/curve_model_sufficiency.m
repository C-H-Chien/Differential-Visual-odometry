
% -- Verify that the camera center model C(t) is a sufficient model using KITTI dataset --
clc; clear all; close all;

% -- define parameters --
datasetName = 'KITTI';
start_fr = 1;
end_fr = 20;
num_fr = end_fr-start_fr+1;

% -- parameters for segements plottings --
show_fr_num = end_fr - start_fr;
remain = num_fr - show_fr_num;
buffer_size = round((remain)/2);
start_t = start_fr+buffer_size;
end_t = end_fr-buffer_size;

if strcmp(datasetName, 'KITTI')
    % -- choose sequence --
    sequence = '00';
    
    % -- image data directories --
    imgFolder = '/home/chchien/datasets/KITTI/sequences-gray/';
    imgDir = strcat(imgFolder, sequence);

    % -- read real poses from KITTI ground truths --
    poseFolder = '/home/chchien/BrownU/research/Differential-Visual-Odometry/';
    dataset = 'KITTI-gt/poses/';
    postfix = '.txt';
    
    pose_FileName = strcat(poseFolder, dataset, sequence, postfix);
    [all_R, all_T] = read_groundtruths_kitti(pose_FileName);
    
    % -- read time of each pose --
    timeFile = '/times.txt';
    times_FileName = fullfile(imgDir, timeFile);
    timesFileRd = fopen(times_FileName, 'r');
    ldata = textscan(timesFileRd, '%s', 'CollectOutput', true);
    line = string(ldata{1});
    time_per_fr = str2double(line);
    
    % -- segment and initializations --
    r = all_T(:,start_fr:end_fr)';
    cg_vel = zeros(size(r,1), 3);       % -- velocities --
    cg_sped = zeros(size(r,1), 1);      % -- speed (magnitude of velocity) --
    cg_mag_acc = zeros(size(r,1), 1);   % -- magnitude of accelerations --
    cg_curv = zeros(size(r,1), 1);      % -- curvatures --
    cg_tors = zeros(size(r,1), 1);      % -- torsions --
    cg_mag_trip = zeros(size(r,1), 1);  % -- magnitude of triple primes --
    cg_curv_prime = zeros(size(r,1), 1);% -- derivative of curvature --

    sigma = 2.2;
    data_range = -10:1:10;
    ker_size = data_range';
    [G1, ~] = gaussian_derivatives(0, ker_size, sigma, 0);
    Gfp = gradient(G1);
    Gfp_2 = gradient(Gfp);
    Gfp_3 = gradient(Gfp_2);

    r_conv(:,1) = conv(r(:,1)', G1, 'same')';
    r_conv(:,2) = conv(r(:,2)', G1, 'same')';
    r_conv(:,3) = conv(r(:,3)', G1, 'same')';

    vel_conv(:,1) = conv(r_conv(:,1)', Gfp./time_per_fr(2,1), 'same')';
    vel_conv(:,2) = conv(r_conv(:,2)', Gfp./time_per_fr(2,1), 'same')';
    vel_conv(:,3) = conv(r_conv(:,3)', Gfp./time_per_fr(2,1), 'same')';

    fd_x = vel_conv(:,1);
    fd_y = vel_conv(:,2);
    fd_z = vel_conv(:,3);
    acc_conv(:,1) = conv(fd_x, Gfp./time_per_fr(2,1), 'same')';
    acc_conv(:,2) = conv(fd_y, Gfp./time_per_fr(2,1), 'same')';
    acc_conv(:,3) = conv(fd_z, Gfp./time_per_fr(2,1), 'same')';

    dd_x = acc_conv(:,1);
    dd_y = acc_conv(:,2);
    dd_z = acc_conv(:,3);
    trip_conv(:,1) = conv(dd_x, Gfp./time_per_fr(2,1), 'same')';
    trip_conv(:,2) = conv(dd_y, Gfp./time_per_fr(2,1), 'same')';
    trip_conv(:,3) = conv(dd_z, Gfp./time_per_fr(2,1), 'same')';

    for i = 1:size(r, 1)
        cg_sped(i, 1) = norm(vel_conv(i,:));
        cg_mag_acc(i, 1) = norm(acc_conv(i,:));
        cg_mag_trip(i,1) = norm(trip_conv(i,:));
    end

    % -- initialization --
    unit_tangent = zeros(size(vel_conv));
    unit_normal = zeros(size(vel_conv));
    unit_binormal = zeros(size(vel_conv));
    C_double_prime = zeros(size(vel_conv));
    C_triple_prime = zeros(size(vel_conv));
    curv_kappa = zeros(size(unit_tangent, 1), 1);
    tors_tau = zeros(size(unit_tangent, 1), 1);
    acc_prime = zeros(size(unit_tangent, 1), 1);
    acc_mag = zeros(size(unit_tangent, 1), 1);
    
    % -- a) compute a(t) and a'(t) --
    acc_mag(:,1) = conv(cg_sped, Gfp./time_per_fr(2,1), 'same')';
    acc_prime(:,1) = conv(acc_mag, Gfp./time_per_fr(2,1), 'same')';
    
    % -- b) compute unit tangent vector T(t) --
    for vi = 1:size(vel_conv, 1)
        unit_tangent(vi,:) = vel_conv(vi,:) ./ cg_sped(vi, 1);
    end
    
    % -- c) compute unit normal vector N(t) and C''(t)--
    unit_tangent_prime(:,1) = conv(unit_tangent(:,1), Gfp./time_per_fr(2,1), 'same')';
    unit_tangent_prime(:,2) = conv(unit_tangent(:,2), Gfp./time_per_fr(2,1), 'same')';
    unit_tangent_prime(:,3) = conv(unit_tangent(:,3), Gfp./time_per_fr(2,1), 'same')';
    for vi = 1:size(unit_tangent_prime, 1)
        unit_normal(vi,:) = unit_tangent_prime(vi,:) ./ norm(unit_tangent_prime(vi,:));
        %test_T_dot_N(vi,:) = dot(unit_tangent(vi,:), unit_normal(vi,:));
        
        % -- C''(t) --
        C_double_prime(vi,:) = acc_mag(vi,1) * unit_tangent(vi,:) + cg_sped(vi,1) * unit_tangent_prime(vi,:);
    end
    
    % -- d) compute unit binormal vector B(t) --
    for bi = 1:size(unit_normal, 1)
        unit_binormal(bi, :) = cross(unit_tangent(bi,:), unit_normal(bi,:));
        unit_binormal(bi, :) = unit_binormal(bi, :) ./ norm(unit_binormal(bi, :));
    end
    
    % -- e) curvature --
    for i = 1:size(r, 1)
        cross_prod_va = cross(vel_conv(i,:), acc_conv(i,:));
        % numer = sign(cross_prod_va(1,2)) * norm(cross_prod_va);
        numer = norm(cross_prod_va);
        denom = norm(vel_conv(i,:))^3;
        cg_curv(i,1) = numer / denom;
    end
    
    % -- f) compute T''(t) and C'''(t) --
    test_C_triple_prime(:,1) = conv(C_double_prime(:,1), Gfp./time_per_fr(2,1), 'same')';
    test_C_triple_prime(:,2) = conv(C_double_prime(:,2), Gfp./time_per_fr(2,1), 'same')';
    test_C_triple_prime(:,3) = conv(C_double_prime(:,3), Gfp./time_per_fr(2,1), 'same')';
    
    unit_tangent_double_prime(:,1) = conv(unit_tangent_prime(:,1), Gfp./time_per_fr(2,1), 'same')';
    unit_tangent_double_prime(:,2) = conv(unit_tangent_prime(:,2), Gfp./time_per_fr(2,1), 'same')';
    unit_tangent_double_prime(:,3) = conv(unit_tangent_prime(:,3), Gfp./time_per_fr(2,1), 'same')';
    for vi = 1:size(unit_tangent_prime, 1)
        C_triple_prime(vi,:) = acc_prime(vi,1) * unit_tangent(vi,:) + 2*acc_mag(vi,1) * unit_tangent_prime(vi,:) ...
            + cg_sped(vi,1) * unit_tangent_double_prime(vi,:);
    end
    
    % -- g) torsion --
    for i = 1:size(r,1)
        tors_tau(i,1) = dot(C_triple_prime(i,:), unit_binormal(i,:)) / (cg_sped(i, 1)^2);
        %tors_tau(i,1) = dot(test_C_triple_prime(i,:), unit_binormal(i,:)) / (curv_kappa(i,1) * (cg_sped(i, 1)^2));
    end


    % -- plot KITTI trajectory --    
    figure;
    plot3(all_T(1,start_t:end_t), all_T(2,start_t:end_t), all_T(3,start_t:end_t), 'bo-');
    hold on;
    plot3(all_T(1,start_t), all_T(2,start_t), all_T(3,start_t), 'go');
    text(all_T(1,start_t), all_T(2,start_t), all_T(3,start_t), 'start', 'FontSize', 10, 'Color', 'g');
    hold on;
    plot3(all_T(1,end_t), all_T(2,end_t), all_T(3,end_t), 'ro');
    text(all_T(1,end_t), all_T(2,end_t), all_T(3,end_t), 'end', 'FontSize', 10, 'Color', 'r');
    xlabel("x (m)");
    ylabel("y (m)");
    zlabel("z (m)");
    axis equal;
    set(gcf,'color','w');
    hold on;
end

% -- construt the model for C(s) --
% -- redefining notations --
k_0 = cg_curv;
tau_0 = tors_tau;
T_0 = unit_tangent;
N_0 = unit_normal;
B_0 = unit_binormal;

omega_0 = zeros(size(cg_curv));
for i = 1:size(cg_curv, 1)
    omega_0(i,1) = sqrt(k_0(i,1)^2 + tau_0(i,1)^2);
end

%s = zeros(size(r, 1), 1);
% -- first-order approximation --
s = cg_sped*time_per_fr(2,1) + cg_mag_acc*0.5*(time_per_fr(2,1)^2) + cg_mag_trip*(time_per_fr(2,1)^3)/6;
C_0 = [all_T(1,start_t); all_T(2,start_t); all_T(3,start_t)];
C = zeros(3, size(start_t:end_t, 2));
C(:, 1) = C_0;
for i = buffer_size+1:num_fr-buffer_size-1
    C_1 = [(tau_0(i,1)/omega_0(i,1))^2*s(i,1) + (k_0(i,1)^2/omega_0(i,1)^3)*sin(omega_0(i,1)*s(i,1)); ...
           (k_0(i,1)/(omega_0(i,1)^2))*(1-cos(omega_0(i,1)*s(i,1))); ...
           (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^2)*s(i,1) - (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^3)*sin(omega_0(i,1)*s(i,1))];
    
    C(:, i) = C(:,i-1) + C_1(1,1)*T_0(i,:)' + C_1(2,1)*N_0(i,:)' + C_1(3,1)*B_0(i,:)';
end
plot3(C(1,buffer_size:num_fr-buffer_size-1), C(2,buffer_size:num_fr-buffer_size-1), C(3,buffer_size:num_fr-buffer_size-1), 'ro-');

% ----------- PLOTTINGS -------------
% -- plot Gaussian convolved speed --
figure;
subplot(1,6,1);
plot(start_t:end_t, cg_sped(buffer_size:num_fr-buffer_size-1, 1), 'b-', 'DisplayName', 'computed speed');
xlim([start_t, end_t]);
title({'speed'});
xlabel("points");
ylabel("speed (m/s)");

% -- plot Gaussian convolved acceleration --
subplot(1,6,2);
plot(start_t:end_t, cg_mag_acc(buffer_size:num_fr-buffer_size-1, 1), 'b-');
xlim([start_t, end_t]);
%ylim([0, 3.5]);
title('acceleration');
xlabel("points");
ylabel("acceleration (m/s^2)");

% -- plot Gaussian convolved curvature --
subplot(1,6,3);
plot(start_t:end_t, cg_curv(buffer_size:num_fr-buffer_size-1, 1), 'b-');
xlim([start_t, end_t]);
%ylim([-0.15, 0.15]);
title('curvature');
xlabel("points");
ylabel("curvature (1/m)");

% -- plot Gaussian convolved torsion --
subplot(1,6,4);
% plot(start_t:end_t, cg_tors(buffer_size:num_fr-buffer_size-1, 1), 'b-');
% hold on;
plot(start_t:end_t, tors_tau(buffer_size:num_fr-buffer_size-1, 1), 'b-');
xlim([start_t, end_t]);
ylim([-1, 1]);
title('torsion');
xlabel("points");
ylabel("torsion");

% -- plot thrid derivative of curve --
subplot(1,6,5);
plot(start_t:end_t, cg_mag_trip(buffer_size:num_fr-buffer_size-1, 1), 'b-');
xlim([start_t, end_t]);
%ylim([0, 10]);
title({'derivative of'},{'acceleration'});
xlabel("points");
ylabel("derivative of acceleration (m/s^3)");

% -- plot derivative of curvature --
subplot(1,6,6);
plot(start_t:end_t, cg_curv_prime(buffer_size:num_fr-buffer_size-1, 1), 'b-');
xlim([start_t, end_t]);
%ylim([-0.5, 0.5]);
title({'derivative of'},{'curvature'});
xlabel("points");
ylabel('derivative of curvature (1/m*s)');

set(gcf,'color','w');

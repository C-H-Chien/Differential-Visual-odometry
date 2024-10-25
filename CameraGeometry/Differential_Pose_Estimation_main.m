
% -- do sliding window over a trajectory and plot the optimized poses 
% accompanied by KITTI ground truth trajectory --
clc; clear all; close all;

% -- image data directories --
imgFolder = '/home/chchien/datasets/KITTI/sequences-gray/00';
category_l = '/image_0/';
category_r = '/image_1/';

% -- intrinsic matrix --
K = [718.856, 0, 607.193; 0, 718.86, 185.216; 0, 0, 1];

% -- define parameters --
start_fr = 725;
plot_start_fr = start_fr;
covisible_frames = 5;
plot_all = 0;
plot_tracks = 20;
plot_trajectory = 1;
sliding_window_length = 10;
use_monocular = 1;
use_ORB = 0;
e_3 = [0,0,1];
debug = 0;

% -- read time of each pose --
timeFile = '/times.txt';
times_FileName = fullfile(imgFolder, timeFile);
timesFileRd = fopen(times_FileName, 'r');
ldata = textscan(timesFileRd, '%s', 'CollectOutput', true);
line = string(ldata{1});
time_per_fr = str2double(line);

% -- read real poses from KITTI ground truths --
poseFolder = '/home/chchien/BrownU/research/Differential-Visual-Odometry/';
dataset = 'KITTI-gt/poses/';
sequence = '00.txt';
pose_FileName = fullfile(poseFolder, dataset, sequence);
[all_R, all_T] = read_groundtruths_kitti(pose_FileName);

% -- initialize all poses --
trajectory_R = all_R;
trajectory_T = all_T;
pg_trajectory_T = all_T;

% -- set camera coordinate of R and T --
R_cc = zeros(3,3,covisible_frames);
T_cc = zeros(3,1,covisible_frames);

% -- slide the window --
for w = 1:sliding_window_length+1

    % -- initialization --
    cov_indx = 1;
    if w > 1
        start_fr = start_fr+1;
    end
    fprintf("\n=== Sliding Window Starting on Frame %d ===\n", start_fr+1);
    
    % -- perform inside a sliding window --
    for i = start_fr-1:start_fr+covisible_frames-1
        imgName = num2str(i,'%06.f.png');
        fr1_FileName = fullfile(imgFolder, category_l, imgName);
        frR_FileName = fullfile(imgFolder, category_r, imgName);

        if i == start_fr - 1    % ???
            % -- read and initialize an image --
            img_prev = imInit_in_SURF(fr1_FileName);

            % -- detect and extract features --
            if use_ORB
                keyPointsObj_prev = detectORBFeatures(img_prev);
                [f_prev_Features, vpts_prev] = extractFeatures(img_prev, keyPointsObj_prev);
                f_prev = f_prev_Features.Features;
            else
                keyPointsObj_prev = detectSURFFeatures(img_prev);
                [f_prev, vpts_prev] = extractFeatures(img_prev, keyPointsObj_prev);
            end

            continue;
        elseif i == start_fr
            % -- read and initialize an image --
            ref_img = imInit_in_SURF(fr1_FileName);

            % -- detect and extract features --
            if use_ORB
                keyPointsObj = detectORBFeatures(ref_img);
                [f1_Features, vpts1] = extractFeatures(ref_img, keyPointsObj);
                f1 = f1_Features.Features;
            else
                keyPointsObj = detectSURFFeatures(ref_img);
                [f1, vpts1] = extractFeatures(ref_img, keyPointsObj);
            end
            f_size = size(keyPointsObj.Location,1);

            % -- initialize empty covisible track arrays based on feature size --
            covis_tracks = zeros(f_size, 5, covisible_frames);            % -- observation x, observation y, depth, gamma x, gamma y --
            covis_observe_tracks_u = zeros(f_size, covisible_frames);     % -- feature tracks, x --
            covis_observe_tracks_v = zeros(f_size, covisible_frames);     % -- feature tracks, y --
            covis_gamma_tracks_u = zeros(f_size, covisible_frames);       % -- gamma(t) tracks, x --
            covis_gamma_tracks_v = zeros(f_size, covisible_frames);       % -- gamma(t) tracks, y --

            % -- store the locations of features --
            covis_tracks(:,1, cov_indx) =  keyPointsObj.Location(:,1);
            covis_tracks(:,2, cov_indx) =  keyPointsObj.Location(:,2);

            % -- compute and store the depths of the features if not using monocular camera --
            if use_monocular == 0
                ref_img_r = imInit_in_SURF(frR_FileName);
                % -- compute the disparity from stereo imgs using SGM --
                disparityRange = [0 64];
                disparityMap = disparitySGM(ref_img, ref_img_r, 'DisparityRange',disparityRange,'UniquenessThreshold',20);
                disparityMap(isnan(disparityMap)) = -1;
                depth = (0.54*K(1,1)) ./ disparityMap;
                for d = 1:f_size
                    x = round(keyPointsObj.Location(d,1));
                    y = round(keyPointsObj.Location(d,2));
                    if depth(y, x) > 0
                        covis_tracks(d,3, cov_indx) = depth(y, x);
                    end
                end
            end

            % -- store gamma_0, the first point is the same between observation and gamma track --
            covis_tracks(:,4, cov_indx) = covis_tracks(:,1, cov_indx);
            covis_tracks(:,5, cov_indx) = covis_tracks(:,2, cov_indx);

            % -- retreive R_0 and T_0 at world coordinate --
            R_0 = trajectory_R(:,:,i+1); % -- 1 index difference between index of trajectory_R and index of frame --
            T_0 = trajectory_T(:,i+1);
            
            % -- make R_0 and T_0 be at the camera coordinate --
            R_cc(:,:,cov_indx) = eye(3);
            T_cc(:,:,cov_indx) = [0; 0; 0];

            cov_indx = cov_indx + 1;
            continue;
        else
            % -- read and initialize an image --
            I2 = imInit_in_SURF(fr1_FileName);

            % -- detect and extract features --
            if use_ORB
                keyPointsObj_I2 = detectORBFeatures(I2);
                [f2_Features, vpts2] = extractFeatures(I2, keyPointsObj_I2);
                f2 = f2_Features.Features;
            else
                keyPointsObj_I2 = detectSURFFeatures(I2);
                [f2, vpts2] = extractFeatures(I2, keyPointsObj_I2);
            end
            
            % -- find indices and 2d points of correspondences --
            indexPairs = matchFeatures(f1,f2) ;
            matchedPoints1 = vpts1(indexPairs(:,1));
            matchedPoints2 = vpts2(indexPairs(:,2));

            % -- delete features in covis_tracks and flows in covis_observe_tracks_u, covis_observe_tracks_v that are not matched --
            covis_tracks = covis_tracks(indexPairs(:,1),:,:);
            covis_observe_tracks_u = covis_observe_tracks_u(indexPairs(:,1),:,:);   % -- observation, x --
            covis_observe_tracks_v = covis_observe_tracks_v(indexPairs(:,1),:,:);   % -- observation, y --
            f2 = f2(indexPairs(:,2),:);
            vpts2 = vpts2(indexPairs(:,2),:);

            covis_gamma_tracks_u = covis_gamma_tracks_u(indexPairs(:,1),:,:);       % -- gamma, x --
            covis_gamma_tracks_v = covis_gamma_tracks_v(indexPairs(:,1),:,:);       % -- gamma, y --

            % -- store matched features and the displacements --
            covis_tracks(:,1, cov_indx) = matchedPoints2.Location(:,1);
            covis_tracks(:,2, cov_indx) = matchedPoints2.Location(:,2);
            covis_observe_tracks_u(:, cov_indx) = matchedPoints2.Location(:,1) - matchedPoints1.Location(:,1);
            covis_observe_tracks_v(:, cov_indx) = matchedPoints2.Location(:,2) - matchedPoints1.Location(:,2);

            % -- make R and T be at the camera coordinate --
            R = trajectory_R(:,:,i+1);
            T = trajectory_T(:,i+1);
            R_cc(:,:,cov_indx) = R_0' * R;
            T_cc(:,:,cov_indx) = R_0'*(T - T_0);
            
            if i == start_fr+1 && use_monocular == 1
                % -- compute the derivatives of R and T --
                [Rprime, Tprime, time_difference] = pose_derivatives(all_R, all_T, time_per_fr, start_fr+1);
                
                % -- compute rho_0; need a try to exclude the previous frame before start_fr --
                [covis_tracks, covis_observe_tracks_u, covis_observe_tracks_v, f2, vpts2, eq_rho_0, new_indexPairs] = get_rho(keyPointsObj_prev, f2, f_prev, vpts2, vpts_prev, Rprime, Tprime, time_difference, covis_tracks, covis_observe_tracks_u, covis_observe_tracks_v, K);
                covis_tracks(:,3, 1) = eq_rho_0(:,1);
                covis_gamma_tracks_u = covis_gamma_tracks_u(new_indexPairs(:,1),:,:);
                covis_gamma_tracks_v = covis_gamma_tracks_v(new_indexPairs(:,1),:,:);
            end

            % -- compute gamma_t --
            pix_0 = ones(3,1);
            for g = 1:size(covis_tracks, 1)
                rho_0 = covis_tracks(g,3, 1);
                if rho_0 > 0  
                    pix_0(1, 1) = covis_tracks(g,1,1);
                    pix_0(2, 1) = covis_tracks(g,2,1);
                    gamma_0 = K \ pix_0;
                    numer = rho_0.*(R_cc(:,:,cov_indx)'*gamma_0) - T_cc(:,:,cov_indx);
                    denom = rho_0.*(e_3*(R_cc(:,:,cov_indx)'*gamma_0)) - e_3*T_cc(:,:,cov_indx);
                    
                    gamma_t = numer ./ denom;
                    gamma_t_pix = (K * gamma_t);
                    covis_tracks(g,4:5, cov_indx) = gamma_t_pix(1:2,1);
                end
            end

            % -- store displacements of gamma_t --
            covis_gamma_tracks_u(:, cov_indx) = covis_tracks(:,4, cov_indx) - covis_tracks(:,4, cov_indx-1);
            covis_gamma_tracks_v(:, cov_indx) = covis_tracks(:,5, cov_indx) - covis_tracks(:,5, cov_indx-1);

            % -- replace f1 by f2 and vpts1 by vpts2 for the next round --
            f1 = f2;
            vpts1 = vpts2;
            cov_indx = cov_indx + 1;
        end
    end

    % -- delete invalid tracks and outliers (mismatched feature points and gamma(t) out of img boundary) --
    [covis_tracks, covis_observe_tracks_u, covis_observe_tracks_v, covis_gamma_tracks_u, covis_gamma_tracks_v] = remove_invalid_tracks(covis_tracks, covisible_frames, covis_observe_tracks_u, covis_observe_tracks_v, covis_gamma_tracks_u, covis_gamma_tracks_v);

    % -- compute the error per frame and the error over all frames --
    euclidean_dist = 0;
    err_over_fr = 0;
    err_per_fr = zeros(1,covisible_frames);
    for fr = 1:covisible_frames
        euclidean_dist = 0;
        for i = 1:size(covis_tracks, 1)
            x_err = (covis_tracks(i,1,fr) - covis_tracks(i,4,fr))^2;
            y_err = (covis_tracks(i,2,fr) - covis_tracks(i,5,fr))^2;
            euclidean_dist = euclidean_dist + sqrt(x_err + y_err);
        end
        err_per_fr(1,fr) = euclidean_dist / size(covis_tracks, 1);
        err_over_fr = err_over_fr + err_per_fr(1,fr);
    end

    fprintf("The average error of all frames of gamma(t) is %f\n", err_over_fr/covisible_frames);

    
    % -- analytic appraoch for computing rho_0 --
%     [analytic_covis_tracks, analytic_covis_u, analytic_covis_v, err_per_point, neg_indices] = analytic_on_rho(covis_tracks, R_cc, T_cc, K, 0, 0);
%     covis_tracks(neg_indices,:,:) = [];
%     covis_observe_tracks_u(neg_indices,:,:) = [];
%     covis_observe_tracks_v(neg_indices,:,:) = [];
%     covis_gamma_tracks_u(neg_indices,:,:) = [];
%     covis_gamma_tracks_v(neg_indices,:,:) = [];
    
    % -- analytic apporoach version 2 --
    [analytic_covis_tracks_v2, analytic_covis_u, analytic_covis_v, err_per_point_v2, neg_indices] = ...
        analytic_on_rho_v2(covis_tracks, R_cc, T_cc, K);
    covis_tracks(neg_indices,:,:) = [];
    updated_covis_tracks = analytic_covis_tracks_v2(:,1:2,:);
    covis_tracks_update_rho = covis_tracks;
    covis_tracks_update_rho(:,3,covisible_frames) = analytic_covis_tracks_v2(:,3,covisible_frames);
    covis_tracks_update_rho(:,4:5,:) = analytic_covis_tracks_v2(:,1:2,:);
    
    
    % -- optimize rho_0 using gradient descent --
    %disp("optimizing rhos ...");
    %lr = 0.1;
    %iter_num = 15000;
    %[covis_tracks_update_rho, updated_covis_tracks] = grad_optimization_on_rho(covis_tracks, lr, R_cc, T_cc, K, iter_num);
    %[covis_tracks_update_rho, updated_covis_tracks] = gd_optimization_on_rho_v2(covis_tracks, lr, R_cc, T_cc, K, iter_num);
    updated_covis_u = zeros(size(updated_covis_tracks,1), covisible_frames);
    updated_covis_v = zeros(size(updated_covis_tracks,1), covisible_frames);
    for i = 1:covisible_frames
        if i > 1
            updated_covis_u(:,i) = updated_covis_tracks(:,1,i) - updated_covis_tracks(:,1,i-1);
            updated_covis_v(:,i) = updated_covis_tracks(:,2,i) - updated_covis_tracks(:,2,i-1);
        end
    end

%     err_per_point_gd = zeros(size(updated_covis_tracks, 1), 1);
%     for i = 1:size(updated_covis_tracks, 1)
%         euclidean_dist = 0;
%         for fr = 1:covisible_frames
%             x_err = (covis_tracks(i,1,fr) - updated_covis_tracks(i,1,fr))^2;
%             y_err = (covis_tracks(i,2,fr) - updated_covis_tracks(i,2,fr))^2;
%             euclidean_dist = euclidean_dist + sqrt(x_err + y_err);
%         end
%         err_per_point_gd(i,1) = euclidean_dist / covisible_frames;
%     end
%     
%     % -- collecting data --
%     collect_rhos = zeros(size(covis_tracks, 1),4);
%     collect_rhos(:,1) = covis_tracks(:,3,1);
%     collect_rhos(:,2) = covis_tracks_update_rho(:,3,covisible_frames);
%     collect_rhos(:,3) = analytic_covis_tracks(:,3,covisible_frames);
%     collect_rhos(:,4) = analytic_covis_tracks_v2(:,3,covisible_frames);
%     err_wrt_rho = zeros(size(covis_tracks, 1),3);
%     err_wrt_rho(:,1) = err_per_point_gd(:,1);
%     err_wrt_rho(:,2) = err_per_point(:,1);
%     err_wrt_rho(:,3) = err_per_point_v2(:,1);
    
    
%     indices = find(collect_rhos(:,3)>100);
%     collect_rhos(indices,:,:) = [];
%     investigate_covis_tracks = covis_tracks;
%     investigate_covis_tracks(indices,:,:) = [];
    
    % -- INVESTIGATION: plot error function w.r.t.rho_0 --
    if debug
        figure;
        plot_error_function_rho(covis_tracks, R_cc, T_cc, K, 38);
    end

%     indices = find(covis_tracks_update_rho(:,3,covisible_frames)>150);
%     covis_tracks(indices,:,:) = [];
%     covis_observe_tracks_u(indices, :) = [];
%     covis_observe_tracks_v(indices, :) = [];
%     covis_gamma_tracks_u(indices, :) = [];
%     covis_gamma_tracks_v(indices, :) = [];
%     covis_tracks_update_rho(indices,:,:) = [];
%     updated_covis_u(indices,:,:) = [];
%     updated_covis_v(indices,:,:) = [];
%     updated_covis_tracks(indices,:,:) = [];

    % -- compute the error for each track before pose optimization --
    err_per_fr = zeros(size(covis_tracks,1),2);
    for i = 1:size(covis_tracks, 1)
        euclidean_dist = 0;
        for fr = 1:covisible_frames
            x_err = (covis_tracks(i,1,fr) - updated_covis_tracks(i,1,fr))^2;
            y_err = (covis_tracks(i,2,fr) - updated_covis_tracks(i,2,fr))^2;
            euclidean_dist = euclidean_dist + sqrt(x_err + y_err);
        end
        err_per_fr(i,1) = euclidean_dist / covisible_frames;
    end
    
    % -- optimize pose using Levenburg-Marquadt --
    disp("optimizing poses ...");
    covis_tracks_update_rho(:,4:5,:) = updated_covis_tracks(:,1:2,:);
    [covis_tracks_with_updated_poses, updated_R, updated_T, delta_R, delta_T] = grad_optimization_on_pose(covis_tracks_update_rho, R_cc, T_cc, K, 5);
    updated_pose_covis_u = zeros(size(covis_tracks_with_updated_poses,1), covisible_frames);
    updated_pose_covis_v = zeros(size(covis_tracks_with_updated_poses,1), covisible_frames);
    for i = 1:covisible_frames
        if i > 1
            updated_pose_covis_u(:,i) = covis_tracks_with_updated_poses(:,1,i) - covis_tracks_with_updated_poses(:,1,i-1);
            updated_pose_covis_v(:,i) = covis_tracks_with_updated_poses(:,2,i) - covis_tracks_with_updated_poses(:,2,i-1);
        end
    end
    
    % -- compute the error for each track after pose optimization --
    for i = 1:size(covis_tracks, 1)
        euclidean_dist = 0;
        for fr = 1:covisible_frames
            x_err = (covis_tracks(i,1,fr) - covis_tracks_with_updated_poses(i,1,fr))^2;
            y_err = (covis_tracks(i,2,fr) - covis_tracks_with_updated_poses(i,2,fr))^2;
            euclidean_dist = euclidean_dist + sqrt(x_err + y_err);
        end
        err_per_fr(i,2) = euclidean_dist / covisible_frames;
    end

    [pg_updated_R, pg_updated_T] = pose_graph_optimization(start_fr+1, all_R, all_T, updated_R, updated_T, covis_tracks, 500);
    %-*- transfer from camera coordinate to world coordinate -*-
    for f = 1:covisible_frames
       pg_updated_R(:,:,f) = R_0 * pg_updated_R(:,:,f);
       pg_updated_T(:,:,f) = R_0*pg_updated_T(:,:,f) + T_0;
    end

    % -- transfer from camera coordinate to world coordinate --
    store_updated_R = updated_R;
    store_updated_T = updated_T;
    for f = 1:covisible_frames
        store_updated_R(:,:,f) = R_0 * updated_R(:,:,f);
        store_updated_T(:,:,f) = R_0*updated_T(:,:,f) + T_0;
    end

    new_update_T = scaling(start_fr+1, store_updated_T, all_T, covis_tracks_update_rho);
    
    pg_new_update_T = scaling(start_fr+1, pg_updated_T, all_T, covis_tracks_update_rho);
    
    % -- store the optimized R and T --
    %trajectory_R(:,:,start_fr+1:start_fr+covisible_frames) = store_updated_R;
    %trajectory_T(:,start_fr+1:start_fr+covisible_frames) = new_update_T;
    %trajectory_T(:,start_fr:start_fr+covisible_frames-1) = updated_T;
    
    trajectory_R(:,:,start_fr+1:start_fr+covisible_frames) = pg_updated_R;
    trajectory_T(:,start_fr+1:start_fr+covisible_frames) = pg_new_update_T;
    %trajectory_T(:,start_fr+1:start_fr+covisible_frames) = pg_updated_T;
end

% -- plot KITTI trajectory and the differential estimations --
if plot_trajectory
    figure;
    plot(all_T(1,plot_start_fr+1:start_fr+covisible_frames), all_T(3,plot_start_fr+1:start_fr+covisible_frames), 'go', 'DisplayName', 'KITTI Ground Truths');
    hold on;
    plot(trajectory_T(1,plot_start_fr+1:start_fr+covisible_frames), trajectory_T(3,plot_start_fr+1:start_fr+covisible_frames), 'ro', 'DisplayName', 'Differential Pose Estimates');
    %hold on;
    %plot(new_update_T(1,1:covisible_frames), new_update_T(3,1:covisible_frames), 'bo', 'DisplayName', 'Scaled Differential Pose Estimates');
    %hold on;
    %plot(pg_trajectory_T(1,plot_start_fr+1:start_fr+covisible_frames), pg_trajectory_T(3,plot_start_fr+1:start_fr+covisible_frames), 'bo', 'DisplayName', 'Pose-Graph Optimized Estimates');
    legend;
    xlabel("x (m)");
    ylabel("z (m)");
    set(gcf,'color','w');
    axis equal;
end

% -- select number of plot_tracks tracks by random permutation --
%rng(2);
if plot_tracks > size(covis_tracks,1)
    plot_tracks = size(covis_tracks,1);
end
rnd_tracks = randperm(size(covis_tracks,1), plot_tracks);

% -- construct colors of the tracks --
base_colors = ["r+", "g+", "c+", "b+", "y+", "m+"];
base_q_colors = ["r", "g", "c", "b", "y", "m+"];
base_colors_dg = ["ro", "go", "co", "bo", "yo", "mo"];
n = ceil(plot_tracks / size(base_colors,2));
colors = repmat(base_colors, 1, n);
q_colors = repmat(base_q_colors, 1, n);
colors_dg = repmat(base_colors_dg, 1, n);

% -- Plot the covisible tracks and gamma(t) across the covisible frames --
if ~plot_trajectory
    figure;
    imshow(ref_img);
    for i = 1:covisible_frames
        if plot_all
            % -- plote SURF feature tracks --
            hold on;
            plot(covis_tracks(:,1,i), covis_tracks(:,2,i), 'g+');
            if i >= 2
                hold on;
                quiver(covis_tracks(:,1,i-1), covis_tracks(:,2,i-1), covis_observe_tracks_u(:,i), covis_observe_tracks_v(:,i), 0, 'g');
            end

            % -- plot gamma(t) --
            hold on;
            plot(covis_tracks(:,4,i), covis_tracks(:,5,i), 'ro');
            if i >= 2
                hold on;
                quiver(covis_tracks(:,4,i-1), covis_tracks(:,5,i-1), covis_gamma_tracks_u(:,i), covis_gamma_tracks_v(:,i), 0, 'r');
            end

            % -- plot tracks computed by updated rhos --
            hold on;
            plot(updated_covis_tracks(:,1,i), updated_covis_tracks(:,2,i), 'm*');
            if i >= 2
                hold on;
                quiver(updated_covis_tracks(:,1,i-1), updated_covis_tracks(:,2,i-1), updated_covis_u(:,i), updated_covis_v(:,i), 0, 'm');
            end

            % -- plot tracks computed by updated pose --
            hold on;
            plot(covis_tracks_with_updated_poses(:,1,i), covis_tracks_with_updated_poses(:,2,i), 'cs');
            if i >= 2
                hold on;
                quiver(covis_tracks_with_updated_poses(:,1,i-1), covis_tracks_with_updated_poses(:,2,i-1), updated_pose_covis_u(:,i), updated_pose_covis_v(:,i), 0, 'c');
            end

            % -- plot tracks computed by updated pose and repeated updated rhos --
    %         hold on;
    %         plot(new_updated_covis_tracks(:,1,i), new_updated_covis_tracks(:,2,i), 'yd');
    %         if i >= 2
    %             hold on;
    %             quiver(new_updated_covis_tracks(:,1,i-1), new_updated_covis_tracks(:,2,i-1), new_updated_covis_u(:,i), new_updated_covis_v(:,i), 0, 'y');
    %         end
        else        
            for j = 1:plot_tracks
%                 % -- plote SURF feature tracks --
%                 hold on;
%                 plot(covis_tracks(j,1,i), covis_tracks(j,2,i), 'g+');
%                 if i >= 2
%                     hold on;
%                     quiver(covis_tracks(j,1,i-1), covis_tracks(j,2,i-1), covis_observe_tracks_u(j,i), covis_observe_tracks_v(j,i), 0, 'g');
%                 end

                % -- plot gamma(t) --
                hold on;
                plot(covis_tracks(j,4,i), covis_tracks(j,5,i), 'ro');
                if i >= 2
                    hold on;
                    quiver(covis_tracks(j,4,i-1), covis_tracks(j,5,i-1), covis_gamma_tracks_u(j,i), covis_gamma_tracks_v(j,i), 0, 'r');
                end

%                 % -- plot tracks computed by updated rhos --
%                 hold on;
%                 plot(updated_covis_tracks(j,1,i), updated_covis_tracks(j,2,i), 'm*');
%                 if i >= 2
%                     hold on;
%                     quiver(updated_covis_tracks(j,1,i-1), updated_covis_tracks(j,2,i-1), updated_covis_u(j,i), updated_covis_v(j,i), 0, 'm');
%                 end
%                 
%                 % -- plot tracks computed by analytic rhos --
%                 hold on;
%                 plot(analytic_covis_tracks(j,1,i), analytic_covis_tracks(j,2,i), 'yd');
%                 if i >= 2
%                     hold on;
%                     quiver(analytic_covis_tracks(j,1,i-1), analytic_covis_tracks(j,2,i-1), analytic_covis_u(j,i), analytic_covis_v(j,i), 0, 'y');
%                 end

                % -- plot tracks computed by updated pose --
%                 hold on;
%                 plot(covis_tracks_with_updated_poses(j,1,i), covis_tracks_with_updated_poses(j,2,i), 'cs');
%                 if i >= 2
%                     hold on;
%                     quiver(covis_tracks_with_updated_poses(j,1,i-1), covis_tracks_with_updated_poses(j,2,i-1), updated_pose_covis_u(j,i), updated_pose_covis_v(j,i), 0, 'c');
%                 end

                % -- plot tracks computed by updated pose and repeated updated rhos --
    %             hold on;
    %             plot(new_updated_covis_tracks(j,1,i), new_updated_covis_tracks(j,2,i), 'yd');
    %             if i >= 2
    %                 hold on;
    %                 quiver(new_updated_covis_tracks(j,1,i-1), new_updated_covis_tracks(j,2,i-1), new_updated_covis_u(j,i), new_updated_covis_v(j,i), 0, 'y');
    %             end
    
                % -- plot gamma(t) with scaled T and scaled rho_0 --
%                 hold on;
%                 plot(scaled_covis_tracks(j,1,i), scaled_covis_tracks(j,2,i), 'yd');
%                 if i >= 2
%                     hold on;
%                     quiver(scaled_covis_tracks(j,1,i-1), scaled_covis_tracks(j,2,i-1), scaled_covis_u(j,i), scaled_covis_v(j,i), 0, 'y');
%                 end
            end
        end
    end
end

% figure;
% imshow(ref_img);
% hold on;
% plot(covis_tracks(:,1,1), covis_tracks(:,2,1), 'g+');
% text(covis_tracks(:,1,1), covis_tracks(:,2,1), string(1:size(covis_tracks, 1)), 'Color', 'm', 'FontSize', 10);

% -- the legend on the graph --
% x = linspace(0,10);
% y1 = sin(x);
% y2 = sin(0.9*x);
% y3 = sin(0.8*x);
% y4 = sin(0.7*x);
% 
% plot(x,y1,'g-+','DisplayName','SURF Tracks')
% hold on
% plot(x,y2,'r-o','DisplayName','\gamma(t) by \rho_0 computed in eq. (54)')
% plot(x,y3,'m-*','DisplayName','\gamma(t) by gradient descent optimized \rho_0')
% plot(x,y4,'y-d','DisplayName','\gamma(t) by analytically solved \rho_0')
% hold off
% 
% lgd = legend;
% lgd.FontSize = 14;

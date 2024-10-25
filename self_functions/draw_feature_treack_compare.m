clear;
close all;

datasetName = 'KITTI';
sequenceName = '08';
category_left = '/image_1/'; %the images we are actually using
[all_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, category_left);

start_fr = 69; %> must greater than 1
projectedImgName = imgFileName(start_fr,1);

covisible_feature_tracks_sift = load('seq8_SIFT_original_tracks.mat', 'covisible_feature_tracks').covisible_feature_tracks;
reprojected_feature_tracks_sift = load('seq8_SIFT_reprojected_tracks.mat', 'reprojected_feature_tracks').reprojected_feature_tracks;

covisible_feature_tracks_surf = load('seq8_SURF_original_tracks.mat', 'covisible_feature_tracks').covisible_feature_tracks;
reprojected_feature_tracks_surf = load('seq8_SURF_reprojected_tracks.mat', 'reprojected_feature_tracks').reprojected_feature_tracks;

covisible_feature_tracks_orb = load('seq8_orb_original_tracks.mat', 'covisible_feature_tracks').covisible_feature_tracks;
reprojected_feature_tracks_orb = load('seq8_orb_reprojected_tracks.mat', 'reprojected_feature_tracks').reprojected_feature_tracks;

covisible_frames = 25;
colors = lines(6);

figure;
imshow(projectedImgName);

hold on;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIFT %%%%%%%%%%%%%%%%%%%
num_tracks_sift = size(covisible_feature_tracks_sift, 1);
% -- compute u_vec and v_vec of SIFT feature tracks --
u_vec_feature_sift = zeros(size(covisible_feature_tracks_sift, 1), covisible_frames);
v_vec_feature_sift = zeros(size(covisible_feature_tracks_sift, 1), covisible_frames);
u_vec_feature_reprojected_sift = zeros(size(reprojected_feature_tracks_sift, 1), covisible_frames);
v_vec_feature_reprojected_sift = zeros(size(reprojected_feature_tracks_sift, 1), covisible_frames);

num_tracks_surf = size(covisible_feature_tracks_surf, 1);
% -- compute u_vec and v_vec of SURF feature tracks --
u_vec_feature_surf = zeros(size(covisible_feature_tracks_surf, 1), covisible_frames);
v_vec_feature_surf = zeros(size(covisible_feature_tracks_surf, 1), covisible_frames);
% -- compute u_vec and v_vec of SURF reprojected feature tracks --
u_vec_feature_reprojected_surf = zeros(size(reprojected_feature_tracks_surf, 1), covisible_frames);
v_vec_feature_reprojected_surf = zeros(size(reprojected_feature_tracks_surf, 1), covisible_frames);

num_tracks_orb = size(covisible_feature_tracks_orb, 1);
% -- compute u_vec and v_vec of feature tracks --
u_vec_feature_orb = zeros(size(covisible_feature_tracks_orb, 1), covisible_frames);
v_vec_feature_orb = zeros(size(covisible_feature_tracks_orb, 1), covisible_frames);
% -- compute u_vec and v_vec of reprojected feature tracks --
u_vec_feature_reprojected_orb = zeros(size(reprojected_feature_tracks_orb, 1), covisible_frames);
v_vec_feature_reprojected_orb = zeros(size(reprojected_feature_tracks_orb, 1), covisible_frames);

for j = 2:covisible_frames

    u_vec_feature_sift(:, j) = covisible_feature_tracks_sift(:, 1, j) - covisible_feature_tracks_sift(:, 1, j-1);
    v_vec_feature_sift(:, j) = covisible_feature_tracks_sift(:, 2, j) - covisible_feature_tracks_sift(:, 2, j-1);
    u_vec_feature_reprojected_sift(:, j) = reprojected_feature_tracks_sift(:, 1, j) - reprojected_feature_tracks_sift(:, 1, j-1);
    v_vec_feature_reprojected_sift(:, j) = reprojected_feature_tracks_sift(:, 2, j) - reprojected_feature_tracks_sift(:, 2, j-1);

    u_vec_feature_surf(:, j) = covisible_feature_tracks_surf(:, 1, j) - covisible_feature_tracks_surf(:, 1, j-1);
    v_vec_feature_surf(:, j) = covisible_feature_tracks_surf(:, 2, j) - covisible_feature_tracks_surf(:, 2, j-1);
    u_vec_feature_reprojected_surf(:, j) = reprojected_feature_tracks_surf(:, 1, j) - reprojected_feature_tracks_surf(:, 1, j-1);
    v_vec_feature_reprojected_surf(:, j) = reprojected_feature_tracks_surf(:, 2, j) - reprojected_feature_tracks_surf(:, 2, j-1);

    u_vec_feature_orb(:, j) = covisible_feature_tracks_orb(:, 1, j) - covisible_feature_tracks_orb(:, 1, j-1);
    v_vec_feature_orb(:, j) = covisible_feature_tracks_orb(:, 2, j) - covisible_feature_tracks_orb(:, 2, j-1);
    u_vec_feature_reprojected_orb(:, j) = reprojected_feature_tracks_orb(:, 1, j) - reprojected_feature_tracks_orb(:, 1, j-1);
    v_vec_feature_reprojected_orb(:, j) = reprojected_feature_tracks_orb(:, 2, j) - reprojected_feature_tracks_orb(:, 2, j-1);

end

   

% Plotting feature tracks with different colors
for i = 1:covisible_frames
    % figure;
    % imshow(imgFileName(start_fr-1+i,1));
    % hold on;  

    for k = 1:num_tracks_sift
        % Check if track should be plotted (not all zeros)
        if any(covisible_feature_tracks_sift(k, :, i) ~= 0)
            plot(covisible_feature_tracks_sift(k, 1, i), covisible_feature_tracks_sift(k, 2, i), 'g+', 'Color', colors(1, :));
            plot(reprojected_feature_tracks_sift(k, 1, i), reprojected_feature_tracks_sift(k, 2, i), 'o-', 'Color', colors(2, :));
            
            if i >= 2
                quiver(covisible_feature_tracks_sift(k, 1, i-1), covisible_feature_tracks_sift(k, 2, i-1), u_vec_feature_sift(k, i), v_vec_feature_sift(k, i), 0, 'Color', colors(1, :));
                quiver(reprojected_feature_tracks_sift(k, 1, i-1), reprojected_feature_tracks_sift(k, 2, i-1), u_vec_feature_reprojected_sift(k, i), v_vec_feature_reprojected_sift(k, i), 0, 'Color', colors(2, :));
            end
        end
    end

    hold on;  % Move hold on outside loop for better performance

    for k = 1:num_tracks_surf
        % Check if track should be plotted (not all zeros)
        if any(covisible_feature_tracks_surf(k, :, i) ~= 0)
            plot(covisible_feature_tracks_surf(k, 1, i), covisible_feature_tracks_surf(k, 2, i), 'g+', 'Color', colors(3, :));
            plot(reprojected_feature_tracks_surf(k, 1, i), reprojected_feature_tracks_surf(k, 2, i), 'o-', 'Color', colors(4, :));
            if i >= 2
                quiver(covisible_feature_tracks_surf(k, 1, i-1), covisible_feature_tracks_surf(k, 2, i-1), u_vec_feature_surf(k, i), v_vec_feature_surf(k, i), 0, 'Color', colors(3, :));
                quiver(reprojected_feature_tracks_surf(k, 1, i-1), reprojected_feature_tracks_surf(k, 2, i-1), u_vec_feature_reprojected_surf(k, i), v_vec_feature_reprojected_surf(k, i), 0, 'Color', colors(4, :));
            end
        end
    end

    hold on;

    for k = 1:num_tracks_orb
        % Check if track should be plotted (not all zeros)
        if any(covisible_feature_tracks_orb(k, :, i) ~= 0)
            plot(covisible_feature_tracks_orb(k, 1, i), covisible_feature_tracks_orb(k, 2, i), 'g+', 'Color', colors(5, :));
            plot(reprojected_feature_tracks_orb(k, 1, i), reprojected_feature_tracks_orb(k, 2, i), 'o-', 'Color', colors(6, :));

            if i >= 2
                quiver(covisible_feature_tracks_orb(k, 1, i-1), covisible_feature_tracks_orb(k, 2, i-1), u_vec_feature_orb(k, i), v_vec_feature_orb(k, i), 0, 'Color', colors(5, :));
                quiver(reprojected_feature_tracks_orb(k, 1, i-1), reprojected_feature_tracks_orb(k, 2, i-1), u_vec_feature_reprojected_orb(k, i), v_vec_feature_reprojected_orb(k, i), 0, 'Color', colors(6, :));
            end
        end
    end

    hold on;

    % Manually create dummy plots for legend
    dummy_sift_original = plot(NaN, NaN, 'g+', 'Color', colors(1, :));
    dummy_sift_reprojected = plot(NaN, NaN, 'o-', 'Color', colors(2, :));
    dummy_surf_original = plot(NaN, NaN, 'g+', 'Color', colors(3, :));
    dummy_surf_reprojected = plot(NaN, NaN, 'o-', 'Color', colors(4, :));
    dummy_orb_original = plot(NaN, NaN, 'g+', 'Color', colors(5, :));
    dummy_orb_reprojected = plot(NaN, NaN, 'o-', 'Color', colors(6, :));
    
    % Add legend to show color association
    legend([dummy_sift_original, dummy_sift_reprojected, dummy_surf_original, dummy_surf_reprojected, dummy_orb_original, dummy_orb_reprojected], ...
           {'SIFT Original', 'SIFT Reprojected', 'SURF Original', 'SURF Reprojected', 'ORB Original', 'ORB Reprojected'}, ...
           'Location', 'Best');

    pause(1);

end




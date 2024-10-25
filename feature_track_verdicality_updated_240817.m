%%%%%%%%%%%%%%%%%%%%%%%%% Initial Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%> Feature track veridicality
clc; clear all; close all;
% Setup VLFeat
run('/users/zqiwu/Desktop/vlfeat-0.9.21/toolbox/vl_setup.m');
%%%%%%%%%%%%%%%%%%%%%%%%% Initial Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%> define parameters

%%%%%%%%%%% KITTI   %%%%%%%%%%%%%%%%%%
datasetName = 'EuRoC';
sequenceName = 'MH_01_easy';
K = [718.856, 0, 607.193; 0, 718.86, 185.216; 0, 0, 1]; %> intrinsic matrix
category_left = '/image_1/'; %the images we are actually using
category_right = '/image_0/';

%{
%%%%%%%%%%% EuRoC %%%%%%%%%%%%%%%%%%%
datasetName = 'EuRoC';
sequenceName = 'MH_01_easy';
category_left = '/cam0/'; %the images we are actually using
category_right = '/cam1/';
K = [458.654, 0, 367.215; 0, 457.296, 248.375; 0, 0, 1]; %> intrinsic matrix
%}
featureType = "ORB";%"vlfeat-SIFT";
filePostFix = '.txt';

K_ = reshape(K', [1,9]);
invK = inv(K);

%> vlfeat parameters
vlfeat_params = [];
vlfeat_params.SiftThreshPeak = 10;                              % SIFT parameter
vlfeat_params.SiftThreshEdge = 20;                             % SIFT parameter
vlfeat_params.ratioTestRatio = 1.33;                            % SIFT parameter, 1/0.75 = 1.333

%> Settings
PARAMS.INLIER_REPROJ_ERR_THRESH = 2;

PARAMS.DISPLAY_MATCHES_FOR_DEBUG     = 0;
PARAMS.DRAW_INLIER_RATIO_TRACKLENGTH = 0;
PARAMS.DRAW_TRACKS = 1;
PARAMS.DRAW_EPIPOLAR_LINE = 0;

%> Function Macro for constructing a skew-symmetric matrix
skew_T = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];

%%%%%%%%%%%%%%%%%%%%%%%%% Dataset and Image Handling %%%%%%%%%%%%%%%%%%%

%> read all R, T, time stamp, and image directories of each frame, both
%> left and right at 05 both have 2761 images
[all_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, category_left);
%[~, ~, ~, imgFileName_Right] = readDataset(datasetName, sequenceName, category_right);

seq_start = 303; %> must greater than 1
seq_end = 304; 

%%%%%%%%%%%%%%%%%%%%%%%% Feature Tracking Initialization %%%%%%%%%%%%%%%%%%%%%%%%
%> initialize tracks
covisible_frames = 25;
use_monocular = 0;
if featureType == 'vlfeat-SIFT'
    collection_featureType = 'SIFT';
else
    collection_featureType = featureType;
end

currentFileDirectory = pwd + "/";
datasequenceName = datasetName +  "-seq"  + sequenceName  +"/";
folderName = strcat(currentFileDirectory, datasequenceName);
status = mkdir(folderName);

if status == 1
    disp(['Folder "', folderName, '" was created successfully.']);
else
    disp(['Failed to create folder "', folderName, '".']);
end

fileName = 'num-of-tracks-per-length-';
fullOutputFileName = strcat(folderName, fileName, collection_featureType, filePostFix);
disp(fullOutputFileName);
printToFile = fopen(fullOutputFileName, 'w');

fileName = 'num-of-verdical-tracks-per-length-';
fullOutputFileName = strcat(folderName, fileName, collection_featureType, filePostFix);
disp(fullOutputFileName);
verdical_track_printToFile = fopen(fullOutputFileName, 'w');

reproj_error_fileName = 'reproj_error_';
reproj_error_fullOutputFileName = strcat(folderName, reproj_error_fileName, collection_featureType, filePostFix);
disp(reproj_error_fullOutputFileName);
reproj_error_printToFile = fopen(reproj_error_fullOutputFileName, 'w');

inlier_ratio_fileName = 'verdical_Ratio_';
inlier_ratio_fullOutputFileName = strcat(folderName, inlier_ratio_fileName, collection_featureType, filePostFix);
disp(inlier_ratio_fullOutputFileName);
inlier_ratio_printToFile = fopen(inlier_ratio_fullOutputFileName, 'w');


%> loop from minimal covisible frames to maximal covisible frames
veridical_ratio = zeros(covisible_frames-1,1);
window_idx = 0;

%> collect numbers of N covisible frames
collect_current_num_of_track_length = zeros(covisible_frames-1, 1);
collect_current_verdical_num_of_track_length = zeros(covisible_frames-1, 1);
collect_total_num_of_track_length = zeros(size(seq_start:seq_end, 2), covisible_frames-1);


%> set camera coordinate of R and T
R_cc = zeros(3,3,covisible_frames);
T_cc = zeros(3,1,covisible_frames);
round = 0;
common_inliers = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Processing Loop %%%%%%%%%%%%%%%%%%%%%%%%%%

% The outer loop iterates over the entire sequence of images
for start_fr = seq_start:seq_end
    window_idx = window_idx + 1;
    fprintf(". ");
    if mod(start_fr, 100) == 0
        fprintf("round %d ", round);
        fprintf("\n");
        round = round + 1;
    end

    reprojection_error = ones(covisible_frames-1, 1);

    Rs = [];
    Ts = [];

    % The inner loop processes a window of frames of size 25 for feature tracking.
    for i = start_fr:start_fr+covisible_frames-1
        %> detect and extract features from 25 frames, starting from start_fr

        if i == start_fr
            [featureDescriptions, featurePoints, ~, imgCols] = self_feature_extractor(imgFileName(i,1), featureType, i, start_fr, datasetName, K);
            cov_indx = 1;
            f1_des = featureDescriptions;
            f1_pts = featurePoints;

            %> initialize covisible feature tracks
            covisible_feature_tracks = zeros(size(f1_pts, 1), 2, covisible_frames);
            reprojected_feature_tracks = zeros(size(f1_pts, 1), 2, covisible_frames);
            track_veridicality = ones(size(f1_pts, 1), covisible_frames);


            %> Matched feature points are used to update the feature tracks 
            if strcmp(featureType, 'vlfeat-SIFT') || strcmp(featureType, 'TAPIR') || strcmp(featureType, 'TEEATO')
                covisible_feature_tracks(:,:,cov_indx) = f1_pts(:,1:2);
            else
                covisible_feature_tracks(:,:,cov_indx) = f1_pts.Location;
            end

            %> initialize covisible gamma tracks
            covisible_gamma_tracks = ones(size(f1_pts, 1), 3, covisible_frames);
            covisible_gamma_tracks(:,1:2,cov_indx) = covisible_feature_tracks(:,1:2,cov_indx);

            %> retreive R_0 and T_0 at world coordinate
            R_0 = all_R(:,:,i);
            T_0 = all_T(:,i);
            R1 = R_0;
            T1 = T_0;

            Rs = [Rs, reshape(R1', [1,9])];
            Ts = [Ts; T1];

            %> make R_0 and T_0 be at the camera coordinate
            R_cc(:,:,cov_indx) = eye(3);
            T_cc(:,:,cov_indx) = [0; 0; 0];

            cov_indx = cov_indx + 1;
            continue;
        else
            
            [featureDescriptions, featurePoints, ~, ~] = self_feature_extractor(imgFileName(i,1), featureType, i, start_fr, datasetName, K);
            
            f2_des = featureDescriptions;
            f2_pts = featurePoints;

            %> do feature matching with the first image in the covisible frames
            [matchedPoints1, matchedPoints2, matches] = self_get_feature_correspondences ...
            (featureType, f1_pts, f1_des, f2_pts, f2_des, cov_indx, common_inliers);

            % ================== FEATURE TRACKS ===================
            %> Matched feature points are used to update the feature tracks 
            %> delete unmatched features from prior frames
            covisible_feature_tracks = covisible_feature_tracks(matches(:, 1),:,:);
            covisible_gamma_tracks = covisible_gamma_tracks(matches(:, 1),:,:);
            reprojected_feature_tracks = reprojected_feature_tracks(matches(:, 1),:,:);
            track_veridicality = track_veridicality(matches(:, 1),:);
            f2_des = f2_des(matches(:, 2),:);
            f2_pts = f2_pts(matches(:, 2),:);
            %> store the current matched features
            covisible_feature_tracks(:,:,cov_indx) = matchedPoints2;
            covisible_gamma_tracks(:,1:2,cov_indx)  = matchedPoints2;
          

            %> record number of covisible tracks
            collect_current_num_of_track_length(cov_indx-1, 1) = size(covisible_feature_tracks, 1);

            
            %> Forward direction: find inliers of the current frame using 
            %  the epipolar line constructed from the previous frame
            %> (i) Compute the Relative Pose
            R2 = all_R(:,:,i);
            T2 = all_T(:,i);

            Rs = [Rs, reshape(R2', [1,9])];
            Ts = [Ts; T2];

            Rel_R = R2' * R1;
            Rel_T = R2' * (T1 - T2);
            E21 = skew_T(Rel_T) * Rel_R;
            F21 = invK' * E21 * invK;
            x1 = matchedPoints1';
            
            %> (ii) Compute the epipolar line coefficients
            Apixel_21 = F21(1,:) * [x1; ones(1, size(x1, 2))];
            Bpixel_21 = F21(2,:) * [x1; ones(1, size(x1, 2))];
            Cpixel_21 = F21(3,:) * [x1; ones(1, size(x1, 2))];

            %> (iii) Find the inliers
            matchedPoints2 = matchedPoints2';
            A_xi  = Apixel_21.*matchedPoints2(1,:);
            B_eta = Bpixel_21.*matchedPoints2(2,:);
            numerOfDist = abs(A_xi + B_eta + Cpixel_21);
            denomOfDist = Apixel_21.^2 + Bpixel_21.^2;
            denomOfDist = sqrt(denomOfDist);
            dist21 = numerOfDist./denomOfDist;
            inlier_indices_forward = find(dist21 <= PARAMS.INLIER_REPROJ_ERR_THRESH);
            num_inliers_forward = length(inlier_indices_forward);
            
            %> Backward direction: find inliers of the previous frame using
            %  the epipolar line constructed by the current frame
            %> (i) Compute the Relative Pose
            Rel_R = R1' * R2;
            Rel_T = R1' * (T2 - T1);
            E12 = skew_T(Rel_T) * Rel_R;
            F12 = invK' * E12 * invK;
            x2 = matchedPoints2;
            
            %> (ii) Compute the epipolar line coefficients
            Apixel_12 = F12(1,:) * [x2; ones(1, size(x1, 2))];
            Bpixel_12 = F12(2,:) * [x2; ones(1, size(x1, 2))];
            Cpixel_12 = F12(3,:) * [x2; ones(1, size(x1, 2))];
            
            %> (iii) Find the inliers
            matchedPoints1 = matchedPoints1';
            A_xi  = Apixel_12.*matchedPoints1(1,:);
            B_eta = Bpixel_12.*matchedPoints1(2,:);
            numerOfDist = abs(A_xi + B_eta + Cpixel_12);
            denomOfDist = Apixel_12.^2 + Bpixel_12.^2;
            denomOfDist = sqrt(denomOfDist);
            dist12 = numerOfDist./denomOfDist;
            inlier_indices_backward = find(dist12 <= PARAMS.INLIER_REPROJ_ERR_THRESH);
            num_inliers_backward = length(inlier_indices_backward);


            % Find common inliers between forward and backward directions
            % Update the feature tracks to only include common inliers
            common_inliers = intersect(inlier_indices_forward, inlier_indices_backward);
            % covisible_feature_tracks = covisible_feature_tracks(common_inliers, :, :);
            % covisible_gamma_tracks = covisible_gamma_tracks(common_inliers, :, :);
            % reprojected_feature_tracks = reprojected_feature_tracks(common_inliers, :, :);
            all_indices = 1:size(track_veridicality, 1);
            non_common_inliers = setdiff(all_indices, common_inliers);
            track_veridicality(non_common_inliers, :) = 0;
           

            if PARAMS.DISPLAY_MATCHES_FOR_DEBUG == 1
                %> Display settings
                disp_Img_SideBySide = 0;
                disp_Match_Lines    = 0;
                disp_Epipolar_Lines = 1;
                disp_With_Text      = 0;
                disp_Num_Inliers = min(20, num_inliers_forward);
                img1 = imread(imgFileName(i-1,1));
                img2 = imread(imgFileName(i,  1));
                %veridical_point_1 = 
                % self_draw_epipolar_lines(img1, img2, matchedPoints1, matchedPoints2, F21, PARAMS.INLIER_REPROJ_ERR_THRESH, ...
                %                          disp_Img_SideBySide, disp_Match_Lines, disp_Epipolar_Lines, disp_Num_Inliers, disp_With_Text);
                self_draw_epipolar_lines(img1, img2, matchedPoints1, matchedPoints2, F21, PARAMS.INLIER_REPROJ_ERR_THRESH, ...
                                          disp_Img_SideBySide, disp_Match_Lines, disp_Epipolar_Lines, disp_Num_Inliers, disp_With_Text);
                
            end



            
            [feature_track_num,~] = size(covisible_feature_tracks);
            feature_track_original = [];
            Feature_Track = [];
            tot_reproj_err = 0;
            feature_track_num_valid_from_nview = 1;
            for feature_track_idx = 1: feature_track_num
                %only run NView triangulation on veridical tracks
                if(track_veridicality(feature_track_idx, cov_indx) == 1)
                    Ps = squeeze(covisible_gamma_tracks(feature_track_idx, :, 1:cov_indx));
                    Feature_Track = [];
                    for cov = 1:cov_indx
                        p = Ps(:, cov);
                        Feature_Track = [Feature_Track, p(1:2, :)];
                        feature_track_original = [feature_track_original; Ps(1:2, cov)'];
                    end
                    [corrected_features, reproj_err, is_sol_global_optimal, ~] = fast_multiview_triangulation_mex(Feature_Track, K_, Rs, Ts, 0);
                    %reprojected_feature_tracks(feature_track_idx, :, 1:cov_indx) = corrected_features;
                    %> Check if the returned reprojection errors are correct
                    total_err = 0;
                    total_idx = 0;
                    valid = 1;
                    for ci = 1:cov_indx
                        % if (abs(norm(Feature_Track(:,ci) - correced_features(:,ci)) - reproj_err(ci)) < 1e-5 && abs(reproj_err(ci)) < 4)
                        if reproj_err(ci) < 3
                            total_err = total_err + reproj_err(ci);
                            total_idx = total_idx+1;
                        else
                            valid = 0;
                        end
                    end
                    if valid==1 
                        reprojected_feature_tracks(feature_track_idx, :, 1:cov_indx) = corrected_features;
                    else
                        reprojected_feature_tracks(feature_track_idx, :, cov_indx) = zeros(2, 1);
                        track_veridicality(feature_track_idx, :) = 0;
                    end
    
                    %> Confirm that all feature points on the feature track has
                    %  reprojection error lesser than 3 pixels. If not, the
                    %  reprojected error should not be added up to tot_reproj_err.
                    if total_idx == cov_indx
                        tot_reproj_err = tot_reproj_err + total_err/total_idx;
                        feature_track_num_valid_from_nview = feature_track_num_valid_from_nview + 1;
                    end
                end
            end

            reprojection_error(cov_indx-1, 1) = tot_reproj_err / feature_track_num_valid_from_nview;
            collect_current_verdical_num_of_track_length(cov_indx-1, 1) = sum(track_veridicality(:, cov_indx)==1);
            if cov_indx==2
                veridical_ratio(cov_indx-1, window_idx) = 0;
            else
                veridical_ratio(cov_indx-1, window_idx) = sum(track_veridicality(:, cov_indx)==1) / size(track_veridicality, 1);
            end
           

            % if ~strcmp(featureType, 'TAPIR') && ~strcmp(featureType, 'TEEATO')
            %     f2_des = f2_des(common_inliers, :);
            %     f1_des = f2_des;
            % end

            %f2_pts = f2_pts(common_inliers, :);
            % Switch f2 to f1
            f1_pts = f2_pts;
            f1_des = f2_des;
            
            R1 = R2;
            T1 = T2;
            
            cov_indx = cov_indx + 1;

        end
    end
    
    collect_total_num_of_track_length(start_fr-seq_start+1, :) = collect_current_num_of_track_length(:,1)';
    
    % Write to the file simultaneously
    for i = 1:covisible_frames-1
        wr_number = collect_current_num_of_track_length(i,1);
        wr_verdical_number = collect_current_verdical_num_of_track_length(i,1);
        wr_ratio = veridical_ratio(i,window_idx);
        wr_err = reprojection_error(i, 1);
        if isnan(wr_ratio)
            wr_ratio = 0;
        end
        if isnan(wr_err)
            wr_err = 0;
        end
        wr_number = num2str(wr_number);
        wr_verdical_number = num2str(wr_verdical_number);
        wr_ratio = num2str(wr_ratio);
        wr_err = num2str(wr_err);


        fprintf(printToFile, wr_number);
        fprintf(printToFile, '\t');
        fprintf(verdical_track_printToFile, wr_verdical_number);
        fprintf(verdical_track_printToFile, '\t');
        fprintf(inlier_ratio_printToFile, wr_ratio);
        fprintf(inlier_ratio_printToFile, '\t');
        fprintf(reproj_error_printToFile, wr_err);
        fprintf(reproj_error_printToFile, '\t');

    end
    fprintf(printToFile, '\n');
    fprintf(verdical_track_printToFile, '\n');
    fprintf(inlier_ratio_printToFile, '\n');
    fprintf(reproj_error_printToFile, '\n');
end

fprintf("\n");

disp_Img_SideBySide = 0;
disp_Match_Lines    = 0;
disp_Epipolar_Lines = 1;
disp_With_Text      = 0;
%disp_Num_Inliers = min(20, num_inliers_forward);
img1 = imread(imgFileName(17,1));
img2 = imread(imgFileName(18,1));
% self_draw_epipolar_lines(img1, img2, matchedPoints1, matchedPoints2, F21, PARAMS.INLIER_REPROJ_ERR_THRESH, ...
%                          disp_Img_SideBySide, disp_Match_Lines, disp_Epipolar_Lines, disp_Num_Inliers, disp_With_Text);
% R1 = all_R(:,:,7);
% T1 = all_T(:,7);
% R2 = all_R(:,:,8);
% T2 = all_T(:,8);
% Rel_R = R2' * R1;
% Rel_T = R2' * (T1 - T2);
% E21 = skew_T(Rel_T) * Rel_R;
% F21 = invK' * E21 * invK;
% self_draw_epipolar_lines(img1, img2, covisible_feature_tracks(3, :, 7)', covisible_feature_tracks(3, :, 8)', F21, PARAMS.INLIER_REPROJ_ERR_THRESH, ...
%                           disp_Img_SideBySide, disp_Match_Lines, disp_Epipolar_Lines, disp_Num_Inliers, disp_With_Text);

% Draw feature tracks and gamma tracks
if PARAMS.DRAW_TRACKS
    disp("drawing tracks");
    self_drawColoredFeatureTracks_5_tracks(imgFileName(start_fr,1), covisible_frames, covisible_feature_tracks(:,:, :), reprojected_feature_tracks(:,:, :), -1);
    %self_drawColoredFeatureTracks_5_tracks_with_img(imgFileName,covisible_frames, start_fr, covisible_feature_tracks(2,:, :), reprojected_feature_tracks(2,:, :), -1);
end

fclose(printToFile);
fclose(inlier_ratio_printToFile);
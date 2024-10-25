%> Feature track veridicality
clc; clear all; close all;

%> define parameters
datasetName = 'KITTI';
sequenceName = '05';
category_left = '/image_0/';
category_right = '/image_1/';
featureType = 'ORB';
K = [718.856, 0, 607.193; 0, 718.86, 185.216; 0, 0, 1]; %> intrinsic matrix
invK = inv(K);

%> vlfeat parameters
vlfeat_params = [];
vlfeat_params.SiftThreshPeak = 10;                              % SIFT parameter
vlfeat_params.SiftThreshEdge = 20;                             % SIFT parameter
vlfeat_params.ratioTestRatio = 1.6;                            % SIFT parameter, 1/0.75 = 1.333

%> read all R, T, time stamp, and image directories of each frame
[all_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, category_left);
[~, ~, ~, imgFileName_Right] = readDataset(datasetName, sequenceName, category_right);

seq_start = 2; %> must greater than 1
seq_end = 2730;

%> initialize tracks
covisible_frames = 25;
%min_covisible_frames = 6;
%max_covisible_frames = 6;
e_3 = [0,0,1];
use_monocular = 0;
do_RansacEssentialMatrix = 1;
do_writeResToFile = 1;
if do_writeResToFile
    writeFileName = 'total_num_of_tracks.txt';
    currentFileDirectory = '/home/chchien/BrownU/research/Differential-Visual-Odometry/Code-22-Feb/Verifications/';
    fullOutputFileName = fullfile(currentFileDirectory, writeFileName);
    printToFile = fopen(fullOutputFileName, 'w');
end

%> draw settings
drawInlierRatio_TrackLength = 0;
drawTracks = 0;
drawEpipolarLine = 0;
debug_show = 0;

%> measure inlier ratio
reprojectionErrThr = 1;
inlier_thr = 1;

%> loop from minimal covisible frames to maximal covisible frames
bidirectional_inlier_ratio = zeros(covisible_frames-1, 2);
window_idx = 0;

%> collect numbers of N covisible frames
collect_current_num_of_track_length = zeros(covisible_frames-1, 1);
collect_total_num_of_track_length = zeros(size(seq_start:seq_end, 2), covisible_frames-1);

if do_RansacEssentialMatrix
    %> set camera coordinate of R and T
    R_cc = zeros(3,3,covisible_frames);
    T_cc = zeros(3,1,covisible_frames);
    %> essential matrix in camera coordinate
    E = zeros(3,3,covisible_frames-1);
    revE = zeros(3,3,covisible_frames-1);
end
round = 0;
for start_fr = seq_start:seq_end
    window_idx = window_idx + 1;
    fprintf(". ");
    if mod(start_fr, 100) == 0
        fprintf("round %d ", round);
        fprintf("\n");
        round = round + 1;
    end
    for i = start_fr:start_fr+covisible_frames-1
        %> detect and extract features
        [featureDescriptions, featurePoints, ~, imgCols] = self_feature_extractor(imgFileName(i,1), featureType);

        if i == start_fr
            cov_indx = 1;
            f1_des = featureDescriptions;
            f1_pts = featurePoints;

            %> initialize covisible feature tracks
            covisible_feature_tracks = zeros(size(f1_pts, 1), 2, covisible_frames);
            if strcmp(featureType, 'vlfeat-SIFT')
                covisible_feature_tracks(:,:,cov_indx) = f1_pts(:,1:2);
            else
                covisible_feature_tracks(:,:,cov_indx) = f1_pts.Location;
            end

            %> initialize covisible gamma tracks
            covisible_gamma_tracks = zeros(size(f1_pts, 1), 3, covisible_frames);
            covisible_gamma_tracks(:,1:2,cov_indx) = covisible_feature_tracks(:,1:2,cov_indx);

            if do_RansacEssentialMatrix
                %> retreive R_0 and T_0 at world coordinate
                R_0 = all_R(:,:,i);
                T_0 = all_T(:,i);
                C_0 = -R_0' * T_0;

                %> make R_0 and T_0 be at the camera coordinate
                R_cc(:,:,cov_indx) = eye(3);
                T_cc(:,:,cov_indx) = [0; 0; 0];
            end

            cov_indx = cov_indx + 1;
            continue;
        else
            [featureDescriptions, featurePoints, ~, ~] = self_feature_extractor(imgFileName(i,1), featureType);
            f2_des = featureDescriptions;
            f2_pts = featurePoints;

            %> do feature matching
            if strcmp(featureType, 'vlfeat-SIFT')
                [matches, scores] = vl_ubcmatch(f1_des', f2_des', vlfeat_params.ratioTestRatio);
                matchedPoints1 = f1_pts(matches(1, :), 1:2);
                matchedPoints2 = f2_pts(matches(2, :), 1:2);
            else
                indexPairs = matchFeatures(f1_des,f2_des);
                matchedPoints1 = f1_pts(indexPairs(:,1));   %> not useful in current scope
                matchedPoints2 = f2_pts(indexPairs(:,2));
            end

            % ================== FEATURE TRACKS ===================
            if strcmp(featureType, 'vlfeat-SIFT')
                %> delete unmatched features from prior frames
                covisible_feature_tracks = covisible_feature_tracks(matches(1, :),:,:);
                f2_des = f2_des(matches(2, :),:);
                f2_pts = f2_pts(matches(2, :),:);
                %> store the current matched features
                covisible_feature_tracks(:,:,cov_indx) = f2_pts(:,1:2);
            else
                %> delete unmatched features from prior frames
                covisible_feature_tracks = covisible_feature_tracks(indexPairs(:,1),:,:);
                f2_des = f2_des(indexPairs(:,2),:);
                f2_pts = f2_pts(indexPairs(:,2),:);
                %> store the current matched features
                covisible_feature_tracks(:,:,cov_indx) = matchedPoints2.Location;
            end

            %> record number of covisible tracks
            collect_current_num_of_track_length(cov_indx-1, 1) = size(covisible_feature_tracks, 1);

            if size(covisible_feature_tracks, 1) < 5 && window_idx > 1
                feature_size_sufficiency = 0;
                bidirectional_inlier_ratio(cov_indx-1, window_idx) = min(bidirectional_inlier_ratio(cov_indx-1, 1:window_idx-1));
            elseif size(covisible_feature_tracks, 1) < 5 && window_idx <= 1
                bidirectional_inlier_ratio(cov_indx-1, window_idx) = 1 - 1/size(covisible_feature_tracks,1);
                feature_size_sufficiency = 0;
            else
                feature_size_sufficiency = 1;
            end

            if do_RansacEssentialMatrix && feature_size_sufficiency
                %> *[Important Notes]*
                %> Should use the ground truth absolute pose to compute the
                %  relative pose, and then use fundamental matrix to decide
                %  whether features on the feature track are inliers or not
                %> Oct 30 2023: Compute relative pose of the current
                %  camera w.r.t the first camera in the local window
                R = all_R(:,:,i);
                T = all_T(:,i);
                C = -R' * T;
                Rel_R = R_0' * R;
                Rel_T = R_0' * (C_0 - C);
                Rel_Tx = [0 -Rel_T(3) Rel_T(2);
                        Rel_T(3) 0  -Rel_T(1);
                        -Rel_T(2) Rel_T(1) 0];

                %> Compute the Fundamental Matrix
                E = Rel_Tx * Rel_R;
                F = invK' * E * invK;
                
                %> *[TODO]* Use F to decide features on the feature track
                %  are inliers or not.

%                 %> Below is old version which is not the right way to do
%                 matchImg1 = ones(3, size(covisible_feature_tracks, 1));
%                 matchImg2 = ones(3, size(covisible_feature_tracks, 1));
%                 matchImg1(1:2, :) = covisible_feature_tracks(:,1:2,cov_indx-1)';
%                 matchImg2(1:2, :) = covisible_feature_tracks(:,1:2,cov_indx)';
%                 [finalE, inlierIndx_forward] = self_Ransac4Essential(5000, reprojectionErrThr, matchImg1, matchImg2, K, debug_show);
%                 E(:,:,cov_indx-1) = finalE;
%                 inlierRatio_forward = size(inlierIndx_forward,2) / size(covisible_feature_tracks, 1);
% 
%                 %> Reverse the matching and find the inliers
%                 [finalE, inlierIndx_backward] = self_Ransac4Essential(5000, reprojectionErrThr, matchImg2, matchImg1, K, debug_show);
%                 revE(:,:,cov_indx-1) = finalE;
%                 inlierRatio_backward = size(inlierIndx_backward,2) / size(covisible_feature_tracks, 1);
% 
%                 bidirectional_inlier_ratio(cov_indx-1, window_idx) = min(inlierRatio_forward, inlierRatio_backward);
            end
            %> switch f2 to f1
            f1_des = f2_des;
            f1_pts = f2_pts;
            cov_indx = cov_indx + 1;
        end
    end
    
    collect_total_num_of_track_length(start_fr-seq_start+1, :) = collect_current_num_of_track_length(:,1)';
    
    %> write to the file simultneously
    if do_writeResToFile
        for i = 1:size(collect_current_num_of_track_length, 1)
            wr_number = num2str(collect_current_num_of_track_length(i,1));
            fprintf(printToFile, wr_number);
            fprintf(printToFile, '\t');
        end
        fprintf(printToFile, '\n');
    end

end
fprintf("\n");

%> I commented out the below code 2 years ago...

%fprintf("end-point inlier ratio: %f\n", min(inlierRatio_forward, inlierRatio_backward));

%     % ================== DISTANCE TO EPIPOLAR LINE =====================
%     %> compute coefficients of an epipolar line
%     A = zeros(size(E, 3), size(covisible_feature_tracks, 1));
%     B = zeros(size(E, 3), size(covisible_feature_tracks, 1));
%     C = zeros(size(E, 3), size(covisible_feature_tracks, 1));
%     feat1 = ones(3,size(covisible_feature_tracks, 1));
%     feat2 = feat1;
%     dist12 = zeros(size(E, 3), size(feat1, 2));
%     dist21 = zeros(size(E, 3), size(feat1, 2));
%     denomOfDist = zeros(size(E, 3), size(feat1, 2));
%     numerOfDist = zeros(size(E, 3), size(feat1, 2));
%     A_ep = zeros(size(E, 3), size(feat1, 2));
%     B_it = zeros(size(E, 3), size(feat1, 2));
%     
%     %> loop over all essential matrices
%     for p = 1:covisible_frames-1
%         
%         feat1(1,:) = covisible_feature_tracks(:,1,p)';
%         feat1(2,:) = covisible_feature_tracks(:,2,p)';
%         feat2(1,:) = covisible_feature_tracks(:,1,p+1)';
%         feat2(2,:) = covisible_feature_tracks(:,2,p+1)';
%         
%         calE = inv(K)' * E(:,:,p) * inv(K);
%         A(p, :) = calE(1,:) *feat1;
%         B(p, :) = calE(2,:) *feat1;
%         C(p, :) = calE(3,:) *feat1;
%         
%         %> compute the distance from the features to the epipolar line
%         for k = 1 : size(feat1, 2)
%             A_ep(:,k) = A(:,k).*feat2(1, k);    %> epsilon
%             B_it(:,k) = B(:,k).*feat2(2, k);    %> ita
%         end
%         numerOfDist = abs(A_ep + B_it + C);
%         denomOfDist = A.^2 + B.^2;
%         denomOfDist = sqrt(denomOfDist);
%         dist12 = numerOfDist./denomOfDist;
%     end
%     
%     %> reverse directional epipolar line    
%     %> loop over all essential matrices
%     for p = 1:covisible_frames-1
%         feat1(1,:) = covisible_feature_tracks(:,1,p)';
%         feat1(2,:) = covisible_feature_tracks(:,2,p)';
%         feat2(1,:) = covisible_feature_tracks(:,1,p+1)';
%         feat2(2,:) = covisible_feature_tracks(:,2,p+1)';
%         
%         calE = inv(K)' * revE(:,:,p) * inv(K);
%         A(p, :) = calE(1,:) *feat2;
%         B(p, :) = calE(2,:) *feat2;
%         C(p, :) = calE(3,:) *feat2;
%         
%         %> compute the distance from the features to the epipolar line
%         for k = 1 : size(feat1, 2)
%             A_ep(:,k) = A(:,k).*feat1(1, k);    %> epsilon
%             B_it(:,k) = B(:,k).*feat1(2, k);    %> ita
%         end
%         numerOfDist = abs(A_ep + B_it + C);
%         denomOfDist = A.^2 + B.^2;
%         denomOfDist = sqrt(denomOfDist);
%         dist21 = numerOfDist./denomOfDist;
%     end
%     % =================================================================
%     
%     %> use distance to epipolar line to decide inlier/outlier
%     num_tracks = size(covisible_feature_tracks, 1);
%     forward_inlier_ratio = 0;
%     backward_inlier_ratio = 0;
%      
% %     for p = 1:size(dist12,1)
% %         inlier_idx_forward = find(dist12(p,:)<inlier_thr);
% %         inlier_idx_backward = find(dist21(p,:)<inlier_thr);
% %         %consistency_inlier_num = size(intersect(inlier_idx_forward, inlier_idx_backward), 2);
% %         %disp(consistency_inlier_num/num_tracks);
% %         %disp(size(inlier_idx_forward, 2)/num_tracks);
% %         %disp(size(inlier_idx_backward, 2)/num_tracks);
% %         
% %         forward_inlier_ratio = forward_inlier_ratio + size(inlier_idx_forward, 2)/num_tracks;
% % %         backward_inlier_ratio = backward_inlier_ratio + size(inlier_idx_backward, 2)/num_tracks;
% %     end
%     
%     inlier_idx_forward = find(dist12(end,:)<inlier_thr);
%     inlier_idx_backward = find(dist21(end,:)<inlier_thr);
%     %disp(size(inlier_idx_forward, 2)/size(forward_inlierIndx,2));
%     %disp(size(inlier_idx_backward, 2)/size(backward_inlierIndx,2));
%     
% %     total_inlier_ratio = forward_inlier_ratio + backward_inlier_ratio;
% %     total_inlier_ratio = total_inlier_ratio / (2*(covisible_frames-1));
% %     disp(total_inlier_ratio);
%     
%     if drawEpipolarLine
%         img1 = imread(imgFileName(start_fr,1));
%         img2 = imread(imgFileName(start_fr+1,1));
%         
%         %> compute the coefficients of the epipolar line function
%         draw_feature_index = 1;
%         a_img2 = E(1,:,1)*feat1(:, draw_feature_index);
%         b_img2 = E(2,:,1)*feat1(:, draw_feature_index);
%         c_img2 = E(3,:,1)*feat1(:, draw_feature_index);
% 
%         yMin_img2 = -c_img2./b_img2;
%         yMax_img2 = (-c_img2 - a_img2*imgCols) ./ b_img2;
% 
%         figure, img4 = [img1, img2];  
%         colormap('gray');
%         imagesc(img4);
% 
%         hold on;
%         line([imgCols, imgCols+imgCols], ...
%              [round(yMin_img2), round(yMax_img2)], 'Color', 'r');
%         plot(feat1(1,draw_feature_index), feat1(2,draw_feature_index), 'r+');
%         plot(feat2(1,draw_feature_index)+imgCols, feat2(2,draw_feature_index), 'r+');
%         hold off;
%     end


%> draw inlier ratio v.s. track lengths
% if drawInlierRatio_TrackLength
%     figure;
%     set_tracks = min_covisible_frames:max_covisible_frames;
%     plot(set_tracks, inlier_ratio(:,1), 'bo-');
%     xlabel('track length');
%     ylabel('inlier ratio');
%     legend(featureType);
%     set(gcf,'color','w');
% end

%> draw feature tracks and gamma tracks
if drawTracks
    self_drawFeatureTracks(imgFileName(start_fr,1), covisible_frames, covisible_feature_tracks, covisible_gamma_tracks, -1);
end

fclose(printToFile);

function [matchedPoints1, matchedPoints2, bidirectional_consistent_matches] = ...
    self_get_feature_correspondences(featureType, f1_pts, f1_des, f2_pts, f2_des, cov_indx, common_inliers)

    %> do feature matching with the first image in the covisible frames
    if strcmp(featureType, 'vlfeat-SIFT') 
        %> Perform forward directional feature matching
        [matches_1to2, ~] = vl_ubcmatch(f1_des', f2_des');

        %> make sure that the other direction (img2 -> img1) of matching is 
        %  consistent with the current direction (img1 -> img2)
        [matches_2to1, ~] = vl_ubcmatch(f2_des', f1_des');

        %> Find consistent matches from both directions
        matches_1to2 = matches_1to2';
        matches_2to1([1 2], :) = matches_2to1([2 1], :);    %> swap rows
        matches_2to1 = sortrows(matches_2to1');
        bidirectional_consistent_matches = intersect(matches_1to2, matches_2to1, 'rows');
        bidirectional_consistent_matches = bidirectional_consistent_matches';

        matchedPoints1 = f1_pts(bidirectional_consistent_matches(1, :), 1:2);
        matchedPoints2 = f2_pts(bidirectional_consistent_matches(2, :), 1:2);
        bidirectional_consistent_matches = bidirectional_consistent_matches';
    else
        if strcmp(featureType, 'TAPIR') || strcmp(featureType, 'TEEATO')
            matchedPoints1 = f1_pts; 
            if cov_indx==2
                matchedPoints2 = f2_pts;
            else
                matchedPoints2 = f2_pts(common_inliers, :);
            end
        else
            indexPairs_1to2 = matchFeatures(f1_des,f2_des);
            indexPairs_2to1 = matchFeatures(f2_des,f1_des);

            indexPairs_2to1(:, [1 2]) = indexPairs_2to1(:, [2 1]);    %> swap rows
            indexPairs_2to1 = sortrows(indexPairs_2to1);
            bidirectional_consistent_matches = intersect(indexPairs_1to2, indexPairs_2to1, 'rows');

            matchedPoints1 = f1_pts(bidirectional_consistent_matches(:,1)).Location;   %> not useful in current scope
            matchedPoints2 = f2_pts(bidirectional_consistent_matches(:,2)).Location;
        end
    end
end
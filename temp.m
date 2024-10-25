clear;
close all;

R1 = [1, 0, 0; 0, 1, 0; 0, 0, 1];
T1 = [0; 0; 0];
R2 = [1, 0.0005, -0.0021; -0.0005, 1, -0.0012; 0.0021, 0.0012, 1];
T2 = [-0.0469; -0.0284; 0.8587];
R3 = [1, 0.0010, -0.0041; -0.0011, 1, -0.0023; 0.0041, 0.0023, 1];
T3 = [-0.0937; -0.0568; 1.7163];

K = [718.856, 0, 607.193; 0, 718.86, 185.216; 0, 0, 1];
invK = inv(K);
%> Function Macro for constructing a skew-symmetric matrix
skew_T = @(T)[0, -T(3,1), T(2,1); T(3,1), 0, -T(1,1); -T(2,1), T(1,1), 0];

matchedPoints1 = [130.5744, 343.3321, 1]';
matchedPoints2 = [89.9844, 360.8625, 1]';
matchedPoints3 = [111.3772, 351.2311, 1]';

% 
% Rel_R = R2' * R1;
% Rel_T = R2' * (T1 - T2);
% E21 = skew_T(Rel_T) * Rel_R;
% F21 = invK' * E21 * invK;
% x1 = matchedPoints1';
% 
% %> (ii) Compute the epipolar line coefficients
% Apixel_21 = F21(1,:) * [x1; ones(1, size(x1, 2))];
% Bpixel_21 = F21(2,:) * [x1; ones(1, size(x1, 2))];
% Cpixel_21 = F21(3,:) * [x1; ones(1, size(x1, 2))];
% 
% %> (iii) Find the inliers
% matchedPoints2 = matchedPoints2';
% A_xi  = Apixel_21.*matchedPoints2(1,:);
% B_eta = Bpixel_21.*matchedPoints2(2,:);
% numerOfDist = abs(A_xi + B_eta + Cpixel_21);
% denomOfDist = Apixel_21.^2 + Bpixel_21.^2;
% denomOfDist = sqrt(denomOfDist);
% dist21 = numerOfDist./denomOfDist;
% 
% %> Backward direction: find inliers of the previous frame using
% %  the epipolar line constructed by the current frame
% %> (i) Compute the Relative Pose
% Rel_R = R1' * R2;
% Rel_T = R1' * (T2 - T1);
% E12 = skew_T(Rel_T) * Rel_R;
% F12 = invK' * E12 * invK;
% x2 = matchedPoints2;
% 
% %> (ii) Compute the epipolar line coefficients
% Apixel_12 = F12(1,:) * [x2; ones(1, size(x1, 2))];
% Bpixel_12 = F12(2,:) * [x2; ones(1, size(x1, 2))];
% Cpixel_12 = F12(3,:) * [x2; ones(1, size(x1, 2))];
% 
% %> (iii) Find the inliers
% matchedPoints1 = matchedPoints1';
% A_xi  = Apixel_12.*matchedPoints1(1,:);
% B_eta = Bpixel_12.*matchedPoints1(2,:);
% numerOfDist = abs(A_xi + B_eta + Cpixel_12);
% denomOfDist = Apixel_12.^2 + Bpixel_12.^2;
% denomOfDist = sqrt(denomOfDist);
% dist12 = numerOfDist./denomOfDist;




%%%%%%%%%%%%%%%%%%%% NVIEW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_ = reshape(K', [1,9]);
Rs = [reshape(R1', [1,9]), reshape(R3', [1,9])];
Ts = [T1; T3];
Feature_Track = [matchedPoints1(1:2,:), matchedPoints3(1:2, :)];

debug = 0;
[corrected_features, reproj_errs] = fast_multiview_triangulation_mex(Feature_Track, K_, Rs, Ts, debug);

%> Check if the returned reprojection erros are correct and making senses
for i = 1:2
    assert(abs(norm(Feature_Track(:,i) - corrected_features(:,i)) - reproj_errs(i)) < 1e-8);
    assert(abs(reproj_errs(i)) < 5);
end

%> average reprojection error
avg_reproj_err = mean(reproj_errs);


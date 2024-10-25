%% Use gradient descent optimization to fit a EuRoC dataset trajectory from a caetesian parametrized circular helix curve

clc; clear all; close all;

% -- define parameters --
datasetName = 'EuRoC';
sequenceName = 'MH_03_medium/mav0/';

% -- read R, T, and time stamp of each frame --
[every_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, '');

% -- sample the frames --
keyframe_sampling = 1;
keyfr_idx = 1;
for i = 1:size(every_T,2)
    if mod(i-1, keyframe_sampling) == 0
        all_T(:,keyfr_idx) = every_T(:,i);
        time_per_keyfr(1,keyfr_idx) = time_per_fr(1,i);
        keyfr_idx = keyfr_idx + 1;
    end
end

start_fr = 1905;
end_fr = size(all_T,2)-50;
seq_end_fr = size(all_T,2)-50;
%seq_end_fr = 2000;
window_size = 9;

run_curve_model_segment_wise = 1;

run_segment_wise = 1;
debugging = 0;

% -- dynamic model params and geometry model params --
dynamicParams = [];
geometryParams = [];

% -- drawings and printings --
print_segment_wise = run_segment_wise;
drawErrChangeAlongSequence = ~run_segment_wise;

drawSegementStartFr = start_fr;
drawCurveAndTrajectory = run_segment_wise;
drawGeometryModelFit = run_segment_wise;
drawTrajectory = 0;
drawOptimizedFrenetFrame = 0;
drawFrenetFrameOnTrajectory = 1;
drawFrenetFrameFromVelAndAcc = 0;

if drawTrajectory
    % -- plot trajectory --
    figure;
    p_ground_truth = plot3(all_T(1,start_fr:end_fr), all_T(2,start_fr:end_fr), all_T(3,start_fr:end_fr), 'bo--');
    hold on;
    plot3(all_T(1,start_fr), all_T(2,start_fr), all_T(3,start_fr), 'ko');
    text(all_T(1,start_fr), all_T(2,start_fr), all_T(3,start_fr), 'start', 'FontSize', 10, 'Color', 'k');
    hold on;
    plot3(all_T(1,end_fr), all_T(2,end_fr), all_T(3,end_fr), 'ko');
    text(all_T(1,end_fr), all_T(2,end_fr), all_T(3,end_fr), 'end', 'FontSize', 10, 'Color', 'k');
    xlabel("x (m)");
    ylabel("y (m)");
    zlabel("z (m)");
    axis equal;
    set(gcf,'color','w');
end


if run_segment_wise
    end_fr = start_fr;
else
    end_fr = seq_end_fr;
end

collection_of_local_avg_errs = zeros(size(1:end_fr-start_fr+1, 2), 1);
collect_err_percentage = zeros(size(1:end_fr-start_fr+1, 2), 2);

% -- loop over all start_fr to end_fr to compute frenet frames --
record_err_idx = 1;
for fr_idx = start_fr:end_fr

    [dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(fr_idx, fr_idx+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);

    vel = dynamicParams.vel_vec(3,:)';
    acc = dynamicParams.acc_vec(3,:)';
    vec_acc_cross_prod = cross(vel, acc);
    init_binormal = vec_acc_cross_prod / norm(vec_acc_cross_prod);
    init_tangent = vel / norm(vel);
    tangent_binormal_cross_prod = -cross(init_tangent, init_binormal);
    init_normal = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);
    
    init_FrenetFrame.T = init_tangent;
    init_FrenetFrame.N = init_normal;
    init_FrenetFrame.B = init_binormal;

    % -- valid trajectory points --
    trajectory_pts = r(:,3:end);
    
    [collect_opt_curve_pts, collect_opt_geometry, collectFrenetFrame, avg_err] = self_approx_by_cartesian_parametrization(trajectory_pts, init_FrenetFrame, datasetName);

    
end

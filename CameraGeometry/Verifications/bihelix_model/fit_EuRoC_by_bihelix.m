%% Use gradient descent optimization to fit a EuRoC dataset trajectory from a bihelical model

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
forward_tangent = 0;
for fr_idx = start_fr:end_fr

    [dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(fr_idx, fr_idx+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);

    if forward_tangent
        vel = (dynamicParams.vel_vec(3,:)' + dynamicParams.vel_vec(4,:)').*0.5;
        acc = (dynamicParams.acc_vec(3,:)' + dynamicParams.acc_vec(4,:)').*0.5;
    else
        vel = dynamicParams.vel_vec(3,:)';
        acc = dynamicParams.acc_vec(3,:)';
    end
    
    vec_acc_cross_prod = cross(vel, acc);
    init_binormal = vec_acc_cross_prod / norm(vec_acc_cross_prod);
    init_tangent = vel / norm(vel);
    tangent_binormal_cross_prod = -cross(init_tangent, init_binormal);
    init_normal = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);
    
    init_FrenetFrame.T = init_tangent;
    init_FrenetFrame.N = init_normal;
    init_FrenetFrame.B = init_binormal;

    % -- valid trajectory points and scale up --
    trajectory_pts = r(:,3:end);
    
    % compute arcLength
    cumulative_arcLength = 0;
    arcLength = zeros(size(trajectory_pts, 2),1);
    for i = 2:size(trajectory_pts, 2)
        cumulative_arcLength = cumulative_arcLength + norm(trajectory_pts(:,i)-trajectory_pts(:,i-1));
        arcLength(i,1) = cumulative_arcLength;
    end

    gdParams.lr = 5;
    %gdParams.lr = 0.5;
    gdParams.num_of_iter = 15000;
    gdParams.num_of_first_success = 100;
    gdParams.num_of_second_success = 200;
    
    min_avg_err = 1000;
    %j_pts = 4;
    err_type = "avg_err";
    %err_type = "end_pt_err";
    opt_frenet_frame = 0;
    common_tau = 1;
    
    if opt_frenet_frame
         % -- convert initial frenet frame to unit sphere angles --
         [theta1, phi1, theta2, phi2] = self_cvt_frenetframe_to_unit_sphere_angle(init_tangent, init_normal);
    end
    
    for j_pts = 3:size(trajectory_pts, 2)-2
        %fetch_idx_helix1 = round(j_pts/2)+1;
        %fetch_idx_helix2 = round(((size(trajectory_pts, 2)-j_pts)/2) + j_pts)+1;
%         initParams.init_k0 = geometryParams.curvature(3,1) / 10;
%         initParams.init_tau0 = geometryParams.torsion(3,1);
%         initParams.init_k1 = geometryParams.curvature(end,1) / 10;
%         initParams.init_tau1 = geometryParams.torsion(end,1);
        
        initParams.init_k0 = 40;
        initParams.init_tau0 = 10;
        initParams.init_k1 = 30;
        initParams.init_tau1 = -10;
        
        if opt_frenet_frame
            gdParams.lr_frenetframe = 0.0001;

            initParams.theta1 = theta1;
            initParams.phi1 = phi1;
            initParams.theta2 = theta2;
            initParams.phi2 = phi2;
           
            if common_tau
                [opt_res, collect_gd_errors, collect_learning_rate, err_percentage] = self_GD_min_curve_err_common_tau_FF(gdParams, initParams, trajectory_pts, j_pts, err_type);
                
                avg_err = collect_gd_errors(end,1);
                if avg_err < min_avg_err
                    min_avg_err = avg_err;
                    min_err_percentage = err_percentage;
                    j_pt_candidate = j_pts;
                    opt_params.k0 = opt_res.k0;
                    opt_params.tau0 = opt_res.tau0;
                    opt_params.k1 = opt_res.k1;
                    opt_params.tau1 = opt_res.tau0;

                    opt_theta1 = opt_res.theta1;
                    opt_phi1 = opt_res.phi1;
                    opt_theta2 = opt_res.theta2;
                    opt_phi2 = opt_res.phi2;
                    % -- vi) convert unit sphere angles to frenet frame --
                    [opt_T, opt_N, opt_B] = self_cvt_unit_sphere_angle_to_frenetframe(opt_theta1, opt_phi1, opt_theta2, opt_phi2);

                    opt_params.T = opt_T;
                    opt_params.N = opt_N;
                    opt_params.B = opt_B;
                    final_gd_errors = collect_gd_errors;
                    final_lr = collect_learning_rate;
                end
            else
                [opt_res, collect_gd_errors, collect_learning_rate, err_percentage] = self_GD_min_curve_err_w_var_frenet_frame(gdParams, initParams, trajectory_pts, j_pts, err_type);

                avg_err = collect_gd_errors(end,1);
                if avg_err < min_avg_err
                    min_avg_err = avg_err;
                    min_err_percentage = err_percentage;
                    j_pt_candidate = j_pts;
                    opt_params.k0 = opt_res.k0;
                    opt_params.tau0 = opt_res.tau0;
                    opt_params.k1 = opt_res.k1;
                    opt_params.tau1 = opt_res.tau1;

                    opt_theta1 = opt_res.theta1;
                    opt_phi1 = opt_res.phi1;
                    opt_theta2 = opt_res.theta2;
                    opt_phi2 = opt_res.phi2;
                    % -- vi) convert unit sphere angles to frenet frame --
                    [opt_T, opt_N, opt_B] = self_cvt_unit_sphere_angle_to_frenetframe(opt_theta1, opt_phi1, opt_theta2, opt_phi2);

                    opt_params.T = opt_T;
                    opt_params.N = opt_N;
                    opt_params.B = opt_B;
                    final_gd_errors = collect_gd_errors;
                    final_lr = collect_learning_rate;
                end
            end
        else
            %init_FrenetFrame.N = [0.4283; 0.9019; 0.0554];
            %init_FrenetFrame.B = cross(init_FrenetFrame.T, init_FrenetFrame.N);
            if common_tau
                [opt_res, collect_gd_errors, collect_learning_rate, err_percentage] = self_GD_min_curve_err_common_tau(gdParams, initParams, init_FrenetFrame, trajectory_pts, j_pts, err_type);
                
                avg_err = collect_gd_errors(end,1);
                if avg_err < min_avg_err
                    min_avg_err = avg_err;
                    min_err_percentage = err_percentage;
                    j_pt_candidate = j_pts;
                    opt_params.k0 = opt_res.k0;
                    opt_params.tau0 = opt_res.tau0;
                    opt_params.k1 = opt_res.k1;
                    opt_params.tau1 = opt_res.tau0;
                    opt_params.T = opt_res.T;
                    opt_params.N = opt_res.N;
                    opt_params.B = opt_res.B;
                    final_gd_errors = collect_gd_errors;
                    final_lr = collect_learning_rate;
                end
            else
                [opt_res, collect_gd_errors, collect_learning_rate, err_percentage] = self_GD_min_curve_err(gdParams, initParams, init_FrenetFrame, trajectory_pts, j_pts, err_type);    

                avg_err = collect_gd_errors(end,1);
                if avg_err < min_avg_err
                    min_avg_err = avg_err;
                    min_err_percentage = err_percentage;
                    j_pt_candidate = j_pts;
                    opt_params.k0 = opt_res.k0;
                    opt_params.tau0 = opt_res.tau0;
                    opt_params.k1 = opt_res.k1;
                    opt_params.tau1 = opt_res.tau1;
                    opt_params.T = opt_res.T;
                    opt_params.N = opt_res.N;
                    opt_params.B = opt_res.B;
                    final_gd_errors = collect_gd_errors;
                    final_lr = collect_learning_rate;
                end
            end
        end
    end

    [opt_curve, ~] = self_generateBihelixFromModel(opt_params.k0, opt_params.tau0, opt_params.k1, opt_params.tau1, opt_params.T, opt_params.N, opt_params.B, trajectory_pts, j_pt_candidate);
    
    % -- print out the optimal values --
    fprintf("Optimal Values:");
    fprintf("\nk0: ");
    fprintf(string(opt_params.k0));
    fprintf("\ntau0: ");
    fprintf(string(opt_params.tau0));
    fprintf("\nk1: ");
    fprintf(string(opt_params.k1));
    fprintf("\ntau1:");
    fprintf(string(opt_params.tau1));
    fprintf("\njuntion point: ");
    fprintf(string(j_pt_candidate));
    if strcmp(err_type, "end_pt_err")
        % -- compute the average error --
        energy_val = 0;
        for i = 1:size(trajectory_pts,2)
            energy_val = energy_val + norm(trajectory_pts(:,i)-opt_curve(:,i));
        end
        energy_val = energy_val / size(trajectory_pts,2);
        err_percentage = (energy_val / arcLength(end,1))*100;
    
        fprintf("\nend-pt err: ");
        fprintf(string(min_avg_err));
        fprintf("\navg-pt err: ");
        fprintf(string(energy_val));
        fprintf("\nerr percentage: ");
        fprintf(string(err_percentage));
    else
        fprintf("\navg-pt err: ");
        fprintf(string(min_avg_err));
        fprintf("\nerr percentage: ");
        fprintf(string(min_err_percentage));
    end

    fprintf("\n");
    
    % -- PLOTTINGS --
    r1 = trajectory_pts';
    r2 = opt_curve';
    figure;
    plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'ground truth');
    hold on;
    plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', 'optimal curve');
    hold on;
    
    drawTangent = opt_params.T;
    drawNormal = opt_params.N;
    drawBinormal = opt_params.B;
    scale_factor = 0.005;
    drawTangent = drawTangent.*scale_factor;
    drawNormal = drawNormal.*scale_factor;
    drawBinormal = drawBinormal.*scale_factor;

    % -- draw frenet frame at the initial point --
    endpoint_T = r1(1,:) + drawTangent';
    endpoint_N = r1(1,:) + drawNormal';
    endpoint_B = r1(1,:) + drawBinormal';

    line([r1(1,1), endpoint_T(1,1)],[r1(1,2), endpoint_T(1,2)],[r1(1,3), endpoint_T(1,3)],'Color','g','LineStyle','-', 'DisplayName', 'tangent');
    hold on;
    line([r1(1,1), endpoint_N(1,1)],[r1(1,2), endpoint_N(1,2)],[r1(1,3), endpoint_N(1,3)],'Color','c','LineStyle','-', 'DisplayName', 'normal');
    hold on;
    line([r1(1,1), endpoint_B(1,1)],[r1(1,2), endpoint_B(1,2)],[r1(1,3), endpoint_B(1,3)],'Color','m','LineStyle','-', 'DisplayName', 'binormal');
    hold on;
    
    legend;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis equal;
    set(gcf,'color','w');
end

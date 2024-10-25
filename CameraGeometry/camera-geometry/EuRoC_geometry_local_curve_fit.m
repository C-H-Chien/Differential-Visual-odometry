
% -- Verify that the camera center model C(t) is a sufficient model using KITTI dataset --
clc; clear all; close all;

% -- define parameters --
datasetName = 'EuRoC';
sequenceName = 'MH_02_easy/mav0/';
verification_class = 'geometry';

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

start_fr = 1281;
end_fr = size(all_T,2)-50;
seq_end_fr = size(all_T,2)-50;
%seq_end_fr = 2000;
window_size = 9;

if strcmp(verification_class, 'geometry')
    run_curve_model_segment_wise = 1;
end

run_segment_wise = 1;
opt_FrenetFrame = 1;
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
drawFrenetFrameOnTrajectory = 1;


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

if strcmp(verification_class, 'geometry')

    if run_segment_wise
        end_fr = start_fr;
    else
        end_fr = seq_end_fr;
    end
    
    collect_propT = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_propN = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_propB = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_all_init_traj_pts = zeros(3, size(1:end_fr-start_fr+1, 2));
    
    collect_curvature = zeros(1, size(1:end_fr-start_fr+1, 2));
    collect_initial_params = zeros(2, size(1:end_fr-start_fr+1, 2));
    collect_vel_vector = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_acc_vector = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_tangent = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_binormal = zeros(3, size(1:end_fr-start_fr+1, 2));
    collect_normal = zeros(3, size(1:end_fr-start_fr+1, 2));
    collection_of_local_avg_errs = zeros(size(1:end_fr-start_fr+1, 2), 1);
    collect_err_percentage = zeros(size(1:end_fr-start_fr+1, 2), 2);
    
    % -- loop over all start_fr to end_fr to compute frenet frames --
    record_err_idx = 1;
    pt_of_frenetframe = 3; % first pt = 3
    forward_tangent = 0;
    for fr_idx = start_fr:end_fr

        [dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(fr_idx, fr_idx+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);
        
        if forward_tangent
            collect_vel_vector(:,record_err_idx) = (dynamicParams.vel_vec(4,:)+dynamicParams.vel_vec(3,:))'.*0.5;
            collect_acc_vector(:,record_err_idx) = (dynamicParams.acc_vec(4,:)+dynamicParams.acc_vec(3,:))'.*0.5;
        else
            collect_vel_vector(:,record_err_idx) = dynamicParams.vel_vec(pt_of_frenetframe,:)';
            collect_acc_vector(:,record_err_idx) = dynamicParams.acc_vec(pt_of_frenetframe,:)';
        end
        
        
        vec_acc_cross_prod = cross(collect_vel_vector(:,record_err_idx), collect_acc_vector(:,record_err_idx));
        
        % -- initial Frenet Frame --
        collect_binormal(:,record_err_idx) = vec_acc_cross_prod / norm(vec_acc_cross_prod);
        collect_tangent(:,record_err_idx) = collect_vel_vector(:,record_err_idx) / norm(collect_vel_vector(:,record_err_idx));
        tangent_binormal_cross_prod = -cross(collect_tangent(:,record_err_idx), collect_binormal(:,record_err_idx));
        collect_normal(:,record_err_idx) = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);
        
        % -- ground truth trajectory --
        % 1) ground truth points
        trajectory_pts = r(:,3:end);
        curve_start_pt = r(:,3);
        %trajectory_pts = trajectory_pts*10;
        
        % 2) arc lengths
        cumulative_arc_length = zeros(size(trajectory_pts, 2), 1);
        arcLength = zeros(size(trajectory_pts, 2), 1);
        for i = 2:size(trajectory_pts,2)
            %pose_dist = sqrt((trajectory_pts(1,i)-trajectory_pts(1,i-1))^2 + (trajectory_pts(2,i)-trajectory_pts(2,i-1))^2 + (trajectory_pts(3,i)-trajectory_pts(3,i-1))^2);
            pose_dist = norm(trajectory_pts(:,i)-trajectory_pts(:,i-1));
            cumulative_arc_length(i,1) = cumulative_arc_length(i-1,1) + pose_dist;
            arcLength(i,1) = cumulative_arc_length(i,1);
        end
        
        prop_T = collect_tangent(:,record_err_idx);
        prop_N = collect_normal(:,record_err_idx);
        prop_B = collect_binormal(:,record_err_idx);
        
        % -- initial cuvature and torsion --
        fetch_idx = round((window_size+3)*0.5)+2;
        init_curvature = geometryParams.curvature(fetch_idx,1);
        init_torsion = geometryParams.torsion(fetch_idx,1);
 
        collect_initial_params(1,record_err_idx) = init_curvature;
        collect_initial_params(2,record_err_idx) = init_torsion;
        
        % -- initial Frenet frame --
        init_guess.T = collect_tangent(:,record_err_idx);
        init_guess.N = collect_normal(:,record_err_idx);
        init_guess.B = collect_binormal(:,record_err_idx);

        % -- GRADIENT DESCENT!!!! --
        if opt_FrenetFrame
            % -- Version II: optimize curvature, torsion, and normal vector --
            % -- i) convert initial frenet frame to unit sphere angles --
            [thetaT, phiT, theta_mkortho, phi_mkortho] = self_cvt_frenetframe_to_unit_sphere_angle(prop_T, prop_N);
            % -- ii) initial guesses --
            init_guess.curvature = 30;
            init_guess.torsion = -5;
            init_guess.theta1 = thetaT;
            init_guess.phi1 = phiT;
            init_guess.theta2 = theta_mkortho;
            init_guess.phi2 = phi_mkortho;
            % -- iii) gradient descent parameters --
            %GD_params.lr = 0.05;
            GD_params.lr = 0.05;
            GD_params.lr_frenetframe = 0.5;
            GD_params.iter_nums = 20000;
            GD_params.first_success_counts = 100;
            GD_params.second_success_counts = 300;
            % -- iv) do gradient descent --
            fprintf('Doing Gradient Descent ...\n');
            [updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_frenetframe_curve_fit(init_guess, GD_params, trajectory_pts, arcLength);
            % -- v) get output updates --
            updated_curvature = updated_geometry.k;
            updated_torsion = updated_geometry.tau;
            opt_thetaT = updated_geometry.theta1;
            opt_phiT = updated_geometry.phi1;
            opt_theta_mkortho = updated_geometry.theta2;
            opt_phi_mkortho = updated_geometry.phi2;
            % -- vi) convert unit sphere angles to frenet frame --
            [prop_T, prop_N, prop_B] = self_cvt_unit_sphere_angle_to_frenetframe(opt_thetaT, opt_phiT, opt_theta_mkortho, opt_phi_mkortho);
        else
            % -- Version I: optimize only curvature, torsion, and normal vector --
            init_guess.curvature = 20;
            init_guess.torsion = 2;

            GD_params.lr = 0.01;
            GD_params.iter_nums = 10000;
            GD_params.first_success_counts = 100;
            GD_params.second_success_counts = 300;

            [updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_curve_fit(init_guess, GD_params, trajectory_pts, arcLength);

            updated_curvature = updated_geometry.k;
            updated_torsion = updated_geometry.tau;
            prop_T = updated_geometry.T;
            prop_N = updated_geometry.N;
            prop_B = updated_geometry.B;
        end        
        
        if print_segment_wise && fr_idx == drawSegementStartFr
            fprintf("optimal curvature: %f\n", updated_curvature);
            fprintf("optimal torsion: %f\n", updated_torsion);
            fprintf("average point error: %f\n", collect_err(end,1));
            fprintf("error percentage: %f\n", err_pctge);
        end

        
        % -- draw 
        if fr_idx == drawSegementStartFr
            [curveFromModel, arcLength] = self_generateCurveFromModel(updated_curvature, updated_torsion, prop_T, prop_N, prop_B, trajectory_pts, curve_start_pt, 'peicewise_constant_model');
            
            drawTrajectoryPts = trajectory_pts;
            drawCurve_start_pt = curve_start_pt;
            drawCurveFromModel = curveFromModel;
            drawTangent = prop_T;
            drawNormal = prop_N;
            drawBinormal = prop_B;
            drawInitTangent = init_guess.T;
            drawInitNormal = init_guess.N;
            drawInitBinormal = init_guess.B;
            drawInitCurvature = geometryParams.curvature(3:end,1);
            drawInitTorsion = geometryParams.torsion(3:end,1);
            drawInitCurvature_prime = geometryParams.curvature_derivative;
            drawInitTorsion_prime = geometryParams.torsion_derivative;
            drawTime = t';
            drawUpdated_curvature = updated_curvature;
            drawUpdated_torsion = updated_torsion;
            drawCollectGDerrs = collect_err;
        end
        
        collect_all_init_traj_pts(:,record_err_idx) = trajectory_pts(:,2);
        %record_err_idx = record_err_idx + 1;

        if ~run_segment_wise
            collection_of_local_avg_errs(record_err_idx, 1) = collect_err(end,1);
            collect_err_percentage(record_err_idx, 1) = (collect_err(end,1) / arcLength(end,1))*100;
            collect_err_percentage(record_err_idx, 2) = arcLength(end,1);
        end
        record_err_idx = record_err_idx + 1;
        if mod(record_err_idx, 20)==0
            fprintf('. ');
        end
        if mod(record_err_idx, 500) == 0
            fprintf('\n');
        end
        
        
        if drawGeometryModelFit && fr_idx == drawSegementStartFr
            figure;
            self_drawGeometryModelFit(drawTime, drawInitCurvature, drawInitTorsion, drawInitCurvature_prime(3:end,1), drawInitTorsion_prime(3:end,1), abs(drawUpdated_curvature), drawUpdated_torsion);
            set(gcf,'color','w');
        end
        
        if drawCurveAndTrajectory && fr_idx == drawSegementStartFr
            %[init_curveFromModel, ~] = self_generateCurveFromModel(init_guess.curvature, init_guess.torsion, drawInitTangent, drawInitNormal, drawInitBinormal, drawTrajectoryPts, drawCurve_start_pt, 'peicewise_constant_model');
            [init_curveFromModel, ~] = self_generateCurveFromModel(drawInitCurvature(5,1), drawInitTorsion(5,1), drawInitTangent, drawInitNormal, drawInitBinormal, drawTrajectoryPts, drawCurve_start_pt, 'peicewise_constant_model');
            
            % -- plot trajectory --
            figure;
            plot3(drawTrajectoryPts(1,:), drawTrajectoryPts(2,:), drawTrajectoryPts(3,:), 'bo--', 'DisplayName', 'ground-truth trajectory');
            hold on;
            plot3(drawCurveFromModel(1,:), drawCurveFromModel(2,:), drawCurveFromModel(3,:), 'ro-', 'DisplayName', 'best-fit trajectory');
            %hold on;
            %plot3(init_curveFromModel(1,:), init_curveFromModel(2,:), init_curveFromModel(3,:), 'ko-', 'DisplayName', 'initial trajectory');
  
            % -- scale down --
            scale_factor = (arcLength(end,1)*0.2);
            drawTangent = drawTangent.*scale_factor;
            drawNormal = drawNormal.*scale_factor;
            drawBinormal = drawBinormal.*scale_factor;
            
            % -- draw frenet frame at the initial point --
            endpoint_T = drawTrajectoryPts(:,1) + drawTangent;
            endpoint_N = drawTrajectoryPts(:,1) + drawNormal;
            endpoint_B = drawTrajectoryPts(:,1) + drawBinormal;
            
            
            line([drawTrajectoryPts(1,1), endpoint_T(1,1)],[drawTrajectoryPts(2,1), endpoint_T(2,1)],[drawTrajectoryPts(3,1), endpoint_T(3,1)],'Color','g','LineStyle','-', 'DisplayName', 'tangent');
            hold on;
            line([drawTrajectoryPts(1,1), endpoint_N(1,1)],[drawTrajectoryPts(2,1), endpoint_N(2,1)],[drawTrajectoryPts(3,1), endpoint_N(3,1)],'Color','c','LineStyle','-', 'DisplayName', 'normal');
            hold on;
            line([drawTrajectoryPts(1,1), endpoint_B(1,1)],[drawTrajectoryPts(2,1), endpoint_B(2,1)],[drawTrajectoryPts(3,1), endpoint_B(3,1)],'Color','m','LineStyle','-', 'DisplayName', 'binormal');
            hold on;
            
            xlabel("x (m)");
            ylabel("y (m)");
            zlabel("z (m)");
            axis equal;
            legend;
            set(gcf,'color','w');
        end
    end

    if drawErrChangeAlongSequence
        
        figure;
        subplot(3,1,1)
        plot(start_fr:seq_end_fr, collect_err_percentage(:,2), 'b-');
        xlabel('frame index');
        ylabel({'total','arc length','(m)'});
        xlim([start_fr, seq_end_fr]);
        
        subplot(3,1,2)
        plot(start_fr:seq_end_fr, collection_of_local_avg_errs(:,1), 'b-');
        xlabel('frame index');
        ylabel({'average local','points error','(m)'});
        xlim([start_fr, seq_end_fr]);

        subplot(3,1,3)
        plot(start_fr:seq_end_fr, collect_err_percentage(:,1), 'b-');
        xlabel('frame index');
        ylabel({'error','percentage','(%)'});
        xlim([start_fr, seq_end_fr]);
        set(gcf,'color','w');
    end
end


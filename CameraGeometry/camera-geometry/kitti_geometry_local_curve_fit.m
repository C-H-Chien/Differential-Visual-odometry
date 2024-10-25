
% -- Verify that the camera center model C(t) is a sufficient model using KITTI dataset --
clc; clear all; close all;

% -- define parameters --
datasetName = 'KITTI';
sequenceName = '00';
verification_class = 'geometry';
category_left = '/image_0/'; % -- not used --

% -- read R, T, and time stamp of each frame --
[all_T, all_R, time_per_fr, ~] = readDataset(datasetName, sequenceName, category_left);

start_fr = 2130;
end_fr = size(all_T,2)-50;
seq_end_fr = size(all_T,2)-50;
%seq_end_fr = 2800;
window_size = 9;

if strcmp(verification_class, 'geometry')
    run_curve_model_segment_wise = 1;
end

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
drawGradientErrors = debugging & run_segment_wise;
drawErrorFunction = debugging & run_segment_wise;
drawFrenetFrameOnTrajectory = 1;
drawFrenetFrameFromVelAndAcc = 0;

piecewise_constant = 1;


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
    
    % -- define the initial Frenet-frame --
%     R_0 = all_R(:,:,2);
%     T_0_test = R_0(:,3);
%     N_0_test = R_0(:,1);
%     B_0_test = R_0(:,2);

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
    for fr_idx = start_fr:end_fr

        [dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(fr_idx, fr_idx+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);
        
        collect_vel_vector(:,record_err_idx) = dynamicParams.vel_vec(3,:)';
        collect_acc_vector(:,record_err_idx) = dynamicParams.acc_vec(3,:)';
        vec_acc_cross_prod = cross(collect_vel_vector(:,record_err_idx), collect_acc_vector(:,record_err_idx));
        collect_binormal(:,record_err_idx) = vec_acc_cross_prod / norm(vec_acc_cross_prod);
        collect_tangent(:,record_err_idx) = collect_vel_vector(:,record_err_idx) / norm(collect_vel_vector(:,record_err_idx));
        tangent_binormal_cross_prod = -cross(collect_tangent(:,record_err_idx), collect_binormal(:,record_err_idx));
        collect_normal(:,record_err_idx) = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);
        
        % -- arc lengths --
        trajectory_pts = r(:,3:end);
        curve_start_pt = r(:,3);
        
        % -- propogate the Frenet frame --
        % -- PROP 1) compute arc length from the last frenet frame --
        %pose_dist = sqrt((trajectory_pts(1,1)-r(1,1))^2 + (trajectory_pts(2,1)-r(2,1))^2 + (trajectory_pts(3,1)-r(3,1))^2);
        % -- PROP 1) propogate the frenet frame for one point forward --
        %[prop_T, prop_N, prop_B] = self_propogateFrenetFrame(geometryParams.curvature(1,1), geometryParams.torsion(1,1), lastT, lastN, lastB, pose_dist);
       
        prop_T = collect_tangent(:,record_err_idx);
        prop_N = collect_normal(:,record_err_idx);
        prop_B = collect_binormal(:,record_err_idx);
        
        % -- initial cuvature and torsion --
        fetch_idx = round((window_size+3)*0.5)+2;
        init_curvature = geometryParams.curvature(fetch_idx,1);
        init_torsion = geometryParams.torsion(fetch_idx,1);
        
        % -- do gradient descent to find optimal curvature and torsion --
        learning_rate = 0.05;
        gd_iterations = 5000;
        first_success_counts = 200;
        second_success_counts = 300;
        
        if fr_idx == 2129
            learning_rate = 10;
            gd_iterations = 20000;
            first_success_counts = 100;
            second_success_counts = 150;
        end

        % -- record errors over gradient descent iterations --
        if run_segment_wise
            if fr_idx == drawSegementStartFr
                collect_gd_errors = zeros(gd_iterations, 1);
            end
        else
            collect_gd_errors = zeros(gd_iterations, 1);
        end
        
        collect_initial_params(1,record_err_idx) = init_curvature;
        collect_initial_params(2,record_err_idx) = init_torsion;

%         % -- do gradient descent --
%         previous_curvature = init_curvature;
%         previous_torsion = init_torsion;
%         init_N = prop_N;
%         init_B = prop_B;
%         candidate_theta = 0;
%         for gd_iter = 1:gd_iterations
%             % -- generate a curve from model using propgated frenet frame
%             % and initial curvature and torsion --
%             [curveFromModel, arcLength] = self_generateCurveFromModel(previous_curvature, previous_torsion, prop_T, prop_N, prop_B, trajectory_pts, curve_start_pt, 'peicewise_constant_model');
% 
%             % -- get gradients of curvature and torsion --
%             [gradient_curvature, gradient_torsion, gradient_theta, avg_pts_err] = ...
%                 self_getGeometryParamsGradient(previous_curvature, previous_torsion, candidate_theta, arcLength, init_N, init_B, prop_T, prop_N, prop_B, trajectory_pts, curveFromModel);
% 
%             % -- update the curvature and torsion --
%             updated_curvature = previous_curvature + learning_rate * gradient_curvature;
%             updated_torsion = previous_torsion + learning_rate * gradient_torsion;
%             
%             % -- Frenet frame --
%             candidate_theta = candidate_theta + 0.05*gradient_theta;
%             
%             updated_candidate_N = cos(candidate_theta)*init_N + sin(candidate_theta)*init_B;
%             updated_binormal = cross(prop_T, updated_candidate_N);
%             prop_N = updated_candidate_N;
%             prop_B = updated_binormal;
% 
%             % -- prepare for the next iteration --
%             previous_curvature = updated_curvature;
%             previous_torsion = updated_torsion;
%             
%             % -- record errors during gradient descent iterations --
%             if run_segment_wise
%                 if fr_idx == drawSegementStartFr
%                     collect_gd_errors(gd_iter,1) = avg_pts_err;
%                 end
%             else
%                 collect_gd_errors(gd_iter,1) = avg_pts_err;
%             end
%             
%             % -- dynamically adjusting the learning rate --            
%             if gd_iter > 80
%                 if collect_gd_errors(gd_iter-1,1) < avg_pts_err
%                     collect_gd_errors(end,1) = collect_gd_errors(gd_iter-1,1);
%                     break;
%                 end
%             end
%         end
        
        % -- do gradient descent --
        previous_curvature = init_curvature;
        previous_torsion = init_torsion;
        minErr_curvature = init_curvature;
        minErr_torsion = init_torsion;
        init_N = prop_N;
        init_B = prop_B;
        candidate_theta = 0;
        fix_flag = 0;
        success_count = 0;
        max_avg_err = 1000;
        collect_learning_rate = zeros(gd_iterations, 1);
        for gd_iter = 1:gd_iterations
            % -- generate a curve from model using propgated frenet frame
            % and initial curvature and torsion --
            [curveFromModel, arcLength] = self_generateCurveFromModel(previous_curvature, previous_torsion, prop_T, prop_N, prop_B, trajectory_pts, curve_start_pt, 'peicewise_constant_model');

            % -- get gradients of curvature and torsion --
            [gradient_curvature, gradient_torsion, gradient_theta, avg_pts_err] = ...
                self_getGeometryParamsGradient(previous_curvature, previous_torsion, candidate_theta, arcLength, init_N, init_B, prop_T, prop_N, prop_B, trajectory_pts, curveFromModel);

            collect_learning_rate(gd_iter, 1) = learning_rate;
            
            err_percentage = (avg_pts_err / arcLength(end,1))*100;
            
            % -- if initially the error is already small enough, then the
            % learning rate should be small. This confines the gradient not to jump too far. --
            if (gd_iter == 1) && (avg_pts_err <= 0.1) && (err_percentage < 3)
                learning_rate = 0.00001;
                fix_flag = 1;
            else
                if (gd_iter > 100) && avg_pts_err > 1
                    learning_rate = learning_rate*0.1;
                    collect_gd_errors(gd_iter,1) = avg_pts_err;
                    previous_curvature = init_curvature;
                    previous_torsion = init_torsion;
                    continue;
                end
            end
            
            % -- update the curvature and torsion --
            updated_curvature = previous_curvature + learning_rate * gradient_curvature;
            updated_torsion = previous_torsion + learning_rate * gradient_torsion;
            
            if avg_pts_err < max_avg_err
                minErr_curvature = updated_curvature;
                minErr_torsion = updated_torsion;
                minErr_T = prop_T;
                minErr_N = prop_N;
                minErr_B = prop_B;
                max_avg_err = avg_pts_err;
            end
            
            % -- Frenet frame --
            candidate_theta = candidate_theta + 0.05*gradient_theta;
            
            updated_candidate_N = cos(candidate_theta)*init_N + sin(candidate_theta)*init_B;
            updated_binormal = cross(prop_T, updated_candidate_N);
            prop_N = updated_candidate_N;
            prop_B = updated_binormal;

            if (gd_iter > 1) && (fix_flag == 0)
                if avg_pts_err - collect_gd_errors(gd_iter-1,1) > 0.01 && (avg_pts_err < 0.28)
                    learning_rate = learning_rate*0.5;
                elseif avg_pts_err - collect_gd_errors(gd_iter-1,1) > 0.005 && (avg_pts_err < 0.28)
                    learning_rate = learning_rate*0.95;
                elseif (gd_iter > 1000) && (avg_pts_err - collect_gd_errors(gd_iter-1,1) < 0.0001) && (avg_pts_err < 0.1)
                    collect_gd_errors(end,1) = avg_pts_err;
                    break;
                elseif success_count >= first_success_counts
                    previous_curvature = updated_curvature;
                    previous_torsion = updated_torsion;
                    if avg_pts_err > 0.2
                        learning_rate = learning_rate*3;
                    else
                        learning_rate = learning_rate*1.6;
                    end
                    success_count = 0;
                elseif success_count >= second_success_counts
                    previous_curvature = updated_curvature;
                    previous_torsion = updated_torsion;
                    if avg_pts_err > 0.2
                        learning_rate = learning_rate*10;
                    else
                        learning_rate = learning_rate*8;
                    end
                    success_count = 0;
                else
                    % -- prepare for the next iteration --
                    previous_curvature = updated_curvature;
                    previous_torsion = updated_torsion;
                    success_count = success_count + 1;
                end
                collect_gd_errors(gd_iter,1) = avg_pts_err;
            else
                collect_gd_errors(gd_iter,1) = avg_pts_err;
            end
        end
        
        mask = collect_gd_errors(:,1) ~= 0;
        valid_collect_gd_errors = collect_gd_errors(mask);
        collect_gd_errors(end,1) = min(valid_collect_gd_errors(:,1));
        collect_curvature(1,record_err_idx) = updated_curvature;

        % -- print out the 
        if print_segment_wise && fr_idx == drawSegementStartFr
            fprintf("optimal curvature: %f\n", minErr_curvature);
            fprintf("optimal torsion: %f\n", minErr_torsion);
            fprintf("average point error: %f\n", collect_gd_errors(end,1));
        end
        
        % -- draw 
        if fr_idx == drawSegementStartFr
            [curveFromModel, arcLength] = self_generateCurveFromModel(minErr_curvature, minErr_torsion, minErr_T, minErr_N, minErr_B, trajectory_pts, curve_start_pt, 'peicewise_constant_model');
            
            drawTrajectoryPts = trajectory_pts;
            drawCurve_start_pt = curve_start_pt;
            drawCurveFromModel = curveFromModel;
            drawInitTangent = prop_T;
            drawInitNormal = prop_N;
            drawInitBinormal = prop_B;
            %drawTangent = updated_prop_T;
            %drawNormal = updated_prop_N;
            %drawBinormal = updated_prop_B;
            drawTangent = collect_tangent(:,record_err_idx);
            drawNormal = collect_normal(:,record_err_idx);
            drawBinormal = collect_binormal(:,record_err_idx);
            drawInitCurvature = geometryParams.curvature(3:end,1);
            drawInitTorsion = geometryParams.torsion(3:end,1);
            drawInitCurvature_prime = geometryParams.curvature_derivative;
            drawInitTorsion_prime = geometryParams.torsion_derivative;
            drawTime = t';
            drawUpdated_curvature = updated_curvature;
            drawUpdated_torsion = updated_torsion;
            drawCollectGDerrs = collect_gd_errors;
        end
        
        collect_all_init_traj_pts(:,record_err_idx) = trajectory_pts(:,2);
        %record_err_idx = record_err_idx + 1;

        if ~run_segment_wise
            collection_of_local_avg_errs(record_err_idx, 1) = collect_gd_errors(end,1);
            collect_err_percentage(record_err_idx, 1) = (collect_gd_errors(end,1) / arcLength(end,1))*100;
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
            [init_curveFromModel, ~] = self_generateCurveFromModel(drawInitCurvature(5,1), drawInitTorsion(5,1), drawInitTangent, drawInitNormal, drawInitBinormal, drawTrajectoryPts, drawCurve_start_pt, 'peicewise_constant_model');
            
            % -- plot trajectory --
            figure;
            plot3(drawTrajectoryPts(1,:), drawTrajectoryPts(2,:), drawTrajectoryPts(3,:), 'bo--', 'DisplayName', 'ground-truth trajectory');
            hold on;
            plot3(drawCurveFromModel(1,:), drawCurveFromModel(2,:), drawCurveFromModel(3,:), 'ro-', 'DisplayName', 'best-fit trajectory');
            hold on;
            plot3(init_curveFromModel(1,:), init_curveFromModel(2,:), init_curveFromModel(3,:), 'ko-', 'DisplayName', 'initial trajectory');
            
            % -- draw frenet frame --
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
        
        if drawGradientErrors
            figure;
            plot(1:gd_iterations, collect_gd_errors(:,1), 'b-');
            xlabel("iterations");
            ylabel("average error");
            set(gcf,'color','w');
        end
        
        if drawErrorFunction
            space_curvature = init_curvature / 100;
            space_torsion = init_torsion / 100;
            curvature_family = 0:space_curvature:init_curvature*2;
            torsion_family = 0:space_torsion:init_torsion*2;
            
            collection_avg_pts_err = zeros(size(curvature_family, 2), 1);
            for m = 1:size(curvature_family, 2)
                [curveMember, ~] = self_generateCurveFromModel(curvature_family(1,m), torsion_family(1,m), T_0, N_0, B_0, trajectory_pts, curve_start_pt, 'peicewise_constant_model');
                
                for i = 1:size(curveMember,2)
                    collection_avg_pts_err(m,1) = collection_avg_pts_err(m,1) + sqrt((trajectory_pts(1,i)-curveMember(1,i))^2 + (trajectory_pts(2,i)-curveMember(2,i))^2 + (trajectory_pts(3,i)-curveMember(3,i))^2);
                end
                collection_avg_pts_err(m,1) = collection_avg_pts_err(m,1) / size(curveMember,2);
            end
            
            figure;
            subplot(2,1,1);
            plot(curvature_family, collection_avg_pts_err(:,1), 'b-');
            xlabel('curvature');
            ylabel('error');
            
            subplot(2,1,2);
            plot(torsion_family, collection_avg_pts_err(:,1), 'b-');
            xlabel('torsion');
            ylabel('error');
            set(gcf,'color','w');
        end
    end
    
%     if drawFrenetFrameOnTrajectory
%         
%         % -- plot trajectory --
%         start_draw = 100;
%         end_draw = 120;
%         
%         %start_draw = start_fr;
%         %end_draw = end_fr-start_fr+1;
%         
%         figure;
%         plot3(collect_all_init_traj_pts(1,start_draw:end_draw), collect_all_init_traj_pts(2,start_draw:end_draw), collect_all_init_traj_pts(3,start_draw:end_draw), 'bo-', 'DisplayName', 'ground-truth trajectory');
%         hold on;
%         for i = start_draw:5:end_draw
%             endpoint_T = collect_all_init_traj_pts(:,i) + collect_propT(:,i);
%             endpoint_N = collect_all_init_traj_pts(:,i) + collect_propN(:,i);
%             endpoint_B = collect_all_init_traj_pts(:,i) + collect_propB(:,i);
%             line([collect_all_init_traj_pts(1,i), endpoint_T(1,1)],[collect_all_init_traj_pts(2,i), endpoint_T(2,1)],[collect_all_init_traj_pts(3,i), endpoint_T(3,1)],'Color','m','LineStyle','-', 'DisplayName', 'tangent');
%             hold on;
%             line([collect_all_init_traj_pts(1,i), endpoint_N(1,1)],[collect_all_init_traj_pts(2,i), endpoint_N(2,1)],[collect_all_init_traj_pts(3,i), endpoint_N(3,1)],'Color','c','LineStyle','-', 'DisplayName', 'normal');
%             hold on;
%             line([collect_all_init_traj_pts(1,i), endpoint_B(1,1)],[collect_all_init_traj_pts(2,i), endpoint_B(2,1)],[collect_all_init_traj_pts(3,i), endpoint_B(3,1)],'Color','k','LineStyle','-', 'DisplayName', 'binormal');
%             hold on;
%             axis equal;
%         end
%         xlabel("x (m)");
%         ylabel("y (m)");
%         zlabel("z (m)");
%         axis equal;
%         %legend;
%         set(gcf,'color','w');
%     end
    
    if drawFrenetFrameFromVelAndAcc
        
        % -- plot trajectory --
        start_draw = 10;
        end_draw = 30;
        
        figure;
        plot3(collect_all_init_traj_pts(1,start_draw:end_draw), collect_all_init_traj_pts(2,start_draw:end_draw), collect_all_init_traj_pts(3,start_draw:end_draw), 'bo-', 'DisplayName', 'ground-truth trajectory');
        hold on;
        for i = start_draw:5:end_draw
            endpoint_T = collect_all_init_traj_pts(:,i) + collect_tangent(:,i);
            endpoint_N = collect_all_init_traj_pts(:,i) + collect_normal(:,i);
            endpoint_B = collect_all_init_traj_pts(:,i) + collect_binormal(:,i);
            line([collect_all_init_traj_pts(1,i), endpoint_T(1,1)],[collect_all_init_traj_pts(2,i), endpoint_T(2,1)],[collect_all_init_traj_pts(3,i), endpoint_T(3,1)],'Color','m','LineStyle','-', 'DisplayName', 'tangent');
            hold on;
            line([collect_all_init_traj_pts(1,i), endpoint_N(1,1)],[collect_all_init_traj_pts(2,i), endpoint_N(2,1)],[collect_all_init_traj_pts(3,i), endpoint_N(3,1)],'Color','c','LineStyle','-', 'DisplayName', 'normal');
            hold on;
            line([collect_all_init_traj_pts(1,i), endpoint_B(1,1)],[collect_all_init_traj_pts(2,i), endpoint_B(2,1)],[collect_all_init_traj_pts(3,i), endpoint_B(3,1)],'Color','k','LineStyle','-', 'DisplayName', 'binormal');
            hold on;
            axis equal;
        end
        xlabel("x (m)");
        ylabel("y (m)");
        zlabel("z (m)");
        axis equal;
        %legend;
        set(gcf,'color','w');
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


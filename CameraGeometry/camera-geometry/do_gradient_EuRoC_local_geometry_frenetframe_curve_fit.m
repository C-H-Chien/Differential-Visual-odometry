function [updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_frenetframe_curve_fit(init_guess, GD_params, trajectory_pts, arcLength)
% Gradient Descent Used to Fit a Local Segment of Curve 
% GradientVariables: curvature, torsion, theta1, phi1, theta2, phi2

% gradient descent parameters
learning_rate = GD_params.lr;
learning_rate_ff = GD_params.lr_frenetframe;
gd_iterations = GD_params.iter_nums;
first_success_counts = GD_params.first_success_counts;
second_success_counts = GD_params.second_success_counts;

% initial guesses 
previous_curvature = init_guess.curvature;
previous_torsion = init_guess.torsion;
previous_theta_T = init_guess.theta1;
previous_phi_T = init_guess.phi1;
previous_theta_mkortho = init_guess.theta2;
previous_phi_mkortho = init_guess.phi2;

success_count = 0;
max_avg_err = 1000;
fix_flag = 0;

% -- either "avg-pt-err" or "end-pt-err" --
energy_func_type = "avg-pt-err";

collect_lr = zeros(gd_iterations, 1);
collect_err = zeros(gd_iterations, 1);
for gd_iter = 1:gd_iterations

    if mod(gd_iter, 50) == 0
        fprintf('. ');
    end
    if mod(gd_iter, 500) == 0
        fprintf('\n');
    end
    
        
    [T, N, B] = self_cvt_unit_sphere_angle_to_frenetframe(previous_theta_T, previous_phi_T, previous_theta_mkortho, previous_phi_mkortho);

    % -- generate a helix from model --
    [curveFromModel, ~] = self_generateHelixFromModel(previous_curvature, previous_torsion, T, N, B, arcLength', trajectory_pts(:,1));
    % -- get gradients --
    % -- the energy function accounts only the end-point error --
    gradient_curvature = 0;
    gradient_torsion = 0;
    gradient_thetaT = 0;
    gradient_phiT = 0;
    gradient_thetaMkOrtho = 0;
    gradient_phiMkOrtho = 0;
    
    if strcmp(energy_func_type, "end-pt-err")
        % -- get derivatives w.r.t. variables (curvature, torsion, theta1, phi1, theta2, phi2) --
        Df = self_diff_func_wrt_curve_geometry_frenetframe(previous_curvature, previous_torsion, previous_theta_T, previous_phi_T, previous_theta_mkortho, previous_phi_mkortho, arcLength(end,1));

        end_pt_err_x = trajectory_pts(1,end)-curveFromModel(1,end);
        end_pt_err_y = trajectory_pts(2,end)-curveFromModel(2,end);
        end_pt_err_z = trajectory_pts(3,end)-curveFromModel(3,end);

        % -- x dimension --
        gradient_thetaT = gradient_thetaT + end_pt_err_x*Df.f1theta1;
        gradient_phiT = gradient_phiT + end_pt_err_x*Df.f1phi1;
        gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_x*Df.f1theta2;
        gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_x*Df.f1phi2; 
        gradient_curvature = gradient_curvature + end_pt_err_x*Df.f1k;
        gradient_torsion = gradient_torsion + end_pt_err_x*Df.f1tau;

        % -- y dimension --
        gradient_curvature = gradient_curvature + end_pt_err_y*Df.f2k;
        gradient_torsion = gradient_torsion + end_pt_err_y*Df.f2tau;
        gradient_thetaT = gradient_thetaT + end_pt_err_y*Df.f2theta1;
        gradient_phiT = gradient_phiT + end_pt_err_y*Df.f2phi1;
        gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_y*Df.f2theta2;
        gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_y*Df.f2phi2;

        % -- z dimension --
        gradient_curvature = gradient_curvature + end_pt_err_z*Df.f3k;
        gradient_torsion = gradient_torsion + end_pt_err_z*Df.f3tau;
        gradient_thetaT = gradient_thetaT + end_pt_err_z*Df.f3theta1;
        gradient_phiT = gradient_phiT + end_pt_err_z*Df.f3phi1;
        gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_z*Df.f3theta2;
        gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_z*Df.f3phi2;
        
        % -- return the end-point trajectory-curve point-wise error and error percentage --
        avg_pts_err = 0;
        %for i = 1:size(trajectory_pts,2)
            avg_pts_err = avg_pts_err + norm(trajectory_pts(:,end)-curveFromModel(:,end));
        %end
        %avg_pts_err = avg_pts_err / size(trajectory_pts,2);
    elseif strcmp(energy_func_type, "avg-pt-err")
        for d = 2:size(arcLength,1)
            % -- get derivatives w.r.t. variables (curvature, torsion, theta1, phi1, theta2, phi2) --
            Df = self_diff_func_wrt_curve_geometry_frenetframe(previous_curvature, previous_torsion, previous_theta_T, previous_phi_T, previous_theta_mkortho, previous_phi_mkortho, arcLength(d,1));
            end_pt_err_x = trajectory_pts(1,d)-curveFromModel(1,d);
            end_pt_err_y = trajectory_pts(2,d)-curveFromModel(2,d);
            end_pt_err_z = trajectory_pts(3,d)-curveFromModel(3,d);
            
            % -- x dimension --
            gradient_thetaT = gradient_thetaT + end_pt_err_x*Df.f1theta1;
            gradient_phiT = gradient_phiT + end_pt_err_x*Df.f1phi1;
            gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_x*Df.f1theta2;
            gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_x*Df.f1phi2; 
            gradient_curvature = gradient_curvature + end_pt_err_x*Df.f1k;
            gradient_torsion = gradient_torsion + end_pt_err_x*Df.f1tau;

            % -- y dimension --
            gradient_curvature = gradient_curvature + end_pt_err_y*Df.f2k;
            gradient_torsion = gradient_torsion + end_pt_err_y*Df.f2tau;
            gradient_thetaT = gradient_thetaT + end_pt_err_y*Df.f2theta1;
            gradient_phiT = gradient_phiT + end_pt_err_y*Df.f2phi1;
            gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_y*Df.f2theta2;
            gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_y*Df.f2phi2;

            % -- z dimension --
            gradient_curvature = gradient_curvature + end_pt_err_z*Df.f3k;
            gradient_torsion = gradient_torsion + end_pt_err_z*Df.f3tau;
            gradient_thetaT = gradient_thetaT + end_pt_err_z*Df.f3theta1;
            gradient_phiT = gradient_phiT + end_pt_err_z*Df.f3phi1;
            gradient_thetaMkOrtho = gradient_thetaMkOrtho + end_pt_err_z*Df.f3theta2;
            gradient_phiMkOrtho = gradient_phiMkOrtho + end_pt_err_z*Df.f3phi2;
        end
        
        % -- return the average trajectory-curve point-wise error and error percentage --
        avg_pts_err = 0;
        for i = 1:size(trajectory_pts,2)
            avg_pts_err = avg_pts_err + norm(trajectory_pts(:,i)-curveFromModel(:,i));
        end
        avg_pts_err = avg_pts_err / size(trajectory_pts,2);
    end
    
    
    err_percentage = (avg_pts_err / arcLength(end,1))*100;

    collect_lr(gd_iter, 1) = learning_rate;
  

    % -- if initially the error is already small enough, then the
    % learning rate should be small. This confines the gradient not to jump too far. --
    if (gd_iter == 1) && (avg_pts_err <= 0.1) && (err_percentage < 3)
        learning_rate = 0.00001;
        fix_flag = 1;
    else
        if (gd_iter > 100) && avg_pts_err > 1
            learning_rate = learning_rate*0.1;
            collect_err(gd_iter,1) = avg_pts_err;
            previous_curvature = init_curvature;
            previous_torsion = init_torsion;
            continue;
        end
    end

    % -- update the curvature and torsion --
    updated_curvature = previous_curvature + learning_rate * gradient_curvature;
    updated_torsion = previous_torsion + learning_rate * gradient_torsion;
    updated_thetaT = previous_theta_T + learning_rate_ff * gradient_thetaT;
    updated_phiT = previous_phi_T + learning_rate_ff * gradient_phiT;
    updated_theta_mkortho = previous_theta_mkortho + learning_rate_ff * gradient_thetaMkOrtho;
    updated_phi_mkortho = previous_phi_mkortho + learning_rate_ff * gradient_phiMkOrtho;
    

    if avg_pts_err < max_avg_err
        minErr_curvature = updated_curvature;
        minErr_torsion = updated_torsion;
        minErr_thetaT = updated_thetaT;
        minErr_phiT = updated_phiT;
        minErr_theta_mkortho = updated_theta_mkortho;
        minErr_phi_mkortho = updated_phi_mkortho;
        max_avg_err = avg_pts_err;
    end

    if (gd_iter > 1) && (fix_flag == 0)
        if avg_pts_err - collect_err(gd_iter-1,1) > 0.002 && (avg_pts_err < 0.02)
            learning_rate = learning_rate*0.5;
            learning_rate_ff = learning_rate_ff*0.5;
        elseif avg_pts_err - collect_err(gd_iter-1,1) > 0.0005 && (avg_pts_err < 0.03)
            learning_rate = learning_rate*0.95;
            learning_rate_ff = learning_rate_ff * 0.95;
        elseif (gd_iter > 1000) && (avg_pts_err - collect_err(gd_iter-1,1) < 0.00001) && (avg_pts_err < 0.0001)
            collect_err(end,1) = avg_pts_err;
            break;
        elseif (gd_iter > 1000) && (avg_pts_err - collect_err(gd_iter-1,1) > 0) && (err_percentage > 1)
            learning_rate = learning_rate*0.5;
        elseif success_count >= first_success_counts
            previous_curvature = updated_curvature;
            previous_torsion = updated_torsion;
            previous_theta_T = updated_thetaT;
            previous_phi_T = updated_phiT;
            previous_theta_mkortho = updated_theta_mkortho;
            previous_phi_mkortho = updated_phi_mkortho;
            if avg_pts_err > 0.2
                learning_rate = learning_rate*3;
            else
                learning_rate = learning_rate*1.6;
            end
            learning_rate_ff = learning_rate_ff * 1.1;
            success_count = 0;
        elseif success_count >= second_success_counts
            previous_curvature = updated_curvature;
            previous_torsion = updated_torsion;
            previous_theta_T = updated_thetaT;
            previous_phi_T = updated_phiT;
            previous_theta_mkortho = updated_theta_mkortho;
            previous_phi_mkortho = updated_phi_mkortho;
            if avg_pts_err > 0.2
                learning_rate = learning_rate*10;
            else
                learning_rate = learning_rate*8;
            end
            learning_rate_ff = learning_rate_ff * 1.1;
            success_count = 0;
        else
            % -- prepare for the next iteration --
            previous_curvature = updated_curvature;
            previous_torsion = updated_torsion;
            previous_theta_T = updated_thetaT;
            previous_phi_T = updated_phiT;
            previous_theta_mkortho = updated_theta_mkortho;
            previous_phi_mkortho = updated_phi_mkortho;
            success_count = success_count + 1;
        end
        collect_err(gd_iter,1) = avg_pts_err;
    else
        collect_err(gd_iter,1) = avg_pts_err;
    end
end

mask = collect_err(:,1) ~= 0;
valid_collect_gd_errors = collect_err(mask);
collect_err(end,1) = min(valid_collect_gd_errors(:,1));

% outputs
updated_geometry.k = minErr_curvature;
updated_geometry.tau = minErr_torsion;
updated_geometry.theta1 = minErr_thetaT;
updated_geometry.phi1 = minErr_phiT;
updated_geometry.theta2 = minErr_theta_mkortho;
updated_geometry.phi2 = minErr_phi_mkortho;
err_pctge = (collect_err(end,1) / arcLength(end,1))*100;

if strcmp(energy_func_type, "end-pt-err")
    [T, N, B] = self_cvt_unit_sphere_angle_to_frenetframe(updated_geometry.theta1, updated_geometry.phi1, updated_geometry.theta2, updated_geometry.phi2);

    % -- generate a helix from model --
    [curveFromModel, ~] = self_generateHelixFromModel(updated_geometry.k, updated_geometry.tau, T, N, B, arcLength', trajectory_pts(:,1));
    
    % -- compute average error --
    avg_pts_err = 0;
    for i = 1:size(trajectory_pts,2)
        avg_pts_err = avg_pts_err + norm(trajectory_pts(:,i)-curveFromModel(:,i));
    end
    avg_pts_err = avg_pts_err / size(trajectory_pts,2);
    
    collect_err(end,1) = avg_pts_err;
    err_pctge = (collect_err(end,1) / arcLength(end,1))*100;
end
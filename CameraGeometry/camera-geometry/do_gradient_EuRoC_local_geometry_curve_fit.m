function [updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_curve_fit(init_guess, GD_params, trajectory_pts, arcLength)
% Gradient Descent Used to Fit a Local Segment of Curve 
% Energy Functions: either average point-wise error or end-point error
% 
%

% gradient descent parameters
learning_rate = GD_params.lr;
gd_iterations = GD_params.iter_nums;
first_success_counts = GD_params.first_success_counts;
second_success_counts = GD_params.second_success_counts;

% initial guesses 
previous_curvature = init_guess.curvature;
previous_torsion = init_guess.torsion;
prop_T = init_guess.T;
prop_N = init_guess.N;
prop_B = init_guess.B;
init_N = init_guess.N;
init_B = init_guess.B;

success_count = 0;
max_avg_err = 1000;
fix_flag = 0;
candidate_theta = 0;

energy_func_type = "avg-pt-err";
%energy_func_type = "end-pt-err";

collect_lr = zeros(gd_iterations, 1);
collect_err = zeros(gd_iterations, 1);
for gd_iter = 1:gd_iterations
    % -- generate a curve from model using propgated frenet frame
    % and initial curvature and torsion --
    %[curveFromModel, arcLength] = self_generateCurveFromModel(previous_curvature, previous_torsion, prop_T, prop_N, prop_B, trajectory_pts, trajectory_pts(:,1), 'peicewise_constant_model');
    [curveFromModel, ~] = self_generateHelixFromModel(previous_curvature, previous_torsion, prop_T, prop_N, prop_B, arcLength', trajectory_pts(:,1));

    if strcmp(energy_func_type, "avg-pt-err")
        % -- get gradients of curvature and torsion --
        [gradient_curvature, gradient_torsion, gradient_theta, avg_pts_err] = ...
            self_getGeometryParamsGradient(previous_curvature, previous_torsion, candidate_theta, arcLength, init_N, init_B, prop_T, prop_N, prop_B, trajectory_pts, curveFromModel);
    elseif strcmp(energy_func_type, "end-pt-err")
        % -- the energy function accounts only the end-point error --
        Df = self_diff_func_wrt_curve_geometry(previous_curvature, previous_torsion, arcLength(end,1));
        I = eye(3);
        diff_energy_func_k = 0;
        diff_energy_func_tau = 0;
        
        k0 = previous_curvature;
        tau0 = previous_torsion;
        omega0 = sqrt(k0^2 + tau0^2);
        coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(end,1)));
        dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
        for d = 1:3
            Derivative_k = Df.f1k*(I(d,:)*prop_T) + Df.f2k*(I(d,:)*prop_N) + Df.f3k*(I(d,:)*prop_B);
            Derivative_tau = Df.f1t*(I(d,:)*prop_T) + Df.f2t*(I(d,:)*prop_N) + Df.f3t*(I(d,:)*prop_B);
            diff_energy_func_k = diff_energy_func_k + (trajectory_pts(d,end)-curveFromModel(d,end))*Derivative_k;
            diff_energy_func_tau = diff_energy_func_tau + (trajectory_pts(d,end)-curveFromModel(d,end))*Derivative_tau;
        end
        gradient_curvature = diff_energy_func_k;
        gradient_torsion = diff_energy_func_tau;
        
        gradient_theta = (trajectory_pts(1,end)-curveFromModel(1,end))*dCdtheta(1,1) + (trajectory_pts(2,end)-curveFromModel(2,end))*dCdtheta(2,1) + (trajectory_pts(3,end)-curveFromModel(3,end))*dCdtheta(3,1);
        
        % -- return the average trajectory-curve point-wise error --
        avg_pts_err = 0;
%         for i = 1:size(trajectory_pts,2)
%             avg_pts_err = avg_pts_err + norm(trajectory_pts(:,i)-curveFromModel(:,i)); %sqrt((C_GT(1,i)-C0(1,i))^2 + (C_GT(2,i)-C0(2,i))^2 + (C_GT(3,i)-C0(3,i))^2);
%         end
%         avg_pts_err = avg_pts_err / size(trajectory_pts,2);
        avg_pts_err = avg_pts_err + norm(trajectory_pts(:,end)-curveFromModel(:,end));
    end

    collect_lr(gd_iter, 1) = learning_rate;

    err_percentage = (avg_pts_err / arcLength(end,1))*100;

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

    %collect_gradients(gd_iter, 1) = gradient_curvature;
    %collect_gradients(gd_iter, 2) = gradient_torsion;

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
    if updated_curvature > 0
        candidate_theta = candidate_theta + 0.5*gradient_theta;
    else
        candidate_theta = candidate_theta + pi;
    end

    updated_candidate_N = cos(candidate_theta)*init_N + sin(candidate_theta)*init_B;
    updated_binormal = cross(prop_T, updated_candidate_N);
    prop_N = updated_candidate_N;
    prop_B = updated_binormal;

    if (gd_iter > 1) && (fix_flag == 0)
        if avg_pts_err - collect_err(gd_iter-1,1) > 0.002 && (avg_pts_err < 0.02)
            learning_rate = learning_rate*0.5;
        elseif avg_pts_err - collect_err(gd_iter-1,1) > 0.0005 && (avg_pts_err < 0.03)
            learning_rate = learning_rate*0.95;
        elseif (gd_iter > 1000) && (avg_pts_err - collect_err(gd_iter-1,1) < 0.00001) && (avg_pts_err < 0.0001)
            collect_err(end,1) = avg_pts_err;
            break;
        elseif (gd_iter > 1000) && (avg_pts_err - collect_err(gd_iter-1,1) > 0) && (err_percentage > 1)
            learning_rate = learning_rate*0.5;
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
updated_geometry.T = minErr_T;
updated_geometry.N = minErr_N;
updated_geometry.B = minErr_B;
err_pctge = (collect_err(end,1) / arcLength(end,1))*100;

if strcmp(energy_func_type, "end-pt-err")
    % -- generate a helix from model --
    [curveFromModel, ~] = self_generateHelixFromModel(updated_geometry.k, updated_geometry.tau, updated_geometry.T, updated_geometry.N, updated_geometry.B, arcLength', trajectory_pts(:,1));
    
    % -- compute average error --
    avg_pts_err = 0;
    for i = 1:size(trajectory_pts,2)
        avg_pts_err = avg_pts_err + norm(trajectory_pts(:,i)-curveFromModel(:,i));
    end
    avg_pts_err = avg_pts_err / size(trajectory_pts,2);
    
    collect_err(end,1) = avg_pts_err;
    err_pctge = (collect_err(end,1) / arcLength(end,1))*100;
end

end
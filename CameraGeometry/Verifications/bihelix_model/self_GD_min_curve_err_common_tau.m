function [opt_res, collect_gd_errors, collect_learning_rate, minErr_err_percentage] = self_GD_min_curve_err_common_tau(gdParams, initParams, initFrenetFrame, trajectory_pts, j_idx, err_type)
    % -- do gradient descent --
    learning_rate = gdParams.lr;
    gd_iterations = gdParams.num_of_iter;
    first_success_counts = gdParams.num_of_first_success;
    second_success_counts = gdParams.num_of_second_success;
    collect_gd_errors = zeros(gd_iterations, 1);
    
    k0 = initParams.init_k0;
    tau0 = initParams.init_tau0;
    k1 = initParams.init_k1;
    
    minErr_k0 = k0;
    minErr_tau0 = tau0;
    minErr_k1 = k1;
    
    prop_T = initFrenetFrame.T;
    prop_N = initFrenetFrame.N;
    prop_B = initFrenetFrame.B;
    
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
        [bihelixCurveFromModel, arcLength] = self_generateBihelixFromModel(k0, tau0, k1, tau0, prop_T, prop_N, prop_B, trajectory_pts, j_idx);
        
        % -- get gradients of curvature and torsion --
        [gradient_k0, gradient_tau0, gradient_k1, gradient_theta, energy_val] = ...
        self_getGeometryParamsGradientForBihelix_common_tau(k0, tau0, k1, candidate_theta, arcLength, init_N, init_B, prop_T, prop_N, prop_B, trajectory_pts, bihelixCurveFromModel, j_idx, err_type);
        
        % -- record all learning rates throughout gradient descent --
        collect_learning_rate(gd_iter, 1) = learning_rate;
        
        % -- calculate error percentage --
        err_percentage = (energy_val / arcLength(end,1))*100;

        % -- if initially the error is already small enough, then the
        % learning rate should be small. This confines the gradient not to jump too far. --
        if (gd_iter == 1) && (energy_val <= 0.01) && (err_percentage < 3)
            learning_rate = 0.00001;
            fix_flag = 1;
        else
            if (gd_iter > 100) && energy_val > 1
                learning_rate = learning_rate*0.1;
                collect_gd_errors(gd_iter,1) = energy_val;
                k0 = initParams.init_k0;
                tau0 = initParams.init_tau0;
                k1 = initParams.init_k1;
                continue;
            end
        end

        if energy_val < max_avg_err
            minErr_k0 = k0;
            minErr_tau0 = tau0;
            minErr_k1 = k1;
            minErr_T = prop_T;
            minErr_N = prop_N;
            minErr_B = prop_B;
            minErr_err_percentage = err_percentage;
            max_avg_err = energy_val;
        end
        
        % -- update the curvature and torsion --
        updated_k0 = k0 + learning_rate * gradient_k0;
        updated_tau0 = tau0 + learning_rate * gradient_tau0;
        updated_k1 = k1 + learning_rate * gradient_k1;

        % -- Frenet frame --
        %candidate_theta = candidate_theta + 0.05*gradient_theta;
        candidate_theta = candidate_theta + 0.9*gradient_theta;

        updated_candidate_N = cos(candidate_theta)*init_N + sin(candidate_theta)*init_B;
        updated_binormal = cross(prop_T, updated_candidate_N);
        prop_N = updated_candidate_N;
        prop_B = updated_binormal;

        if (gd_iter > 1) && (fix_flag == 0)
            if energy_val > collect_gd_errors(gd_iter-1,1) && (energy_val < 0.01)
                learning_rate = learning_rate*0.5;
            elseif (gd_iter > 1000) && (energy_val > collect_gd_errors(gd_iter-1,1)) && (err_percentage < 3)
                collect_gd_errors(end,1) = energy_val;
                break;
            elseif energy_val > collect_gd_errors(gd_iter-1,1)
                learning_rate = learning_rate*0.5;
            elseif success_count >= first_success_counts
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                if energy_val > 0.1
                    learning_rate = learning_rate*2;
                else
                    learning_rate = learning_rate*1.2;
                end
                success_count = 0;
            elseif success_count >= second_success_counts
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                if energy_val > 0.1
                    learning_rate = learning_rate*2;
                else
                    learning_rate = learning_rate*1.5;
                end
                success_count = 0;
            else
                % -- prepare for the next iteration --
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                success_count = success_count + 1;
            end
            collect_gd_errors(gd_iter,1) = energy_val;
        else
            collect_gd_errors(gd_iter,1) = energy_val;
        end
    end

    mask = collect_gd_errors(:,1) ~= 0;
    valid_collect_gd_errors = collect_gd_errors(mask);
    collect_gd_errors(end,1) = min(valid_collect_gd_errors(:,1));
    
    % -- output the optimal results --
    opt_res.k0 = minErr_k0;
    opt_res.tau0 = minErr_tau0;
    opt_res.k1 = minErr_k1;
    opt_res.T = minErr_T;
    opt_res.N = minErr_N;
    opt_res.B = minErr_B;
    opt_res.energy_val = max_avg_err;
    
end
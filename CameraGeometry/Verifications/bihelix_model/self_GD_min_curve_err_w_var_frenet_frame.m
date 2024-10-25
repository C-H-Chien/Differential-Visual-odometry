function [opt_res, collect_gd_errors, collect_learning_rate, minErr_err_percentage] = ...
         self_GD_min_curve_err_w_var_frenet_frame(gdParams, initParams, trajectory_pts, j_idx, err_type)
    % -- do gradient descent --
    learning_rate = gdParams.lr;
    learning_rate_ff = gdParams.lr_frenetframe;
    gd_iterations = gdParams.num_of_iter;
    first_success_counts = gdParams.num_of_first_success;
    second_success_counts = gdParams.num_of_second_success;
    collect_gd_errors = zeros(gd_iterations, 1);
    
    k0 = initParams.init_k0;
    tau0 = initParams.init_tau0;
    k1 = initParams.init_k1;
    tau1 = initParams.init_tau1;
    
%     R = [initFrenetFrame.T'; initFrenetFrame.N'; initFrenetFrame.B'];
%     eul_angle = rotm2eul(R);
%     theta3 = eul_angle(1);
%     theta2 = eul_angle(2);
%     theta1 = eul_angle(3);

    theta1 = initParams.theta1;
    phi1 = initParams.phi1;
    theta2 = initParams.theta2;
    phi2 = initParams.phi2;

    minErr_k0 = k0;
    minErr_tau0 = tau0;
    minErr_k1 = k1;
    minErr_tau1 = tau1;
    
    fix_flag = 0;
    success_count = 0;
    max_avg_err = 1000;
    collect_learning_rate = zeros(gd_iterations, 1);
    for gd_iter = 1:gd_iterations
        % -- construct Frenet frame from three angles --
%         FF = eul2rotm([theta3, theta2, theta1]);
%         T0 = FF(1,:)';
%         N0 = FF(2,:)';
%         B0 = FF(3,:)';
        [T0, N0, B0] = self_cvt_unit_sphere_angle_to_frenetframe(theta1, phi1, theta2, phi2);

        % -- generate a curve from model using propgated frenet frame
        % and initial curvature and torsion --
        [bihelixCurveFromModel, arcLength] = self_generateBihelixFromModel(k0, tau0, k1, tau1, T0, N0, B0, trajectory_pts, j_idx);
        
        % -- get gradients of curvature and torsion --
        [gradient_k0, gradient_tau0, gradient_k1, gradient_tau1, gradient_theta1, gradient_phi1, gradient_theta2, gradient_phi2, energy_val] = ...
        self_getGeometryParamsFrenetFrameGradientsForBihelix(k0, tau0, k1, tau1, theta1, phi1, theta2, phi2, arcLength, trajectory_pts, bihelixCurveFromModel, j_idx, err_type);

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
                tau1 = initParams.init_tau1;
                theta1 = initParams.theta1;
                theta2 = initParams.theta2;
                phi1 = initParams.phi1;
                phi2 = initParams.phi2;
                continue;
            end
        end

        if energy_val < max_avg_err
            minErr_k0 = k0;
            minErr_tau0 = tau0;
            minErr_k1 = k1;
            minErr_tau1 = tau1;
            %minErr_FF = eul2rotm([theta3, theta2, theta1]);
            minErr_theta1 = theta1;
            minErr_theta2 = theta2;
            minErr_phi1 = phi1;
            minErr_phi2 = phi2;
            minErr_err_percentage = err_percentage;
            max_avg_err = energy_val;
        end
        
        % -- update the curvature and torsion --
        updated_k0 = k0 + learning_rate * gradient_k0;
        updated_tau0 = tau0 + learning_rate * gradient_tau0;
        updated_k1 = k1 + learning_rate * gradient_k1;
        updated_tau1 = tau1 + learning_rate * gradient_tau1;
        updated_theta1 = theta1 + learning_rate_ff * gradient_theta1;
        updated_theta2 = theta2 + learning_rate_ff * gradient_theta2;
        updated_phi1 = phi1 + learning_rate_ff * gradient_phi1;
        updated_phi2 = phi2 + learning_rate_ff * gradient_phi2;

        if (gd_iter > 1) && (fix_flag == 0)
            if energy_val > collect_gd_errors(gd_iter-1,1) && (energy_val < 0.01)
                learning_rate = learning_rate*0.5;
                learning_rate_ff = learning_rate_ff*0.5;
            elseif (gd_iter > 1000) && (energy_val > collect_gd_errors(gd_iter-1,1)) && (err_percentage < 3)
                collect_gd_errors(end,1) = energy_val;
                break;
            elseif energy_val > collect_gd_errors(gd_iter-1,1)
                learning_rate = learning_rate*0.5;
                learning_rate_ff = learning_rate_ff*0.5;
            elseif success_count >= first_success_counts
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                tau1 = updated_tau1;
                theta1 = updated_theta1;
                theta2 = updated_theta2;
                phi1 = updated_phi1;
                phi2 = updated_phi2;
                if energy_val > 0.1
                    learning_rate = learning_rate*2;
                    learning_rate_ff = learning_rate_ff*2;
                else
                    learning_rate = learning_rate*1.2;
                    learning_rate_ff = learning_rate_ff * 1.2;
                end
                success_count = 0;
            elseif success_count >= second_success_counts
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                tau1 = updated_tau1;
                theta1 = updated_theta1;
                theta2 = updated_theta2;
                phi1 = updated_phi1;
                phi2 = updated_phi2;
                if energy_val > 0.1
                    learning_rate = learning_rate*2;
                    learning_rate_ff = learning_rate_ff*2;
                else
                    learning_rate = learning_rate*1.5;
                    learning_rate_ff = learning_rate_ff*1.5;
                end
                success_count = 0;
            else
                % -- prepare for the next iteration --
                k0 = updated_k0;
                tau0 = updated_tau0;
                k1 = updated_k1;
                tau1 = updated_tau1;
                theta1 = updated_theta1;
                theta2 = updated_theta2;
                phi1 = updated_phi1;
                phi2 = updated_phi2;
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
    opt_res.tau1 = minErr_tau1;
    opt_res.theta1 = minErr_theta1;
    opt_res.theta2 = minErr_theta2;
    opt_res.phi1 = minErr_phi1;
    opt_res.phi2 = minErr_phi2;
    opt_res.energy_val = max_avg_err;
    
end
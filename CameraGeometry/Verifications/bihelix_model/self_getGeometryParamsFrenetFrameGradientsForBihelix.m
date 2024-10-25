function [grad_k0, grad_tau0, grad_k1, grad_tau1, grad_theta1, grad_phi1, grad_theta2, grad_phi2, energy_val] = ...
    self_getGeometryParamsFrenetFrameGradientsForBihelix(k0, tau0, k1, tau1, theta1, phi1, theta2, phi2, arcLength, C_GT, C0, j_idx, err_type)

    % -- initialize all gradients --
    grad_k0 = 0;
    grad_tau0 = 0;
    grad_k1 = 0;
    grad_tau1 = 0;    
    grad_theta1 = 0;
    grad_phi1 = 0;
    grad_theta2 = 0;
    grad_phi2 = 0;

    sj = arcLength(j_idx,1);
    
    if strcmp(err_type, "avg_err")
        % -- For the first helix --
        for i = 1:j_idx
            % -- get arc-length --
            s = arcLength(i,1);

            % -- derivative values --
            Df = self_diff_bihelix_func_wrt_curve_geometry_and_frenet_frame(k0, tau0, k1, tau1, s, sj, theta1, phi1, theta2, phi2);

            % -- gradients of k0, tau0, theta1, theta2, and theta3 --
            grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.f1k0 + (C_GT(2,i)-C0(2,i))*Df.f2k0 + (C_GT(3,i)-C0(3,i))*Df.f3k0;
            grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.f1t0 + (C_GT(2,i)-C0(2,i))*Df.f2t0 + (C_GT(3,i)-C0(3,i))*Df.f3t0;
            grad_theta1 = grad_theta1 + (C_GT(1,i)-C0(1,i))*Df.f1theta1 + (C_GT(2,i)-C0(2,i))*Df.f2theta1 + (C_GT(3,i)-C0(3,i))*Df.f3theta1;
            grad_theta2 = grad_theta2 + (C_GT(1,i)-C0(1,i))*Df.f1theta2 + (C_GT(2,i)-C0(2,i))*Df.f2theta2 + (C_GT(3,i)-C0(3,i))*Df.f3theta2;
            grad_phi1 = grad_phi1 + (C_GT(1,i)-C0(1,i))*Df.f1phi1 + (C_GT(2,i)-C0(2,i))*Df.f2phi1 + (C_GT(3,i)-C0(3,i))*Df.f3phi1;
            grad_phi2 = grad_phi2 + (C_GT(1,i)-C0(1,i))*Df.f1phi2 + (C_GT(2,i)-C0(2,i))*Df.f2phi2 + (C_GT(3,i)-C0(3,i))*Df.f3phi2;
        end
        
        % -- For the second helix --
        for i = j_idx+1:size(arcLength, 1)
            % -- get arc-length --
            s = arcLength(i,1);

            % -- derivative values --
            Df = self_diff_bihelix_func_wrt_curve_geometry_and_frenet_frame(k0, tau0, k1, tau1, s, sj, theta1, phi1, theta2, phi2);

            % -- gradients of k0 and tau0 --
            grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.g1k0 + (C_GT(2,i)-C0(2,i))*Df.g2k0 + (C_GT(3,i)-C0(3,i))*Df.g3k0;
            grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.g1t0 + (C_GT(2,i)-C0(2,i))*Df.g2t0 + (C_GT(3,i)-C0(3,i))*Df.g3t0;

            % -- gradients of k1 and tau1 --
            grad_k1 = grad_k1 + (C_GT(1,i)-C0(1,i))*Df.g1k1 + (C_GT(2,i)-C0(2,i))*Df.g2k1 + (C_GT(3,i)-C0(3,i))*Df.g3k1;
            grad_tau1 = grad_tau1 + (C_GT(1,i)-C0(1,i))*Df.g1t1 + (C_GT(2,i)-C0(2,i))*Df.g2t1 + (C_GT(3,i)-C0(3,i))*Df.g3t1;

            % -- gradients of theta, theta2, and theta3 --
            grad_theta1 = grad_theta1 + (C_GT(1,i)-C0(1,i))*Df.g1theta1 + (C_GT(2,i)-C0(2,i))*Df.g2theta1 + (C_GT(3,i)-C0(3,i))*Df.g3theta1;
            grad_theta2 = grad_theta2 + (C_GT(1,i)-C0(1,i))*Df.g1theta2 + (C_GT(2,i)-C0(2,i))*Df.g2theta2 + (C_GT(3,i)-C0(3,i))*Df.g3theta2;
            grad_phi1 = grad_phi1 + (C_GT(1,i)-C0(1,i))*Df.g1phi1 + (C_GT(2,i)-C0(2,i))*Df.g2phi1 + (C_GT(3,i)-C0(3,i))*Df.g3phi1;
            grad_phi2 = grad_phi2 + (C_GT(1,i)-C0(1,i))*Df.g1phi2 + (C_GT(2,i)-C0(2,i))*Df.g2phi2 + (C_GT(3,i)-C0(3,i))*Df.g3phi2;
        end
        
        % -- [AVG ERR] return the average trajectory-curve point-wise error --
        energy_val = 0;
        for i = 1:size(C_GT,2)
            energy_val = energy_val + norm(C_GT(:,i)-C0(:,i));
        end
        energy_val = energy_val / size(C_GT,2);
    else
        % -- For the second helix --
        i = size(arcLength, 1);

        % -- get arc-length --
        s = arcLength(i,1);

        % -- derivative values --
        Df = self_diff_bihelix_func_wrt_curve_geometry_and_frenet_frame(k0, tau0, k1, tau1, s, sj, theta1, phi1, theta2, phi2);

        % -- gradients of k0 and tau0 --
        grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.g1k0 + (C_GT(2,i)-C0(2,i))*Df.g2k0 + (C_GT(3,i)-C0(3,i))*Df.g3k0;
        grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.g1t0 + (C_GT(2,i)-C0(2,i))*Df.g2t0 + (C_GT(3,i)-C0(3,i))*Df.g3t0;

        % -- gradients of k1 and tau1 --
        grad_k1 = grad_k1 + (C_GT(1,i)-C0(1,i))*Df.g1k1 + (C_GT(2,i)-C0(2,i))*Df.g2k1 + (C_GT(3,i)-C0(3,i))*Df.g3k1;
        grad_tau1 = grad_tau1 + (C_GT(1,i)-C0(1,i))*Df.g1t1 + (C_GT(2,i)-C0(2,i))*Df.g2t1 + (C_GT(3,i)-C0(3,i))*Df.g3t1;

        % -- gradients of theta, theta2, and theta3 --
        grad_theta1 = grad_theta1 + (C_GT(1,i)-C0(1,i))*Df.g1theta1 + (C_GT(2,i)-C0(2,i))*Df.g2theta1 + (C_GT(3,i)-C0(3,i))*Df.g3theta1;
        grad_theta2 = grad_theta2 + (C_GT(1,i)-C0(1,i))*Df.g1theta2 + (C_GT(2,i)-C0(2,i))*Df.g2theta2 + (C_GT(3,i)-C0(3,i))*Df.g3theta2;
        grad_phi1 = grad_phi1 + (C_GT(1,i)-C0(1,i))*Df.g1phi1 + (C_GT(2,i)-C0(2,i))*Df.g2phi1 + (C_GT(3,i)-C0(3,i))*Df.g3phi1;
        grad_phi2 = grad_phi2 + (C_GT(1,i)-C0(1,i))*Df.g1phi2 + (C_GT(2,i)-C0(2,i))*Df.g2phi2 + (C_GT(3,i)-C0(3,i))*Df.g3phi2;

        % -- [END-PT ERR] return the end-point error --
        energy_val = norm(C_GT(:,end)-C0(:,end));
    end
end
function [grad_k0, grad_tau0, grad_k1, grad_theta, energy_val] = ...
    self_getGeometryParamsGradientForBihelix_common_tau(k0, tau0, k1, candidate_theta, arcLength, init_N, init_B, tangent0, normal0, binormal0, C_GT, C0, j_idx, err_type)

    grad_k0 = 0;
    grad_tau0 = 0;
    grad_k1 = 0;
    omega0 = sqrt(k0^2 + tau0^2);
    grad_theta = 0;

    sj = arcLength(j_idx,1);
    
    if strcmp(err_type, "avg_err")
        % -- For the first helix --
        for i = 1:j_idx
            % -- get arc-length --
            s = arcLength(i,1);

            % -- derivative values --
            Df = self_diff_bihelix_func_wrt_curve_geometry_common_tau(k0, tau0, k1, s, sj, tangent0, normal0, binormal0);

            % -- gradients of k0 and tau0 --
            grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.f1k0 + (C_GT(2,i)-C0(2,i))*Df.f2k0 + (C_GT(3,i)-C0(3,i))*Df.f3k0;
            grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.f1t0 + (C_GT(2,i)-C0(2,i))*Df.f2t0 + (C_GT(3,i)-C0(3,i))*Df.f3t0;

            % -- gradient of theta --
            coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(i,1)));
            dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
            grad_theta = grad_theta + (C_GT(1,i)-C0(1,i))*dCdtheta(1,1) + (C_GT(2,i)-C0(2,i))*dCdtheta(2,1) + (C_GT(3,i)-C0(3,i))*dCdtheta(3,1);
        end

        % -- For the second helix --
        for i = j_idx+1:size(arcLength, 1)
            % -- get arc-length --
            s = arcLength(i,1);

            % -- derivative values --
            Df = self_diff_bihelix_func_wrt_curve_geometry_common_tau(k0, tau0, k1, s, sj, tangent0, normal0, binormal0);

            % -- gradients of k0 and tau0 --
            grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.g1k0 + (C_GT(2,i)-C0(2,i))*Df.g2k0 + (C_GT(3,i)-C0(3,i))*Df.g3k0;
            grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.g1t0 + (C_GT(2,i)-C0(2,i))*Df.g2t0 + (C_GT(3,i)-C0(3,i))*Df.g3t0;

            % -- gradients of k1 and tau1 --
            grad_k1 = grad_k1 + (C_GT(1,i)-C0(1,i))*Df.g1k1 + (C_GT(2,i)-C0(2,i))*Df.g2k1 + (C_GT(3,i)-C0(3,i))*Df.g3k1;
        end

        % -- [AVG ERR] return the average trajectory-curve point-wise error --
        energy_val = 0;
        for i = 1:size(C_GT,2)
            energy_val = energy_val + sqrt((C_GT(1,i)-C0(1,i))^2 + (C_GT(2,i)-C0(2,i))^2 + (C_GT(3,i)-C0(3,i))^2);
        end
        energy_val = energy_val / size(C_GT,2);
    else
        % -- For the second helix --
        i = size(arcLength, 1);
        
        % -- get arc-length --
        s = arcLength(i,1);

        % -- derivative values --
        Df = self_diff_bihelix_func_wrt_curve_geometry_common_tau(k0, tau0, k1, s, sj, tangent0, normal0, binormal0);

        % -- gradients of k0 and tau0 --
        grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.g1k0 + (C_GT(2,i)-C0(2,i))*Df.g2k0 + (C_GT(3,i)-C0(3,i))*Df.g3k0;
        grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.g1t0 + (C_GT(2,i)-C0(2,i))*Df.g2t0 + (C_GT(3,i)-C0(3,i))*Df.g3t0;

        % -- gradients of k1 and tau1 --
        grad_k1 = grad_k1 + (C_GT(1,i)-C0(1,i))*Df.g1k1 + (C_GT(2,i)-C0(2,i))*Df.g2k1 + (C_GT(3,i)-C0(3,i))*Df.g3k1;

        % -- [END-PT ERR] return the end-point error --
        energy_val = norm(C_GT(:,end)-C0(:,end));
    end
    
end
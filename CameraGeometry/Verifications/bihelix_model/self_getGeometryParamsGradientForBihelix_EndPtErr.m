function [grad_k0, grad_tau0, grad_k1, grad_tau1, grad_theta, energy_val] = ...
    self_getGeometryParamsGradientForBihelix_EndPtErr(k0, tau0, k1, tau1, candidate_theta, arcLength, init_N, init_B, tangent0, normal0, binormal0, C_GT, C0, junction_pt_idx)

    grad_k0 = 0;
    grad_tau0 = 0;
    grad_k1 = 0;
    grad_tau1 = 0;
    
    omega0 = sqrt(k0^2 + tau0^2);
    %omega1 = sqrt(k1^2 + tau1^2);
    
    grad_theta = 0;
    
    sj = arcLength(junction_pt_idx,1);

    for i = 1:junction_pt_idx
        % -- get arc-length --
        s = arcLength(i,1);

        % -- gradient of theta --
        coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(i,1)));
        dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
        grad_theta = grad_theta + (C_GT(1,i)-C0(1,i))*dCdtheta(1,1) + (C_GT(2,i)-C0(2,i))*dCdtheta(2,1) + (C_GT(3,i)-C0(3,i))*dCdtheta(3,1);
    end
    
%     % -- For the first helix, find gradient on the updated angle for normal vector  --
%     i = junction_pt_idx;
%     s = arcLength(i,1);
%     % -- derivative values --
%     Df = self_diff_bihelix_func_wrt_curve_geometry(k0, tau0, k1, tau1, s, sj, tangent0, normal0, binormal0);
%         
%     % -- gradients of k0 and tau0 --
%     grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.f1k0 + (C_GT(2,i)-C0(2,i))*Df.f2k0 + (C_GT(3,i)-C0(3,i))*Df.f3k0;
%     grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.f1t0 + (C_GT(2,i)-C0(2,i))*Df.f2t0 + (C_GT(3,i)-C0(3,i))*Df.f3t0;
    
    
    % -- For the second helix, find the gradients of 2 curvatures, 2 torsions, at the end point --
    %for i = junction_pt_idx+1:size(arcLength, 1)
    i = size(arcLength, 1);
    % -- get arc-length --
    s = arcLength(i,1);

    % -- derivative values --
    Df = self_diff_bihelix_func_wrt_curve_geometry(k0, tau0, k1, tau1, s, sj, tangent0, normal0, binormal0);

    % -- gradients of k0 and tau0 --
    grad_k0 = grad_k0 + (C_GT(1,i)-C0(1,i))*Df.g1k0 + (C_GT(2,i)-C0(2,i))*Df.g2k0 + (C_GT(3,i)-C0(3,i))*Df.g3k0;
    grad_tau0 = grad_tau0 + (C_GT(1,i)-C0(1,i))*Df.g1t0 + (C_GT(2,i)-C0(2,i))*Df.g2t0 + (C_GT(3,i)-C0(3,i))*Df.g3t0;

    % -- gradients of k1 and tau1 --
    grad_k1 = grad_k1 + (C_GT(1,i)-C0(1,i))*Df.g1k1 + (C_GT(2,i)-C0(2,i))*Df.g2k1 + (C_GT(3,i)-C0(3,i))*Df.g3k1;
    grad_tau1 = grad_tau1 + (C_GT(1,i)-C0(1,i))*Df.g1t1 + (C_GT(2,i)-C0(2,i))*Df.g2t1 + (C_GT(3,i)-C0(3,i))*Df.g3t1;
    %end
    
    
    %num_of_pts1 = 1:junction_pt_idx;
    %num_of_pts2 = junction_pt_idx:size(arcLength, 1);
    
    % -- make the average over the gradient --
    %grad_k0 = grad_k0 / size(num_of_pts1, 2);
    %grad_tau0 = grad_tau0 / size(num_of_pts1, 2);
    %grad_k1 = grad_k1 / size(num_of_pts2, 2);
    %grad_tau1 = grad_tau1 / size(num_of_pts2, 2);
    
    % -- [END-PT ERR] return the average trajectory-curve point-wise error --
    energy_val = 0;
    %energy_val = energy_val + norm(C_GT(:,junction_pt_idx)-C0(:,junction_pt_idx));
    energy_val = energy_val + norm(C_GT(:,i)-C0(:,i));
end
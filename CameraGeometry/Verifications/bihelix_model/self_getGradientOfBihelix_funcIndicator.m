function outx = self_getGradientOfBihelix_funcIndicator(x, arcLength, tangent0, normal0, binormal0, C_GT, junction_pt_idx)

    k0 = x(1);
    tau0 = x(2);
    k1 = x(3);
    tau1 = x(4);
    
    [C0, ~] = self_generateBihelixFromModel(k0, tau0, k1, tau1, tangent0, normal0, binormal0, C_GT, junction_pt_idx);

    gradient_curvature1 = 0;
    gradient_torsion1 = 0;
    gradient_curvature2 = 0;
    gradient_torsion2 = 0;
    dCdk0 = zeros(size(arcLength,1), 3);
    dCdtau0 = zeros(size(arcLength,1), 3);
    
    omega0 = sqrt(k0^2 + tau0^2);
    omega1 = sqrt(k1^2 + tau1^2);
    
    gradient_theta = 0;
    
    % -- gradients of k0, tau0, theta --
    for i = 1:junction_pt_idx
        s = arcLength(i,1);

        coeff_T_k0 = (-2*k0*tau0^2/omega0^4)*s + ((2*k0*omega0^2 - 3*k0^3)/omega0^5)*sin(omega0*s) + (k0^3/omega0^4)*cos(omega0*s)*s;
        coeff_N_k0 = ((tau0^2-k0^2) / omega0^4)*(1-cos(omega0*s)) + (k0^2 / omega0^3)*sin(omega0*s)*s;
        coeff_B_k0 = ((tau0^3-k0^2*tau0) / omega0^4)*s - ((tau0*omega0^2-3*tau0*k0^2) / omega0^5)*sin(omega0*s) - (tau0*k0^2 / omega0^4)*cos(omega0*s)*s;

        coeff_T_tau0 = (2*tau0*k0^2 / omega0^4)*s - (3*tau0*k0^2/omega0^5)*sin(omega0*s) + (k0^2*tau0 / omega0^4)*cos(omega0*s)*s;
        coeff_N_tau0 = (-2*tau0*k0 / omega0^4)*(1-cos(omega0*s)) + (k0*tau0 / omega0^3)*sin(omega0*s)*s;
        coeff_B_tau0 = ((k0^3-k0*tau0^2) / omega0^4) - ((k0*omega0^2-3*k0*tau0^2) / omega0^5)*sin(omega0*s) - (k0*tau0^2/omega0^4)*cos(omega0*s)*s;

        dCdk0(i,:) = coeff_T_k0 * tangent0 + coeff_N_k0 * normal0 + coeff_B_k0 * binormal0;
        dCdtau0(i,:) = coeff_T_tau0 * tangent0 + coeff_N_tau0 * normal0 + coeff_B_tau0 * binormal0;

        % -- Frenet frame --
        %coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(i,1)));
        %dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
        %gradient_theta = gradient_theta + (C_GT(1,i)-C0(1,i))*dCdtheta(1,1) + (C_GT(2,i)-C0(2,i))*dCdtheta(2,1) + (C_GT(3,i)-C0(3,i))*dCdtheta(3,1);
    end
    
    % -- propogate the frenet frame and also get the last curve point --
    s1 = arcLength(junction_pt_idx,1);
    % 1) last curve point
    pts_propogation = [(tau0/omega0)^2*s1 + (k0^2/omega0^3)*sin(omega0*s1); ...
                       (k0/(omega0^2))*(1-cos(omega0*s1)); ...
                       (k0*tau0/omega0^2)*s1 - (k0*tau0/omega0^3)*sin(omega0*s1)];
    last_curve_pt = tangent0*pts_propogation(1,1) + normal0*pts_propogation(2,1) + binormal0*pts_propogation(3,1) + C0;
    % 2) propogated Frenet frame
    M = [ (tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*s1), (k0/omega0)*sin(omega0*s1), (k0*tau0/omega0^2)*(1-cos(omega0*s1));
          (-k0/omega0)*sin(omega0*s1), cos(omega0*s1), (tau0/omega0)*sin(omega0*s1);
          (k0*tau0/omega0^2)*(1-cos(omega0*s1)), (-tau0/omega0)*sin(omega0*s1), (k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*s1)
        ];
    propT = M(1,1)*tangent0 + M(1,2)*normal0 + M(1,3)*binormal0;
    propN = M(2,1)*tangent0 + M(2,2)*normal0 + M(2,3)*binormal0;
    propB = M(3,1)*tangent0 + M(3,2)*normal0 + M(3,3)*binormal0;
    
    % -- gradients of k1, tau1 --
    derivative_idx = 1;
    for i = junction_pt_idx+1:size(arcLength, 1)
        s = arcLength(i,1) - arcLength(junction_pt_idx,1);

        coeff_T_k1 = (-2*k1*tau1^2/omega1^4)*s + ((2*k1*omega1^2 - 3*k1^3)/omega1^5)*sin(omega1*s) + (k1^3/omega1^4)*cos(omega1*s)*s;
        coeff_N_k1 = ((tau1^2-k1^2) / omega1^4)*(1-cos(omega1*s)) + (k1^2 / omega1^3)*sin(omega1*s)*s;
        coeff_B_k1 = ((tau1^3-k1^2*tau1) / omega1^4)*s - ((tau1*omega1^2-3*tau1*k1^2) / omega1^5)*sin(omega1*s) - (tau1*k1^2 / omega1^4)*cos(omega1*s)*s;

        coeff_T_tau1 = (2*tau1*k1^2 / omega1^4)*s - (3*tau1*k1^2/omega1^5)*sin(omega1*s) + (k1^2*tau1 / omega1^4)*cos(omega1*s)*s;
        coeff_N_tau1 = (-2*tau1*k1 / omega1^4)*(1-cos(omega1*s)) + (k1*tau1 / omega1^3)*sin(omega1*s)*s;
        coeff_B_tau1 = ((k1^3-k1*tau1^2) / omega1^4) - ((k1*omega1^2-3*k1*tau1^2) / omega1^5)*sin(omega1*s) - (k1*tau1^2/omega1^4)*cos(omega1*s)*s;
        
        dCdk1(derivative_idx,:) = coeff_T_k1 * propT + coeff_N_k1 * propN + coeff_B_k1 * propB;
        dCdtau1(derivative_idx,:) = coeff_T_tau1 * propT + coeff_N_tau1 * propN + coeff_B_tau1 * propB;
        derivative_idx = derivative_idx + 1;
    end

    % -- calculate gradients --
    % -- curve difference is per dimension of curve --
    % 1) first part
    for i = 1:junction_pt_idx
        gradient_curvature1 = gradient_curvature1 + (C_GT(1,i)-C0(1,i))*dCdk0(i,1) + (C_GT(2,i)-C0(2,i))*dCdk0(i,2) + (C_GT(3,i)-C0(3,i))*dCdk0(i,3);
        gradient_torsion1 = gradient_torsion1 + (C_GT(1,i)-C0(1,i))*dCdtau0(i,1) + (C_GT(2,i)-C0(2,i))*dCdtau0(i,2) + (C_GT(3,i)-C0(3,i))*dCdtau0(i,3);
    end
    % 2) second part
    derivative_idx = 1;
    for i = junction_pt_idx+1:size(arcLength, 1)
        gradient_curvature2 = gradient_curvature2 + (C_GT(1,i)-C0(1,i))*dCdk1(derivative_idx,1) + (C_GT(2,i)-C0(2,i))*dCdk1(derivative_idx,2) + (C_GT(3,i)-C0(3,i))*dCdk1(derivative_idx,3);
        gradient_torsion2 = gradient_torsion2 + (C_GT(1,i)-C0(1,i))*dCdtau1(derivative_idx,1) + (C_GT(2,i)-C0(2,i))*dCdtau1(derivative_idx,2) + (C_GT(3,i)-C0(3,i))*dCdtau1(derivative_idx,3);
        derivative_idx = derivative_idx + 1;
    end
    
    num_of_pts1 = 1:junction_pt_idx;
    num_of_pts2 = junction_pt_idx:size(arcLength, 1);
    
    % -- make the average over the gradient --
%     gradient_curvature1 = gradient_curvature1 / size(num_of_pts1, 2);
%     gradient_torsion1 = gradient_torsion1 / size(num_of_pts1, 2);
%     gradient_curvature2 = gradient_curvature2 / size(num_of_pts2, 2);
%     gradient_torsion2 = gradient_torsion2 / size(num_of_pts2, 2);
    
    outx = [gradient_curvature1, gradient_torsion1, gradient_curvature2, gradient_torsion2];
end
function [gradient_curvature, gradient_torsion, gradient_theta, avg_pts_err] = ...
    self_getGeometryParamsGradient(k0, tau0, candidate_theta, arcLength, init_N, init_B, tangent0, normal0, binormal0, C_GT, C0)

    gradient_curvature = 0;
    gradient_torsion = 0;
    dCdk0 = zeros(size(arcLength,1), 3);
    dCdtau0 = zeros(size(arcLength,1), 3);
    
    omega0 = sqrt(k0^2 + tau0^2);
    
    updated_binormal = cross(tangent0, init_N);
    gradient_theta = 0;
    
    for i = 1:size(arcLength,1)
        s = arcLength(i,1);
        
        coeff_T_k0 = (-2*k0*tau0^2/omega0^4)*s + ((2*k0*omega0^2 - 3*k0^3)/omega0^5)*sin(omega0*s) + (k0^3/omega0^4)*cos(omega0*s)*s;
        coeff_N_k0 = ((tau0^2-k0^2) / omega0^4)*(1-cos(omega0*s)) + (k0^2 / omega0^3)*sin(omega0*s)*s;
        coeff_B_k0 = ((tau0^3-k0^2*tau0) / omega0^4)*s - ((tau0*omega0^2-3*tau0*k0^2) / omega0^5)*sin(omega0*s) - (tau0*k0^2 / omega0^4)*cos(omega0*s)*s;

        coeff_T_tau0 = (2*tau0*k0^2 / omega0^4)*s - (3*tau0*k0^2/omega0^5)*sin(omega0*s) + (k0^2*tau0 / omega0^4)*cos(omega0*s)*s;
        coeff_N_tau0 = (-2*tau0*k0 / omega0^4)*(1-cos(omega0*s)) + (k0*tau0 / omega0^3)*sin(omega0*s)*s;
        coeff_B_tau0 = ((k0^3-k0*tau0^2) / omega0^4) - ((k0*omega0^2-3*k0*tau0^2) / omega0^5)*sin(omega0*s) - (k0*tau0^2/omega0^4)*cos(omega0*s)*s;

        dCdk0(i,:) = coeff_T_k0 * tangent0 + coeff_N_k0 * normal0 + coeff_B_k0 * binormal0;
        dCdtau0(i,:) = coeff_T_tau0 * tangent0 + coeff_N_tau0 * normal0 + coeff_B_tau0 * binormal0;
        
%         gradient_curvature = gradient_curvature + (C_GT(1,i)-C0(1,i))*dCdk0(i,1) + (C_GT(2,i)-C0(2,i))*dCdk0(i,2) + (C_GT(3,i)-C0(3,i))*dCdk0(i,3);
%         gradient_torsion = gradient_torsion + (C_GT(1,i)-C0(1,i))*dCdtau0(i,1) + (C_GT(2,i)-C0(2,i))*dCdtau0(i,2) + (C_GT(3,i)-C0(3,i))*dCdtau0(i,3);

        % -- Frenet frame --
        coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(i,1)));
        dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
        gradient_theta = gradient_theta + (C_GT(1,i)-C0(1,i))*dCdtheta(1,1) + (C_GT(2,i)-C0(2,i))*dCdtheta(2,1) + (C_GT(3,i)-C0(3,i))*dCdtheta(3,1);
    end
    
    % -- curve difference is per dimension of curve --
    for i = 1:size(C_GT,2)
        gradient_curvature = gradient_curvature + (C_GT(1,i)-C0(1,i))*dCdk0(i,1) + (C_GT(2,i)-C0(2,i))*dCdk0(i,2) + (C_GT(3,i)-C0(3,i))*dCdk0(i,3);
        gradient_torsion = gradient_torsion + (C_GT(1,i)-C0(1,i))*dCdtau0(i,1) + (C_GT(2,i)-C0(2,i))*dCdtau0(i,2) + (C_GT(3,i)-C0(3,i))*dCdtau0(i,3);
    end
    
    % -- make the average over the gradient --
    gradient_curvature = gradient_curvature / size(arcLength, 1);
    gradient_torsion = gradient_torsion / size(arcLength, 1);
    
    % -- return the average trajectory-curve point-wise error --
    avg_pts_err = 0;
    for i = 1:size(C_GT,2)
        avg_pts_err = avg_pts_err + sqrt((C_GT(1,i)-C0(1,i))^2 + (C_GT(2,i)-C0(2,i))^2 + (C_GT(3,i)-C0(3,i))^2);
    end
    avg_pts_err = avg_pts_err / size(C_GT,2);
end
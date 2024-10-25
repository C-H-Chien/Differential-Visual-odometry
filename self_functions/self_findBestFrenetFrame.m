function [optT, optN, optB] = self_findBestFrenetFrame(init_N, init_T, C_GT, C0, k0, tau0, model_type)
    omega0 = sqrt(k0^2 + tau0^2);
    
    candidate_N = init_N;
    fixed_T = init_T;
    binormal = cross(init_T, init_N);
    init_B = binormal;
    candidate_theta = 0;

    % -- do gradient descent to optimize candidate_theta that candidate_N
    % rotates around fixed_T --
    gd_iter = 500;
    %learning_rate = 0.02;
    learning_rate = 0.5;
    
    % -- define notations for simplicity --
%     Tx = fixed_T(1,1);
%     Ty = fixed_T(2,1);
%     Tz = fixed_T(3,1);
%     skew_sym_T = [0, Tx, Ty; -Tx, 0, Tz; -Ty, -Tz, 0];
    
    collect_avg_errs = zeros(gd_iter, 1);
    
    for i = 1:gd_iter      
        
        % -- generate a curve from model --
        [C_hat, arcLength] = self_generateCurveFromModel(k0, tau0, fixed_T, candidate_N, binormal, C_GT, C0, model_type);
        
        % -- compute the derivative of R*candidate_N with respect to
        % candidate_theta --
%         dRxdtheta = [-sin(candidate_theta)*(Tx^2+Ty^2), cos(candidate_theta)*Tx-sin(candidate_theta)*Ty*Tz, cos(candidate_theta)*Ty+sin(candidate_theta)*Tx*Tz];
%         dRydtheta = [-cos(candidate_theta)*Tx-sin(candidate_theta)*Ty*Tz, -sin(candidate_theta)*(Tx^2+Tz^2), cos(candidate_theta)-sin(candidate_theta)*Tx*Ty];
%         dRzdtheta = [-cos(candidate_theta)*Ty+sin(candidate_theta)*Tx*Tz, -cos(candidate_theta)*Tz-sin(candidate_theta)*Tx*Ty, -sin(candidate_theta)*(Ty^2+Tz^2)];
%         dRNxdtheta = dRxdtheta * candidate_N;
%         dRNydtheta = dRydtheta * candidate_N;
%         dRNzdtheta = dRzdtheta * candidate_N;
        
        gradient_theta = 0;
        
        for j = 1:size(arcLength, 1)
            % -- compute dE/dtheta --
            coefficient = (k0/omega0^2)*(1-cos(omega0*arcLength(j,1)));
            dCdtheta = coefficient*(-sin(candidate_theta)*init_N + cos(candidate_theta)*init_B);
            
            % -- get the gradients!! --
            % -- curve difference is per dimension of curve --
            gradient_theta = gradient_theta + (C_GT(1,j)-C_hat(1,j))*dCdtheta(1,1) + (C_GT(2,j)-C_hat(2,j))*dCdtheta(2,1) + (C_GT(3,j)-C_hat(3,j))*dCdtheta(3,1);
        end
        
        %gradient_theta = gradient_theta / size(arcLength,1);
        candidate_theta = candidate_theta + learning_rate*gradient_theta;
        
        % -- apply Rodrigue's formula to find rotation matrix that rotates candidate_N around fixed_T --
%         R_NaroundT = eye(3) + sin(candidate_theta)*skew_sym_T + (1-cos(candidate_theta))*(skew_sym_T*skew_sym_T);
%         updated_candidate_N = R_NaroundT * init_N;

        % -- with updated candidate_theta, we can update the normal vector
        % --
        updated_candidate_N = cos(candidate_theta)*init_N + sin(candidate_theta)*init_B;
        
        
        validate_dot_product = dot(fixed_T, updated_candidate_N);
        
        % -- prepare for the next GD iteration --
        candidate_N = updated_candidate_N;
        binormal = cross(fixed_T, candidate_N);
        
        % -- calculate the average trajectory-curve point-wise error --
        avg_pts_err = 0;
        for k = 1:size(C_GT,2)
            avg_pts_err = avg_pts_err + sqrt((C_GT(1,k)-C_hat(1,k))^2 + (C_GT(2,k)-C_hat(2,k))^2 + (C_GT(3,k)-C_hat(3,k))^2);
        end
        avg_pts_err = avg_pts_err / size(C_GT,2);
        collect_avg_errs(i,1) = avg_pts_err;
    end

    % -- outputs --
    optT = init_T;
    optN = updated_candidate_N;
    optB = cross(optT, optN);
end

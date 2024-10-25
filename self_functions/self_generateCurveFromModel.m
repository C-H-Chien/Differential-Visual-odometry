function [C, arcLength] = self_generateCurveFromModel(k_0, tau_0, T_0, N_0, B_0, trajectory_pts, curve_start_pt, modelType)
    
    C = zeros(3, size(trajectory_pts,2));
    C_1 = zeros(3, size(trajectory_pts,2));

    % -- compute arc length used by the model curve --
    cumulative_arc_length = zeros(size(trajectory_pts, 2), 1);
    arcLength = zeros(size(trajectory_pts, 2), 1);
    for i = 2:size(trajectory_pts,2)
        pose_dist = sqrt((trajectory_pts(1,i)-trajectory_pts(1,i-1))^2 + (trajectory_pts(2,i)-trajectory_pts(2,i-1))^2 + (trajectory_pts(3,i)-trajectory_pts(3,i-1))^2);
        cumulative_arc_length(i,1) = cumulative_arc_length(i-1,1) + pose_dist;
        arcLength(i,1) = cumulative_arc_length(i,1);
    end

    for i = 1:size(trajectory_pts,2)
        % -- compute arc length used by the model curve --
        ds = arcLength(i,1);

        if strcmp(modelType, 'peicewise_constant_model')
            omega_0 = sqrt(k_0^2 + tau_0^2);

            % -- curve model --
            C_1(:,i) = [(tau_0/omega_0)^2*ds + (k_0^2/omega_0^3)*sin(omega_0*ds); ...
                       (k_0/(omega_0^2))*(1-cos(omega_0*ds)); ...
                       (k_0*tau_0/omega_0^2)*ds - (k_0*tau_0/omega_0^3)*sin(omega_0*ds)];

            % -- aggregate from the start point --
            C(:,i) = curve_start_pt + C_1(1,i)*T_0 + C_1(2,i)*N_0 + C_1(3,i)*B_0;
        end
    end
end
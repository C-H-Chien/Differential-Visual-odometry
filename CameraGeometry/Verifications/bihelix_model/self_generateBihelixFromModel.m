function [C, arcLength] = self_generateBihelixFromModel(k0, tau0, k1, tau1, T_0, N_0, B_0, trajectory_pts, junction_pt_idx)
    
    C = zeros(3, size(trajectory_pts,2));
    %C_1 = zeros(3, size(trajectory_pts,2));
    epsilon = 0.000001;
    omega0 = sqrt(k0^2 + tau0^2);
    omega1 = sqrt(k1^2 + tau1^2);

    % -- compute arc length used by the model curve --
    cumulative_arc_length = zeros(size(trajectory_pts, 2), 1);
    arcLength = zeros(size(trajectory_pts, 2), 1);
    for i = 2:size(trajectory_pts,2)
        pose_dist = sqrt((trajectory_pts(1,i)-trajectory_pts(1,i-1))^2 + (trajectory_pts(2,i)-trajectory_pts(2,i-1))^2 + (trajectory_pts(3,i)-trajectory_pts(3,i-1))^2);
        cumulative_arc_length(i,1) = cumulative_arc_length(i-1,1) + pose_dist;
        arcLength(i,1) = cumulative_arc_length(i,1);
    end

    % -- the start point of the curve, i.e., the first point of trajectory_pts --
    curve_start_pt = trajectory_pts(:,1);
    
    % -- loop over all arclength, plug in the model to create a bihelix curve -- 
    
    % 1) arc length of the junction point --
    js = arcLength(junction_pt_idx, 1);
    
    % 2) curve point at junction
    C_1j = [(tau0/omega0)^2*js + (k0^2/omega0^3)*sin(omega0*js); ...
           (k0/(omega0^2))*(1-cos(omega0*js)); ...
           (k0*tau0/omega0^2)*js - (k0*tau0/omega0^3)*sin(omega0*js)];
    C_j = curve_start_pt + C_1j(1,1)*T_0 + C_1j(2,1)*N_0 + C_1j(3,1)*B_0;
    
    % 3) Frenet frame at junction
    Mj = [ (tau0/omega0)^2 + (k0/omega0)^2*cos(omega0*js), (k0/omega0)*sin(omega0*js), (k0*tau0/omega0^2)*(1-cos(omega0*js));
           (-k0/omega0)*sin(omega0*js), cos(omega0*js), (tau0/omega0)*sin(omega0*js);
           (k0*tau0/omega0^2)*(1-cos(omega0*js)), (-tau0/omega0)*sin(omega0*js), (k0/omega0)^2 + (tau0/omega0)^2*cos(omega0*js)
         ];
    Tj = Mj(1,1)*T_0 + Mj(1,2)*N_0 + Mj(1,3)*B_0;
    Nj = Mj(2,1)*T_0 + Mj(2,2)*N_0 + Mj(2,3)*B_0;
    Bj = Mj(3,1)*T_0 + Mj(3,2)*N_0 + Mj(3,3)*B_0;
    
    for i = 1:size(trajectory_pts,2)
        % -- arc length of each trajectory point --
        ds = arcLength(i,1);

        % -- bihelix model --
        C_1 = [(tau0/omega0)^2*ds + (k0^2/omega0^3)*sin(omega0*ds); ...
               (k0/(omega0^2))*(1-cos(omega0*ds)); ...
               (k0*tau0/omega0^2)*ds - (k0*tau0/omega0^3)*sin(omega0*ds)];
               
        C_2 = [(tau1/omega1)^2*(ds-js) + (k1^2/omega1^3)*sin(omega1*(ds-js)); ...
               (k1/(omega1^2))*(1-cos(omega1*(ds-js))); ...
               (k1*tau1/omega1^2)*(ds-js) - (k1*tau1/omega1^3)*sin(omega1*(ds-js))];
        
        C(:,i) = max(0, (js-ds+epsilon)/abs(js-ds+epsilon))*(curve_start_pt + C_1(1,1)*T_0 + C_1(2,1)*N_0 + C_1(3,1)*B_0) ...
               + max(0, (ds-js)/abs(ds-js))*(C_j + C_2(1,1)*Tj + C_2(2,1)*Nj + C_2(3,1)*Bj);
    end
end
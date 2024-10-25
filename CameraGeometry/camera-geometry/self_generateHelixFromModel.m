function [C, FrenetFrame] = self_generateHelixFromModel(k_0, tau_0, T_0, N_0, B_0, arcLength, curve_start_pt)
    
    omega_0 = sqrt(k_0^2 + tau_0^2);

    C = zeros(3, size(arcLength,2));
    C_1 = zeros(3, size(arcLength,2));
    
    % -- propogate Frenet frame from the first point to the last point --
    ts = arcLength(1,end);
    M = [ (tau_0/omega_0)^2 + (k_0/omega_0)^2*cos(omega_0*ts), (k_0/omega_0)*sin(omega_0*ts), (k_0*tau_0/omega_0^2)*(1-cos(omega_0*ts));
          (-k_0/omega_0)*sin(omega_0*ts), cos(omega_0*ts), (tau_0/omega_0)*sin(omega_0*ts);
          (k_0*tau_0/omega_0^2)*(1-cos(omega_0*ts)), (-tau_0/omega_0)*sin(omega_0*ts), (k_0/omega_0)^2 + (tau_0/omega_0)^2*cos(omega_0*ts)
        ];
    FrenetFrame.T = M(1,1)*T_0 + M(1,2)*N_0 + M(1,3)*B_0;
    FrenetFrame.N = M(2,1)*T_0 + M(2,2)*N_0 + M(2,3)*B_0;
    FrenetFrame.B = M(3,1)*T_0 + M(3,2)*N_0 + M(3,3)*B_0;

    % -- create the circular helix --
    for i = 1:size(arcLength,2)
        % -- arc length of the point i --
        ds = arcLength(1,i);

        % -- curve model --
        C_1(:,i) = [(tau_0/omega_0)^2*ds + (k_0^2/omega_0^3)*sin(omega_0*ds); ...
                   (k_0/(omega_0^2))*(1-cos(omega_0*ds)); ...
                   (k_0*tau_0/omega_0^2)*ds - (k_0*tau_0/omega_0^3)*sin(omega_0*ds)];

        % -- aggregate from the start point --
        C(:,i) = curve_start_pt + C_1(1,i)*T_0 + C_1(2,i)*N_0 + C_1(3,i)*B_0;
    end
end
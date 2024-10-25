function [C, Frenet_frame] = self_create_propogated_Frenet_helix(s, curvature, torsion, init_T, init_N, init_B, init_C)


k_0 = curvature;
tau_0 = torsion;
omega_0 = sqrt(k_0^2 + tau_0^2);

% -- declare propogated T, N, B --
propT = init_T;
propN = init_N;
propB = init_B;

% -- initial cirve point --
C_0 = init_C;

C = zeros(3, size(s,2));
C(:,1) = C_0;
for i = 2:size(s,2)
    % -- 
    T_0 = propT;
    N_0 = propN;
    B_0 = propB;

    % -- arc length difference --
    ds = s(1,i) - s(1,i-1);

    % -- vector C_1 --
    C_1 = [(tau_0/omega_0)^2*ds + (k_0^2/omega_0^3)*sin(omega_0*ds); ...
           (k_0/(omega_0^2))*(1-cos(omega_0*ds)); ...
           (k_0*tau_0/omega_0^2)*ds - (k_0*tau_0/omega_0^3)*sin(omega_0*ds)];

    % -- matrix M --
    M = [ (tau_0/omega_0)^2 + (k_0/omega_0)^2*cos(omega_0*ds), (k_0/omega_0)*sin(omega_0*ds), (k_0*tau_0/omega_0^2)*(1-cos(omega_0*ds));
          (-k_0/omega_0)*sin(omega_0*ds), cos(omega_0*ds), (tau_0/omega_0)*sin(omega_0*ds);
          (k_0*tau_0/omega_0^2)*(1-cos(omega_0*ds)), (-tau_0/omega_0)*sin(omega_0*ds), (k_0/omega_0)^2 + (tau_0/omega_0)^2*cos(omega_0*ds)
        ];

    % -- generated curve points --
    C(:,i) = T_0*C_1(1,1) + N_0*C_1(2,1) + B_0*C_1(3,1) + C_0;
    C_0 = C(:,i);

    % -- propogated Frenet frame --
    propT = M(1,1)*T_0 + M(1,2)*N_0 + M(1,3)*B_0;
    propN = M(2,1)*T_0 + M(2,2)*N_0 + M(2,3)*B_0;
    propB = M(3,1)*T_0 + M(3,2)*N_0 + M(3,3)*B_0;
end

Frenet_frame.tangent = propT;
Frenet_frame.normal = propN;
Frenet_frame.binormal = propB;

end
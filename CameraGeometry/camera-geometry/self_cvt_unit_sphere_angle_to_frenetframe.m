function [T, N, B] = self_cvt_unit_sphere_angle_to_frenetframe(thetaT, phiT, theta_mkortho, phi_mkortho)

    T = [sin(thetaT)*cos(phiT); sin(thetaT)*sin(phiT); cos(thetaT)];
    make_orthogonal_vec = [sin(theta_mkortho)*cos(phi_mkortho); sin(theta_mkortho)*sin(phi_mkortho); cos(theta_mkortho)];
    N = cross(T, make_orthogonal_vec);
    normalized_N = N / norm(N);
    N = normalized_N;
    B = cross(T, N);

end
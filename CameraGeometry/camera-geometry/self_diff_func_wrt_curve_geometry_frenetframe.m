function Df = self_diff_func_wrt_curve_geometry_frenetframe(k, tau, theta1, phi1, theta2, phi2, s)

    c_thetaT = cos(theta1);
    s_thetaT = sin(theta1);
    c_phi_T = cos(phi1);
    s_phi_T = sin(phi1);
    c_theta2 = cos(theta2);
    s_theta2 = sin(theta2);
    c_phi2 = cos(phi2);
    s_phi2 = sin(phi2);
    
    omega = (k^2 + tau^2)^(1/2);

    % -- derivative w.r.t curvature --
    DF1_k = (2*k^2*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2)^2 - (tau*(c_phi_T^2*c_phi2*c_thetaT^2*c_theta2 - c_phi2*c_theta2 + c_phi_T*c_thetaT*s_thetaT*s_theta2 + c_phi_T*c_thetaT^2*c_theta2*s_phi_T*s_phi2)*(tau^2*sin(s*omega) - 2*k^2*sin(s*omega) - s*tau^2*omega + k^2*s*omega + k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - ((c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2) + (k^2*s*sin(s*omega)*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT))/(k^2 + tau^2)^(3/2) - (k*c_phi_T*c_thetaT*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + 2*s*tau^2*omega - k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2);
    DF2_k = ((c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2) - (2*k^2*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2)^2 + (tau*(c_theta2*s_phi2*c_phi_T^2*c_thetaT^2 - c_phi2*c_theta2*s_phi_T*c_phi_T*c_thetaT^2 - s_phi_T*s_theta2*c_thetaT*s_thetaT + c_theta2*s_phi2*s_thetaT^2)*(tau^2*sin(s*omega) - 2*k^2*sin(s*omega) - s*tau^2*omega + k^2*s*omega + k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (k*c_thetaT*s_phi_T*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + 2*s*tau^2*omega - k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (k^2*s*sin(s*omega)*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT))/(k^2 + tau^2)^(3/2);
    DF3_k = (tau*(c_phi_T*c_thetaT*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT) + c_thetaT*s_phi_T*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT))*(tau^2*sin(s*omega) - 2*k^2*sin(s*omega) - s*tau^2*omega + k^2*s*omega + k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (k*s_thetaT*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + 2*s*tau^2*omega - k^2*s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) + (sin(phi1 - phi2)*c_thetaT*c_theta2*(cos(s*omega) - 1))/(k^2 + tau^2) - (2*k^2*sin(phi1 - phi2)*c_thetaT*c_theta2*(cos(s*omega) - 1))/(k^2 + tau^2)^2 - (k^2*s*sin(phi1 - phi2)*sin(s*omega)*c_thetaT*c_theta2)/(k^2 + tau^2)^(3/2);
 
    Df.f1k = DF1_k;
    Df.f2k = DF2_k;
    Df.f3k = DF3_k;
    
    % -- derivative w.r.t. torsion --
    DF1_tau = (2*k*tau*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2)^2 - (k*(c_phi_T^2*c_phi2*c_thetaT^2*c_theta2 - c_phi2*c_theta2 + c_phi_T*c_thetaT*s_thetaT*s_theta2 + c_phi_T*c_thetaT^2*c_theta2*s_phi_T*s_phi2)*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + s*tau^2*omega - k^2*s*omega + s*tau^2*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) + (k*s*tau*sin(s*omega)*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT))/(k^2 + tau^2)^(3/2) + (k^2*tau*c_phi_T*c_thetaT*(2*s*omega - 3*sin(s*omega) + s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2);
    DF2_tau = (k*(c_theta2*s_phi2*c_phi_T^2*c_thetaT^2 - c_phi2*c_theta2*s_phi_T*c_phi_T*c_thetaT^2 - s_phi_T*s_theta2*c_thetaT*s_thetaT + c_theta2*s_phi2*s_thetaT^2)*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + s*tau^2*omega - k^2*s*omega + s*tau^2*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (2*k*tau*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT)*(cos(s*omega) - 1))/(k^2 + tau^2)^2 + (k^2*tau*c_thetaT*s_phi_T*(2*s*omega - 3*sin(s*omega) + s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (k*s*tau*sin(s*omega)*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT))/(k^2 + tau^2)^(3/2);
    DF3_tau = (k*(c_phi_T*c_thetaT*(c_phi_T*c_thetaT*s_theta2 - c_phi2*c_theta2*s_thetaT) + c_thetaT*s_phi_T*(c_thetaT*s_phi_T*s_theta2 - c_theta2*s_phi2*s_thetaT))*(k^2*sin(s*omega) - 2*tau^2*sin(s*omega) + s*tau^2*omega - k^2*s*omega + s*tau^2*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) + (k^2*tau*s_thetaT*(2*s*omega - 3*sin(s*omega) + s*cos(s*omega)*omega))/(k^2 + tau^2)^(5/2) - (2*k*tau*sin(phi1 - phi2)*c_thetaT*c_theta2*(cos(s*omega) - 1))/(k^2 + tau^2)^2 - (k*s*tau*sin(phi1 - phi2)*sin(s*omega)*c_thetaT*c_theta2)/(k^2 + tau^2)^(3/2);

    Df.f1tau = DF1_tau;
    Df.f2tau = DF2_tau;
    Df.f3tau = DF3_tau;
    
    % -- derivative w.r.t. theta1 --
    DF1_theta1 = (k*(s_phi_T*s_thetaT*s_theta2 + c_thetaT*c_theta2*s_phi2)*(cos(s*omega) - 1))/(k^2 + tau^2) - (c_phi_T*s_thetaT*(k^2*sin(s*omega) + s*tau^2*omega))/(k^2 + tau^2)^(3/2) + (k*tau*c_phi_T*(sin(s*omega) - s*omega)*(s_theta2 - 2*c_thetaT^2*s_theta2 + 2*c_phi_T*c_phi2*c_thetaT*c_theta2*s_thetaT + 2*c_thetaT*c_theta2*s_phi_T*s_phi2*s_thetaT))/(k^2 + tau^2)^(3/2);
    DF2_theta1 = (k*tau*(sin(s*omega) - s*omega)*(s_phi_T*s_theta2 - 2*c_thetaT^2*s_phi_T*s_theta2 + 2*c_thetaT*c_theta2*s_phi2*s_thetaT - 2*c_phi_T^2*c_thetaT*c_theta2*s_phi2*s_thetaT + 2*c_phi_T*c_phi2*c_thetaT*c_theta2*s_phi_T*s_thetaT))/(k^2 + tau^2)^(3/2) - (s_phi_T*s_thetaT*(k^2*sin(s*omega) + s*tau^2*omega))/(k^2 + tau^2)^(3/2) - (k*(c_phi_T*s_thetaT*s_theta2 + c_phi2*c_thetaT*c_theta2)*(cos(s*omega) - 1))/(k^2 + tau^2);
    DF3_theta1 = (c_thetaT*(k^2*sin(s*omega) + s*tau^2*omega))/(k^2 + tau^2)^(3/2) - (k*tau*(sin(s*omega) - s*omega)*(2*c_thetaT*s_thetaT*s_theta2 - c_phi_T*c_phi2*c_theta2 - c_theta2*s_phi_T*s_phi2 + 2*c_phi_T*c_phi2*c_thetaT^2*c_theta2 + 2*c_thetaT^2*c_theta2*s_phi_T*s_phi2))/(k^2 + tau^2)^(3/2) - (k*sin(phi1 - phi2)*c_theta2*s_thetaT*(cos(s*omega) - 1))/(k^2 + tau^2);

    Df.f1theta1 = DF1_theta1;
    Df.f2theta1 = DF2_theta1;
    Df.f3theta1 = DF3_theta1;
    
    % -- derivative w.r.t. phi1 --
    DF1_phi1 = (k*tau*c_thetaT*(sin(s*omega) - s*omega)*(s_phi_T*s_thetaT*s_theta2 + c_thetaT*c_theta2*s_phi2 - 2*c_phi_T^2*c_thetaT*c_theta2*s_phi2 + 2*c_phi_T*c_phi2*c_thetaT*c_theta2*s_phi_T))/(k^2 + tau^2)^(3/2) - (k*c_phi_T*c_thetaT*s_theta2*(cos(s*omega) - 1))/(k^2 + tau^2) - (c_thetaT*s_phi_T*(k^2*sin(s*omega) + s*tau^2*omega))/(k^2 + tau^2)^(3/2);
    DF2_phi1 = (c_phi_T*c_thetaT*(k^2*sin(s*omega) + s*tau^2*omega))/(k^2 + tau^2)^(3/2) - (k*tau*c_thetaT*(sin(s*omega) - s*omega)*(c_phi_T*s_thetaT*s_theta2 - c_phi2*c_thetaT*c_theta2 + 2*c_phi_T^2*c_phi2*c_thetaT*c_theta2 + 2*c_phi_T*c_thetaT*c_theta2*s_phi_T*s_phi2))/(k^2 + tau^2)^(3/2) - (k*c_thetaT*s_phi_T*s_theta2*(cos(s*omega) - 1))/(k^2 + tau^2);
    DF3_phi1 = (k*cos(phi1 - phi2)*c_thetaT*c_theta2*(cos(s*omega) - 1))/(k^2 + tau^2) + (k*tau*sin(phi1 - phi2)*sin(2*theta1)*c_theta2*(sin(s*omega) - s*omega))/(2*(k^2 + tau^2)^(3/2));

    Df.f1phi1 = DF1_phi1;
    Df.f2phi1 = DF2_phi1;
    Df.f3phi1 = DF3_phi1;
    
    % -- derivative w.r.t. theta2 --
    DF1_theta2 = (s_thetaT*(c_phi2*s_thetaT*s_theta2 + c_phi_T*c_thetaT*c_theta2) - c_thetaT*s_phi_T*(c_phi_T*c_thetaT*s_phi2*s_theta2 - c_phi2*c_thetaT*s_phi_T*s_theta2))*((k*s*tau)/(k^2 + tau^2) - (k*tau*sin(s*omega))/(k^2 + tau^2)^(3/2)) - (k*(s_phi2*s_thetaT*s_theta2 + c_thetaT*c_theta2*s_phi_T)*(cos(s*omega) - 1))/(k^2 + tau^2);
    DF2_theta2 = (s_thetaT*(s_phi2*s_thetaT*s_theta2 + c_thetaT*c_theta2*s_phi_T) + c_phi_T*c_thetaT*(c_phi_T*c_thetaT*s_phi2*s_theta2 - c_phi2*c_thetaT*s_phi_T*s_theta2))*((k*s*tau)/(k^2 + tau^2) - (k*tau*sin(s*omega))/(k^2 + tau^2)^(3/2)) + (k*(c_phi2*s_thetaT*s_theta2 + c_phi_T*c_thetaT*c_theta2)*(cos(s*omega) - 1))/(k^2 + tau^2);
    DF3_theta2 = (k*tau*(c_phi_T*c_thetaT*(c_phi2*s_thetaT*s_theta2 + c_phi_T*c_thetaT*c_theta2) + c_thetaT*s_phi_T*(s_phi2*s_thetaT*s_theta2 + c_thetaT*c_theta2*s_phi_T))*(sin(s*omega) - s*omega))/(k^2 + tau^2)^(3/2) - (k*sin(phi1 - phi2)*c_thetaT*s_theta2*(cos(s*omega) - 1))/(k^2 + tau^2);
 
    Df.f1theta2 = DF1_theta2;
    Df.f2theta2 = DF2_theta2;
    Df.f3theta2 = DF3_theta2;
    
    % -- derivative w.r.t. phi2 --
    DF1_phi2 = (k*c_phi2*c_theta2*s_thetaT*(cos(s*omega) - 1))/(k^2 + tau^2) - (k*tau*c_theta2*(sin(s*omega) - s*omega)*(- s_phi2*c_phi_T^2*c_thetaT^2 + c_phi2*s_phi_T*c_phi_T*c_thetaT^2 + s_phi2))/(k^2 + tau^2)^(3/2);
    DF2_phi2 = (k*tau*c_theta2*(sin(s*omega) - s*omega)*(c_phi2*c_phi_T^2*c_thetaT^2 + s_phi_T*s_phi2*c_phi_T*c_thetaT^2 + c_phi2*s_thetaT^2))/(k^2 + tau^2)^(3/2) + (k*c_theta2*s_phi2*s_thetaT*(cos(s*omega) - 1))/(k^2 + tau^2);
    DF3_phi2 = -(k*cos(phi1 - phi2)*c_thetaT*c_theta2*(cos(s*omega) - 1))/(k^2 + tau^2) - (k*tau*sin(phi1 - phi2)*sin(2*theta1)*c_theta2*(sin(s*omega) - s*omega))/(2*(k^2 + tau^2)^(3/2));
 
    Df.f1phi2 = DF1_phi2;
    Df.f2phi2 = DF2_phi2;
    Df.f3phi2 = DF3_phi2;

end
function Df = self_diff_func_wrt_curve_geometry(k, tau, s)

    % -- derivative w.r.t curvature --
    Df1_k = -(k*(k^2*sin(s*(k^2 + tau^2)^(1/2)) - 2*tau^2*sin(s*(k^2 + tau^2)^(1/2)) + 2*s*tau^2*(k^2 + tau^2)^(1/2) - k^2*s*cos(s*(k^2 + tau^2)^(1/2))*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(5/2);
    Df2_k = (2*k^2*(cos(s*(k^2 + tau^2)^(1/2)) - 1))/(k^2 + tau^2)^2 - (cos(s*(k^2 + tau^2)^(1/2)) - 1)/(k^2 + tau^2) + (k^2*s*sin(s*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(3/2);
    Df3_k = -(tau*(tau^2*sin(s*(k^2 + tau^2)^(1/2)) - 2*k^2*sin(s*(k^2 + tau^2)^(1/2)) - s*tau^2*(k^2 + tau^2)^(1/2) + k^2*s*(k^2 + tau^2)^(1/2) + k^2*s*cos(s*(k^2 + tau^2)^(1/2))*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(5/2);

    Df.f1k = Df1_k;
    Df.f2k = Df2_k;
    Df.f3k = Df3_k;

    % -- derivative w.r.t. torsion --
    Df1_tau = (k^2*tau*(2*s*(k^2 + tau^2)^(1/2) - 3*sin(s*(k^2 + tau^2)^(1/2)) + s*cos(s*(k^2 + tau^2)^(1/2))*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(5/2);
    Df2_tau = (2*k*tau*(cos(s*(k^2 + tau^2)^(1/2)) - 1))/(k^2 + tau^2)^2 + (k*s*tau*sin(s*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(3/2);
    Df3_tau = -(k*(k^2*sin(s*(k^2 + tau^2)^(1/2)) - 2*tau^2*sin(s*(k^2 + tau^2)^(1/2)) + s*tau^2*(k^2 + tau^2)^(1/2) - k^2*s*(k^2 + tau^2)^(1/2) + s*tau^2*cos(s*(k^2 + tau^2)^(1/2))*(k^2 + tau^2)^(1/2)))/(k^2 + tau^2)^(5/2);

    Df.f1t = Df1_tau;
    Df.f2t = Df2_tau;
    Df.f3t = Df3_tau;

end
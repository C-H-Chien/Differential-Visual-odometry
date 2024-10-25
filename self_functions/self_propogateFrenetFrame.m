function [prop_T, prop_N, prop_B] = self_propogateFrenetFrame(k0, tau0, T_0, N_0, B_0, ds)
    omega0 = sqrt(k0^2 + tau0^2);
    M = [tau0^2/omega0^2 + (k0^2/omega0^2)*cos(omega0*ds), (k0/omega0)*sin(omega0*ds), (k0*tau0/omega0^2)*(1-cos(omega0*ds));
         -(k0/omega0)*sin(omega0*ds), cos(omega0*ds), (tau0/omega0)*sin(omega0*ds);
         (k0*tau0/omega0^2)*(1-cos(omega0*ds)), -(tau0/omega0)*sin(omega0*ds), (k0^2/omega0^2)+(tau0^2/omega0^2)*cos(omega0*ds)];
    
    prop_T = M(1,1)*T_0 + M(1,2)*N_0 + M(1,3)*B_0;
    prop_N = M(2,1)*T_0 + M(2,2)*N_0 + M(2,3)*B_0;
    prop_B = M(3,1)*T_0 + M(3,2)*N_0 + M(3,3)*B_0;
end
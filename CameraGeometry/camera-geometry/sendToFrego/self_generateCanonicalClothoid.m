function [C, FrenetFrame] = self_generateCanonicalClothoid(k_prime, tau_prime, T_0, N_0, B_0, arcLength, curve_start_pt)
    
    % -- initial point and initial Frenet frame --
    FF.T = T_0;
    FF.N = N_0;
    FF.B = B_0;
    C_pt = curve_start_pt;
    C_collections = [];
    
    % -- constants --
    omega = k_prime^2 + tau_prime^2;
    sq_omega = sqrt(omega);
    
    % -- loop over all arc-lengths --
    for i = 1:size(arcLength,2)
        ds = arcLength(i);
        
        % -- calculate Fresnel integrals --
        fresnelSObj = @(t) sin(sq_omega*t.^2/2);
        fresnelCObj = @(t) cos(sq_omega*t.^2/2);
        FS_val = quadgk(fresnelSObj,0,ds,'abstol',1e-15');
        FC_val = quadgk(fresnelCObj,0,ds,'abstol',1e-15');
        
        % -- calculate curve points --
        curve_pt = (k_prime^2 * FC_val + tau_prime^2 * ds) * T_0 ...
                 + (sq_omega * k_prime * FS_val) * N_0 ...
                 + (k_prime * tau_prime * (ds-FC_val)) * B_0;
        curve_pt = curve_pt / omega;
        
        C_collections = [C_collections, curve_pt];
    end
    
    % -- end arc-length --
    end_ds = arcLength(end);
    
    % -- calculate end Frenet frame --
    C = cos(sq_omega*end_ds^2 / 2);
    S = sin(sq_omega*end_ds^2 / 2);
    FF.T = (k_prime^2*C + tau_prime^2)*T_0 + (sq_omega*k_prime*S)*N_0 + (k_prime*tau_prime*(1-C))*B_0;
    FF.N = (-k_prime*S)*T_0 + (sq_omega*C)*N_0 + (tau_prime*S)*B_0;
    FF.B = k_prime*tau_prime*(1-C)*T_0 - (sq_omega*tau_prime*S)*N_0 + (tau_prime^2*C + k_prime^2)*B_0;
    
    FF.T = FF.T / omega;
    FF.N = FF.N / sq_omega;
    FF.B = FF.B / omega;

    % -- outputs --
    C = C_collections;
    FrenetFrame = FF;
end
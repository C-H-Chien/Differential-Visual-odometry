function [comm_G_curve, comm_G_FrenetFrame, FregoPos] = self_generateGeneralCommutatorClothoid(k0, tau0, k_prime, tau_prime, T_0, N_0, B_0, arcLength, curve_start_pt)
    
    % > generate a general 3d clothoid satisfying commutator condition from
    % a canonical 3d clothoid curve

    % -- initial point and initial Frenet frame --
    FF.T = T_0;
    FF.N = N_0;
    FF.B = B_0;
    C_pt = curve_start_pt;
    C_collections = [];
    s_bar = -k0 / k_prime;
    
    s_bar_equivalent = -tau0 / tau_prime;
    if (s_bar ~= s_bar_equivalent)
        fprintf("Commutator Condition Not Satisfied");
        assert(s_bar ~= s_bar_equivalent);
    end
    
    % > intialize format
    p.dk = -k_prime;
    p.dt = -tau_prime;
    
    p.tx0 = FF.T(1,1);
    p.ty0 = FF.T(2,1);
    p.tz0 = FF.T(3,1);
    
    p.nx0 = FF.N(1,1);
    p.ny0 = FF.N(2,1);
    p.nz0 = FF.N(3,1);
    
    p.bx0 = FF.B(1,1);
    p.by0 = FF.B(2,1);
    p.bz0 = FF.B(3,1);
    
    p.x0 = 0;
    p.y0 = 0;
    p.z0 = 0;
    [pos,t,n,b] = clothoid3d(-s_bar, p);
    
    % > generate the frenet frame and the offset point 
    [pt_C, FrenetFrame] = self_generateCanonicalClothoid(-k_prime, -tau_prime, FF.T, FF.N, FF.B, -s_bar, [0;0;0]);
    [offset_pt, ~]   = self_generateCanonicalClothoid(k_prime, tau_prime, FrenetFrame.T, FrenetFrame.N, FrenetFrame.B, -s_bar, [0;0;0]);
    
    %> check
    p.dk = k_prime;
    p.dt = tau_prime;
    
    p.tx0 = t(1,1);
    p.ty0 = t(2,1);
    p.tz0 = t(3,1);
    
    p.nx0 = n(1,1);
    p.ny0 = n(2,1);
    p.nz0 = n(3,1);
    
    p.bx0 = b(1,1);
    p.by0 = b(2,1);
    p.bz0 = b(3,1);
    
    p.x0 = 0;
    p.y0 = 0;
    p.z0 = 0;
    [pos,~,~,~] = clothoid3d(-s_bar, p);
    
    % =======
    
    
    offset_s = arcLength - s_bar;
    offset_init_pt = curve_start_pt - offset_pt;
    [comm_G_curve, comm_G_FrenetFrame] = self_generateCanonicalClothoid(k_prime, tau_prime, FrenetFrame.T, FrenetFrame.N, FrenetFrame.B, offset_s, offset_init_pt);

    %> check
    p.dk = k_prime;
    p.dt = tau_prime;
    %p.dk = k0;
    %p.dt = tau0;
    
    offset_s = arcLength - s_bar;
    offset_init_pt = curve_start_pt - pos;
    
    p.x0 = offset_init_pt(1,1);
    p.y0 = offset_init_pt(2,1);
    p.z0 = offset_init_pt(3,1);
    [FregoPos,t,n,b] = clothoid3d(offset_s, p);

end
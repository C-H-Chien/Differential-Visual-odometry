function [C, FrenetFrame] = self_generateClothoidFromModel(k0, tau0, k_prime, tau_prime, T_0, N_0, B_0, arcLength, curve_start_pt)
    
    % -- constants --
    alpha = (3-2*sqrt(3))/12;
    beta = (3+2*sqrt(3))/12;
    h = arcLength(1,2) - arcLength(1,1);
    
    FF.T = T_0;
    FF.N = N_0;
    FF.B = B_0;
    C_pt = curve_start_pt;
    C_collections = [];
    
    % -- loop over all arc-lengths --
    for i = 0:size(arcLength,2)-1
        s = i*h;
        v1 = s + (0.5 - sqrt(3)/6)*h;
        v2 = s + (0.5 + sqrt(3)/6)*h;
        c1 = alpha*(h*k_prime*v1 + k0) + beta*(h*k_prime*v2 + k0);
        c2 = alpha*(h*tau_prime*v1 + tau0) + beta*(h*tau_prime*v2 + tau0);
        c = sqrt(c1^2 + c2^2);
        d1 = beta*(h*k_prime*v1 + k0) + alpha*(h*k_prime*v2 + k0);
        d2 = beta*(h*tau_prime*v1 + tau0) + alpha*(h*tau_prime*v2 + tau0);
        d = sqrt(d1^2 + d2^2);
        
        exp_B1 = [c*(c1^2 * cos(h*c)+c2^2),              c^2 * c1*sin(h*c),              c*c1*c2*(1-cos(h*c)),              0; ...
                  -c^2 * c1*sin(h*c),                    c^3 * cos(h*c),                 c^2 * c2*sin(h*c),                 0; ...
                  c*c1*c2*(1-cos(h*c)),                 -c^2 * c2*sin(h*c),              c*(c2^2 * cos(h*c) + c1^2),        0; ...
                  (alpha+beta)*(h*c*c2^2+c1^2*sin(h*c)), (alpha+beta)*c*c1*(1-cos(h*c)), (alpha+beta)*c1*c2*(h*c-sin(h*c)), c^3    ];
        
        exp_B2 = [d*(d1^2 * cos(h*d)+d2^2),              d^2 * d1*sin(h*d),              d*d1*d2*(1-cos(h*d)),              0; ...
                  -d^2 * d1*sin(h*d),                    d^3 * cos(h*d),                 d^2 * d2*sin(h*d),                 0; ...
                  d*d1*d2*(1-cos(h*d)),                 -d^2 * d2*sin(h*d),              d*(d2^2 * cos(h*d) + d1^2),        0; ...
                  (alpha+beta)*(h*d*d2^2+d1^2*sin(h*d)), (alpha+beta)*d*d1*(1-cos(h*d)), (alpha+beta)*d1*d2*(h*d-sin(h*d)), d^3    ];
        
        exp_B1 = exp_B1.*(1/c^3);
        exp_B2 = exp_B2.*(1/d^3);
                
        transition_matrix = exp_B1*exp_B2;
        
        tangent = transition_matrix(1,1).*FF.T + transition_matrix(1,2).*FF.N + transition_matrix(1,3).*FF.B;
        normal = transition_matrix(2,1).*FF.T + transition_matrix(2,2).*FF.N + transition_matrix(2,3).*FF.B;
        binormal = transition_matrix(3,1).*FF.T + transition_matrix(3,2).*FF.N + transition_matrix(3,3).*FF.B;
        curve_pt = transition_matrix(4,1).*FF.T + transition_matrix(4,2).*FF.N + transition_matrix(4,3).*FF.B + transition_matrix(4,4).*C_pt;
        
        FF.T = tangent;
        FF.N = normal;
        FF.B = binormal;
        C_pt = curve_pt;
        
        C_collections = [C_collections, curve_pt];
    end

    % -- outputs --
    C = C_collections;
    FrenetFrame = FF;
    
end

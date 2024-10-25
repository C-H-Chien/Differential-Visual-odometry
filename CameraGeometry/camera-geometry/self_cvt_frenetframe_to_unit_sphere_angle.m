function [thetaT, phiT, theta_mkortho, phi_mkortho] = self_cvt_frenetframe_to_unit_sphere_angle(T, N)

    thetaT = acos(T(3,1));
    phiT = asin(T(2,1)/sin(thetaT));
    
    syms theta2 phi2;
%     eq1 = cos(thetaT)*sin(phiT)*sin(theta2) - cos(theta2)*sin(phi2)*sin(thetaT) == N(1,1);
%     eq2 = cos(phi2)*cos(theta2)*sin(thetaT) - cos(phiT)*cos(thetaT)*sin(theta2) == N(2,1);
%     E = [eq1, eq2];
%     S = solve(E, theta2, phi2);
    
%     TcrossM = [cos(theta2)*sin(phiT)*sin(thetaT) - cos(thetaT)*sin(phi2)*sin(theta2); ...
%                cos(phi2)*cos(thetaT)*sin(theta2) - cos(phiT)*cos(theta2)*sin(thetaT); ...
%                cos(phiT)*sin(phi2)*sin(thetaT)*sin(theta2) - cos(phi2)*sin(phiT)*sin(thetaT)*sin(theta2)
%               ];
    
    % -- function 1 --
    TcrossM = [cos(theta2)*sin(phiT)*sin(thetaT) - cos(thetaT)*sin(phi2)*sin(theta2); ...
               cos(phi2)*cos(thetaT)*sin(theta2) - cos(phiT)*cos(theta2)*sin(thetaT)
              ];
    
    TcM_fcn = matlabFunction(TcrossM, 'Vars',{[theta2, phi2]});
    [th_vct, ~, ~, out_struct] = fsolve(@(in1) TcM_fcn(in1)-N(1:2,1), 2*pi*rand(1,2));
    
    out_msg = string(out_struct.message);
    if ~contains(out_msg, "Equation solved")
        % -- alternative function, funtion 2 --
        TcrossM = [cos(theta2)*sin(phiT)*sin(thetaT) - cos(thetaT)*sin(phi2)*sin(theta2); ...
                   cos(phiT)*sin(phi2)*sin(thetaT)*sin(theta2) - cos(phi2)*sin(phiT)*sin(thetaT)*sin(theta2)
                  ];

        TcM_fcn = matlabFunction(TcrossM, 'Vars',{[theta2, phi2]});
        [th_vct, ~, ~, out_struct] = fsolve(@(in1) TcM_fcn(in1)-N(2:3,1), 2*pi*rand(1,2));
        out_msg = string(out_struct.message);
    end
    
    if ~contains(out_msg, "Equation solved")
        % -- alternative function, funtion 3 --
        TcrossM = [cos(theta2)*sin(phiT)*sin(thetaT) - cos(thetaT)*sin(phi2)*sin(theta2); ...
                   cos(phiT)*sin(phi2)*sin(thetaT)*sin(theta2) - cos(phi2)*sin(phiT)*sin(thetaT)*sin(theta2)
                  ];

        TcM_fcn = matlabFunction(TcrossM, 'Vars',{[theta2, phi2]});
        rhs_N = [N(1,1); N(3,1)];
        th_vct = fsolve(@(in1) TcM_fcn(in1)-rhs_N, 2*pi*rand(1,2));
    end
    
    % -- the solutions --
    theta_mkortho = double(th_vct(1));
    phi_mkortho = double(th_vct(2));
end
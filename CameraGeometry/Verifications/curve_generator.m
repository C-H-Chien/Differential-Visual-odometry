
function [C, speed, a0, curvature, torsion, jerk, alpha_t] = curve_generator(start_fr, end_fr, all_T, time_per_fr, peicewise_constant_model, curve_type, dynamic_type)

    start_time = time_per_fr(1,start_fr);
    r = all_T(:,start_fr:end_fr);
    t = time_per_fr(1,start_fr:end_fr);
    t = t - start_time;

    % -- find the best curve that fits the discrete data --
    % -- C(s) only depends on T, independent of R --
    Cx_t = [r(1,:); t];
    Cy_t = [r(2,:); t];
    Cz_t = [r(3,:); t];

    degree = 3;
    poly_x = polyfit(Cx_t(2,:), Cx_t(1,:), degree);
    poly_y = polyfit(Cy_t(2,:), Cy_t(1,:), degree);
    poly_z = polyfit(Cz_t(2,:), Cz_t(1,:), degree);

    % -- polynomial evaluation --
    %eval_poly_x = poly_x(1)*Cx_t(2,:).^3 + poly_x(2)*Cx_t(2,:).^2 + poly_x(3)*Cx_t(2,:) + poly_x(4);

    % -- curve parametrized by t --
    alpha_t = zeros(size(r,2), 3);
    alpha_t(:,1) = poly_x(1)*Cx_t(2,:).^3 + poly_x(2)*Cx_t(2,:).^2 + poly_x(3)*Cx_t(2,:) + poly_x(4);
    alpha_t(:,2) = poly_y(1)*Cy_t(2,:).^3 + poly_y(2)*Cy_t(2,:).^2 + poly_y(3)*Cy_t(2,:) + poly_y(4);
    alpha_t(:,3) = poly_z(1)*Cz_t(2,:).^3 + poly_z(2)*Cz_t(2,:).^2 + poly_z(3)*Cz_t(2,:) + poly_z(4);

    % -- compute velocity --

    vel = zeros(size(r,2), 3);
    speed = zeros(size(r,2), 1);
    vel(:,1) = 3*poly_x(1)*Cx_t(2,:).^2 + 2*poly_x(2)*Cx_t(2,:) + poly_x(3);
    vel(:,2) = 3*poly_y(1)*Cy_t(2,:).^2 + 2*poly_y(2)*Cy_t(2,:) + poly_y(3);
    vel(:,3) = 3*poly_z(1)*Cz_t(2,:).^2 + 2*poly_z(2)*Cz_t(2,:) + poly_z(3);
%     vel(:,1) = 2*poly_x(1)*Cx_t(2,:) + poly_x(2);
%     vel(:,2) = 2*poly_y(1)*Cy_t(2,:) + poly_y(2);
%     vel(:,3) = 2*poly_z(1)*Cz_t(2,:) + poly_z(2);
    for i = 1:size(vel, 1)
        speed(i,1) = norm(vel(i,:));
    end

    % -- compute acceleration --
    % -- acceleration is the change of speed over time --
    a0 = zeros(size(vel,1), 1);
    for i = 1:size(vel, 1)
        if i == 1
            a0(i,1) = 0;
        else
            a0(i,1) = (speed(i,1) - speed(i-1,1))/(t(1,i)-t(1,i-1));
        end
    end

    % -- acceleration in each direction --
    acc = zeros(size(vel,1), 3);
    mag_acc = zeros(size(vel,1), 1);
    acc(:,1) = 6*poly_x(1)*Cx_t(2,:) + 2*poly_x(2);
    acc(:,2) = 6*poly_y(1)*Cy_t(2,:) + 2*poly_y(2);
    acc(:,3) = 6*poly_z(1)*Cz_t(2,:) + 2*poly_z(2);
%     acc(:,1) = 2*poly_x(1);
%     acc(:,2) = 2*poly_y(1);
%     acc(:,3) = 2*poly_z(1);
    for i = 1:size(acc, 1)
        mag_acc(i, 1) = norm(acc(i,:));
    end

    acc_prime = zeros(size(vel,1), 3);
%     acc_prime(:,1) = 6*poly_x(1);
%     acc_prime(:,2) = 6*poly_y(1);
%     acc_prime(:,3) = 6*poly_z(1);
    for i = 1:size(acc, 1)
        mag_acc_prime(i, 1) = norm(acc_prime(i,:));
    end

    % -- compute jerk which is the change of acceleration over time --
    jerk = zeros(size(vel,1), 1);
    for i = 1:size(vel, 1)
        if i == 1
            jerk(i,1) = 0;
        else
            jerk(i,1) = (a0(i,1) - a0(i-1,1))/(t(1,i)-t(1,i-1));
        end
    end

    % -- curvature and torsion --
    curvature = zeros(size(vel,1), 1);
    torsion = zeros(size(vel,1), 1);

    for i = 1:size(vel,1)
        numer = norm(cross(vel(i,:), acc(i,:)));
        denom = norm(vel(i,:))^3;
        curvature(i,1) = numer / denom;

        numer = dot(cross(vel(i,:), acc(i,:)), acc_prime(i,:));
        denom = (norm(cross(vel(i,:), acc(i,:))))^2;
        torsion(i,1) = -numer / denom;
    end

    % -- reconstruct curve C(s) --
    % -- redefining notations --
    if ~peicewise_constant_model
        k_0 = curvature;
        tau_0 = torsion;

        omega_0 = zeros(size(curvature));
        for i = 1:size(curvature, 1)
            omega_0(i,1) = sqrt(k_0(i,1)^2 + tau_0(i,1)^2);
        end
    end

    s_second_order = zeros(size(r,2), 1);
    s_third_order = s_second_order;
    for i = 1:size(r,2)
%         s_second_order(i,1) = speed(i,1)*t(1,i) + 0.5*a0(i,1)*t(1,i)^2;
%         s_third_order(i,1) = speed(i,1)*t(1,i) + 0.5*a0(i,1)*t(1,i)^2 + (jerk(i,1)*t(1,i)^3)/6;
        
        s_second_order(i,1) = speed(3,1)*t(1,i) + 0.5*a0(3,1)*t(1,i)^2;
        s_third_order(i,1) = speed(3,1)*t(1,i) + 0.5*a0(3,1)*t(1,i)^2 + (jerk(3,1)*t(1,i)^3)/6;
    end

    % -- find T, N, B --
    dsdt = zeros(size(r,2), 1);
    for i = 1:size(r,2)
        dsdt(i,1) = norm(vel(i,:));
    end

    init_i = 1;
    T_0 = vel(init_i,:)./dsdt(init_i,1);
    N_0 = acc(init_i,:)./norm(acc(init_i,:));
    B_0 = cross(T_0, N_0);

    % ================================================================================
    C = zeros(3, size(r,2));
    C_1 = zeros(3, size(r,2));
    for i = init_i:size(r,2)
        if strcmp(dynamic_type, 'second_order')
            ds = s_second_order(i,1);
        elseif strcmp(dynamic_type, 'third_order')
            ds = s_third_order(i,1);
        end

        if peicewise_constant_model
            % -- make curvature and torsion constant across N frames --
            %if mod(i, N+1) == 0 || i == init_i
            if i == init_i
                k_0 = curvature(i+1,1);
                tau_0 = torsion(i+1,1);
                omega_0 = sqrt(k_0^2 + tau_0^2);
            end
            
            C_1(:,i) = [(tau_0/omega_0)^2*ds + (k_0^2/omega_0^3)*sin(omega_0*ds); ...
                       (k_0/(omega_0^2))*(1-cos(omega_0*ds)); ...
                       (k_0*tau_0/omega_0^2)*ds - (k_0*tau_0/omega_0^3)*sin(omega_0*ds)];
            if strcmp(curve_type, 'exact')
                C(:,i) = r(:,1) + C_1(1,i)*T_0' + C_1(2,i)*N_0' + C_1(3,i)*B_0';
            elseif strcmp(curve_type, '1st-order-approx')
                C(:,i) = r(:,1) + ds*T_0' + (0.5*k_0*ds^2)*N_0';
            elseif strcmp(curve_type, '2nd-order-approx')
                C(:,i) = r(:,1) + (ds-(k_0^2*ds^3)/6)*T_0' + (0.5*k_0*ds^2)*N_0' + (k_0*tau_0/6)*ds^3*B_0';
            end
        else
            C_1(:,i) = [(tau_0(i,1)/omega_0(i,1))^2*ds + (k_0(i,1)^2/omega_0(i,1)^3)*sin(omega_0(i,1)*ds); ...
               (k_0(i,1)/(omega_0(i,1)^2))*(1-cos(omega_0(i,1)*ds)); ...
               (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^2)*ds - (k_0(i,1)*tau_0(i,1)/omega_0(i,1)^3)*sin(omega_0(i,1)*ds)];
           
            if strcmp(curve_type, 'exact')
                C(:,i) = r(:,1) + C_1(1,i)*T_0' + C_1(2,i)*N_0' + C_1(3,i)*B_0';
            elseif strcmp(curve_type, '1st-order-approx')
                C(:,i) = r(:,1) + ds*T_0' + (0.5*k_0(i,1)*ds^2)*N_0';
            elseif strcmp(curve_type, '2nd-order-approx')
                C(:,i) = r(:,1) + (ds-(k_0(i,1)^2*ds^3)/6)*T_0' + (0.5*k_0(i,1)*ds^2)*N_0' + (k_0(i,1)*tau_0(i,1)/6)*ds^3*B_0';
            end
        end
    end
end

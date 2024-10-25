function [dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(start_fr, end_fr, all_T, time_per_fr, dim, dynamicParams, geoemtryParams)

    start_time = time_per_fr(1,start_fr+2);
    r = all_T(:,start_fr:end_fr+2);
    t = time_per_fr(1,start_fr:end_fr+2);
    t = t - start_time;
    
    if strcmp(dim, '3d')
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
        acc_prime(:,1) = 6*poly_x(1);
        acc_prime(:,2) = 6*poly_y(1);
        acc_prime(:,3) = 6*poly_z(1);
        for i = 1:size(acc, 1)
            mag_acc_prime(i, 1) = norm(acc_prime(i,:));
        end

        % -- compute jolt which is the change of acceleration over time --
        jolt = zeros(size(vel,1), 1);
        for i = 1:size(vel, 1)
            if i == 1
                jolt(i,1) = 0;
            else
                jolt(i,1) = (a0(i,1) - a0(i-1,1))/(t(1,i)-t(1,i-1));
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
        
        dynamicParams.speed = speed;
        dynamicParams.a0 = a0;
        dynamicParams.jolt = jolt;
        dynamicParams.vel_vec = vel;
        dynamicParams.acc_vec = acc;
        
        geometryParams.curvature = curvature;
        geometryParams.torsion = torsion;
        
        % -- find T, N, B --
        dsdt = zeros(size(r,2), 1);
        for i = 1:size(r,2)
            dsdt(i,1) = norm(vel(i,:));
        end

        init_i = 3;
        % -- tangent --
        T_0 = vel(init_i,:)./dsdt(init_i,1);
        % -- normal --
        %derivative_of_T = acc(init_i,:)./dsdt(init_i,1);
        
        %N_0 = derivative_of_T./norm(derivative_of_T);
        %B_0 = cross(T_0, N_0);
        
        geometryParams.tangent = T_0;
        %geometryParams.normal = N_0;
        %geometryParams.binormal = B_0;
        
        curvature_prime = zeros(size(vel,1), 1);
        torsion_prime = zeros(size(vel,1), 1);
        for i = 1:size(vel, 1)
            if i == 1
                curvature_prime(i,1) = 0;
                torsion_prime(i,1) = 0;
            else
                curvature_prime(i,1) = (curvature(i,1) - curvature(i-1,1))/(t(1,i)-t(1,i-1));
                torsion_prime(i,1) = (torsion(i,1) - torsion(i-1,1))/(t(1,i)-t(1,i-1));
            end
        end
        
        geometryParams.curvature_derivative = curvature_prime;
        geometryParams.torsion_derivative = torsion_prime;
        
    elseif strcmp(dim, '2d')
        % -- find the best curve that fits the discrete data --
        % -- C(s) only depends on T, independent of R --
        Cx_t = [r(1,:); t];
        Cy_t = [r(2,:); t];

        degree = 3;
        poly_x = polyfit(Cx_t(2,:), Cx_t(1,:), degree);
        poly_y = polyfit(Cy_t(2,:), Cy_t(1,:), degree);

        % -- polynomial evaluation --
        %eval_poly_x = poly_x(1)*Cx_t(2,:).^3 + poly_x(2)*Cx_t(2,:).^2 + poly_x(3)*Cx_t(2,:) + poly_x(4);

        % -- curve parametrized by t --
        alpha_t = zeros(size(r,2), 3);
        alpha_t(:,1) = poly_x(1)*Cx_t(2,:).^3 + poly_x(2)*Cx_t(2,:).^2 + poly_x(3)*Cx_t(2,:) + poly_x(4);
        alpha_t(:,2) = poly_y(1)*Cy_t(2,:).^3 + poly_y(2)*Cy_t(2,:).^2 + poly_y(3)*Cy_t(2,:) + poly_y(4);

        % -- compute velocity --
        vel = zeros(size(r,2), 3);
        speed = zeros(size(r,2), 1);
        vel(:,1) = 3*poly_x(1)*Cx_t(2,:).^2 + 2*poly_x(2)*Cx_t(2,:) + poly_x(3);
        vel(:,2) = 3*poly_y(1)*Cy_t(2,:).^2 + 2*poly_y(2)*Cy_t(2,:) + poly_y(3);
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
        for i = 1:size(acc, 1)
            mag_acc(i, 1) = norm(acc(i,:));
        end

        % -- compute jolt which is the change of acceleration over time --
        jolt = zeros(size(vel,1), 1);
        for i = 1:size(vel, 1)
            if i == 1
                jolt(i,1) = 0;
            else
                jolt(i,1) = (a0(i,1) - a0(i-1,1))/(t(1,i)-t(1,i-1));
            end
        end

        % -- curvature and torsion --
        curvature = zeros(size(vel,1), 1);
        torsion = zeros(size(vel,1), 1);

        for i = 1:size(vel,1)
            numer = norm(cross(vel(i,:), acc(i,:)));
            denom = norm(vel(i,:))^3;
            curvature(i,1) = numer / denom;

            torsion(i,1) = 0;
        end
        
        dynamicParams.speed = speed;
        dynamicParams.a0 = a0;
        dynamicParams.jolt = jolt;
        dynamicParams.vel_vec = vel;
        dynamicParams.acc_vec = acc;
        
        geometryParams.curvature = curvature;
        geometryParams.torsion = torsion;
    end
    
end

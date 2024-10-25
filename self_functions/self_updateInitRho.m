function [updated_covis_tracks, err_per_point, neg_indices] = self_updateInitRho(covis_tracks, R_cc, T_cc, K)
    
    % -- compute the optimal rho_0 analytically --
    covis_frame_num = size(covis_tracks, 3);
    updated_covis_tracks = covis_tracks(:,1:3,:);
    e3 = [0,0,1];
    
    point_poly_roots = zeros(size(covis_tracks, 1), 4);

    neg_indices = [];
    val_covis_tracks = covis_tracks;

    % -- loop over all points --
    for i = 1:size(val_covis_tracks, 1)

        % -- compute the derivative of the objective function --
        pix_0 = ones(3,1);
        pix_0(1, 1) = val_covis_tracks(i,1,1);
        pix_0(2, 1) = val_covis_tracks(i,2,1);
        gamma_0 = K \ pix_0;

        % -- use the last two frames of pairs of points to compute --
        fr_t1 = covis_frame_num;
        fr_t2 = covis_frame_num-1;

        % -- compute observations \hat{\gamma(t1)} and \hat{\gamma(t2)} --
        pix_t1 = [val_covis_tracks(i,1,fr_t1); val_covis_tracks(i,2,fr_t1); 1];
        pix_t2 = [val_covis_tracks(i,1,fr_t2); val_covis_tracks(i,2,fr_t2); 1];
        hat_gamma_t1 = K \ pix_t1;
        hat_gamma_t2 = K \ pix_t2;
        
        % -- compute coefficients of the biquadratic function --
        c1c2_coterm = (-e3*T_cc(:,:,fr_t1))*R_cc(:,:,fr_t1)'*gamma_0 + T_cc(:,:,fr_t1)*(e3*(R_cc(:,:,fr_t1)'*gamma_0));
        c1 = ((e3*T_cc(:,:,fr_t1))*hat_gamma_t1 - T_cc(:,:,fr_t1))'*c1c2_coterm;
        c2 = (-e3*(R_cc(:,:,fr_t1)'*gamma_0)*hat_gamma_t1 + R_cc(:,:,fr_t1)'*gamma_0)'*c1c2_coterm;
        
        c3c4_coterm = (-e3*T_cc(:,:,fr_t2))*R_cc(:,:,fr_t2)'*gamma_0 + T_cc(:,:,fr_t2)*(e3*(R_cc(:,:,fr_t2)'*gamma_0));
        c3 = ((e3*T_cc(:,:,fr_t2))*hat_gamma_t2 - T_cc(:,:,fr_t2))'*c3c4_coterm;
        c4 = (-e3*(R_cc(:,:,fr_t2)'*gamma_0)*hat_gamma_t2 + R_cc(:,:,fr_t2)'*gamma_0)'*c3c4_coterm;
        
        p1 = e3*(R_cc(:,:,fr_t1)'*gamma_0);
        p2 = e3*(R_cc(:,:,fr_t2)'*gamma_0);
        p3 = -e3*T_cc(:,:,fr_t1);
        p4 = -e3*T_cc(:,:,fr_t2);
        
        a = c2*(p2^3) + c4*(p1^3);
        b = c1*(p2^3) + 3*c2*(p2^2)*p4 + c3*(p1^3) + 3*c4*(p1^2)*p3;
        c = 3*c1*(p2^2)*p4 + 3*c2*p2*(p4^2) + 3*c3*(p1^2)*p3 + 3*c4*p1*(p3^2);
        d = 3*c1*p2*(p4^2) + c2*(p4^3) + 3*c3*p1*(p3^2) + c4*(p3^3);
        e = c1*p4^3 + c3*p3^3;
      
        % -- solve the biquadratic equation --
        eq_coeff = [a b c d e];
        r = roots(eq_coeff);
        point_poly_roots(i, :) = r';
        
        % -- extract the valid rho_0 (positive and real) among 4 solutions --
        for n = 1:size(r, 2)
            if imag(r(1,n)) == 0 && real(r(1,n)) > 0
                % -- store the analytically solved rho_0 --
                updated_covis_tracks(i,3, covis_frame_num) = r(1,n);
                break;
            end
        end

        rho_0 = updated_covis_tracks(i,3, covis_frame_num);
        % -- use the analytically solved rho_0 to compute new gamma(t) --
        if rho_0 > 0
            for fr = 2:covis_frame_num
                gamma_numer = rho_0.*(R_cc(:,:,fr)'*gamma_0) - T_cc(:,:,fr);
                gamma_denom = rho_0.*(e3*(R_cc(:,:,fr)'*gamma_0)) - e3*T_cc(:,:,fr);

                gamma_t = gamma_numer ./ gamma_denom;
                gamma_t_pix = (K * gamma_t);
                updated_covis_tracks(i,1:2,fr) = gamma_t_pix(1:2,1);
            end
        else
            neg_indices = [neg_indices; i];
        end
    end

    % -- remove points with invalid depths --
    updated_covis_tracks(neg_indices,:,:) = [];
    val_covis_tracks(neg_indices,:,:) = [];

    % -- compute the error between new gamma(t) and feature track --
    err_over_fr = 0;
    err_per_point = zeros(size(val_covis_tracks, 1), 1);
    for i = 1:size(val_covis_tracks, 1)
        euclidean_dist = 0;
        for fr = 1:covis_frame_num
            x_err = (val_covis_tracks(i,1,fr) - updated_covis_tracks(i,1,fr))^2;
            y_err = (val_covis_tracks(i,2,fr) - updated_covis_tracks(i,2,fr))^2;
            euclidean_dist = euclidean_dist + sqrt(x_err + y_err);
        end
        err_per_point(i,1) = euclidean_dist / covis_frame_num;
        err_over_fr = err_over_fr + euclidean_dist / covis_frame_num;
    end

    
end
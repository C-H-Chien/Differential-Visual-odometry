function [covis_tracks, f2, vpts2, rho_0, indexPairs] = ...
    self_get_initial_rho(keyPointsObj_prev, f2, f_prev, vpts2, vpts_prev, Rprime, Tprime, time_difference, covis_tracks, K)
    
    % -- initialization --
    e_1 = [1,0,0];
    e_2 = [0,1,0];
    e_3 = [0,0,1];
    f_size = size(keyPointsObj_prev.Location,1);
    prev_covis_tracks = zeros(f_size, 2);

    % -- store the locations of features on the start_fr previous img --
    prev_covis_tracks(:,1) = keyPointsObj_prev.Location(:,1);
    prev_covis_tracks(:,2) = keyPointsObj_prev.Location(:,2);

    % -- match the previsou (start_fr-1) frame with the second frame --
    indexPairs = matchFeatures(f2,f_prev);
    matchedPoints_f2 = vpts2(indexPairs(:,1));
    matchedPoints_prev = vpts_prev(indexPairs(:,2));
    
    % -- delete features in covis_tracks that are not matched --
    covis_tracks = covis_tracks(indexPairs(:,1),:,:);
    f2 = f2(indexPairs(:,1),:);
    vpts2 = vpts2(indexPairs(:,1),:);
    prev_covis_tracks = prev_covis_tracks(indexPairs(:,2),:);
    
    % -- positive x and y counter --
    pos_rho_0_x = 0;
    pos_rho_0_y = 0;
    
%     figure;
%     plot(prev_covis_tracks(:,2), prev_covis_tracks(:,1), 'r+');
%     hold on;
%     plot(covis_tracks(:,2,2), covis_tracks(:,1,2), 'b+');

    rho_0_from_x = zeros(size(covis_tracks, 1), 1);
    rho_0_from_y = zeros(size(covis_tracks, 1), 1);
    for i = 1:size(covis_tracks, 1)
        % -- compute gamma from feature point pixel coordinates on the start frame (start_fr) --
        pix_f1 = ones(3,1);
        pix_f1(1, 1) = covis_tracks(i,1,1);
        pix_f1(2, 1) = covis_tracks(i,2,1);
        gamma_start_fr = K \ pix_f1;
        
        % -- compute gamma from feature point pixel coordinates on the second frame (f2) --
        pix_f2 = ones(3,1);
        pix_f2(1, 1) = covis_tracks(i,1,2);
        pix_f2(2, 1) = covis_tracks(i,2,2);
        gamma_1 = K \ pix_f2;
        
        % -- compute gamma from feature point pixel coordinates on the previous frame --
        pix_prev = ones(3,1);
        pix_prev(1,1) = prev_covis_tracks(i,1);
        pix_prev(2,1) = prev_covis_tracks(i,2);
        gamma_prev = K \ pix_prev;

        % -- compute the derivative of gamma_0 using central difference (NOT USED) --
        gamma_x_prime = (gamma_1(1,1) - gamma_prev(1,1)) / time_difference;
        gamma_y_prime = (gamma_1(2,1) - gamma_prev(2,1)) / time_difference;
        gamma_prime_0 = zeros(3,1);
        
        % -- rho_0 can be computed from either u or v coordinate of
        % features, so I compute both and select rho_0 which gives maximal
        % number of valid depths: --        
        % -- u-coordinate: compute rho_0 for the feature point on the first frame --
        rho_numer = -e_1*Tprime + (e_3*Tprime)*gamma_start_fr(1,1);
        rho_denom = (e_3*Rprime'*gamma_start_fr)*gamma_start_fr(1,1) - e_1*Rprime'*gamma_start_fr + gamma_x_prime;
        rho_0_from_x(i,1) = rho_numer./rho_denom;
        if rho_0_from_x(i,1) > 0
            pos_rho_0_x = pos_rho_0_x + 1;
        end

        % -- v-coordinate: compute rho_0 for the feature point on the first frame using --
        rho_numer = -e_2*Tprime + (e_3*Tprime)*gamma_start_fr(2,1);
        rho_denom = (e_3*Rprime'*gamma_start_fr)*gamma_start_fr(2,1) - e_2*Rprime'*gamma_start_fr + gamma_y_prime;
        rho_0_from_y(i,1) = rho_numer./rho_denom;
        if rho_0_from_y(i,1) > 0
            pos_rho_0_y = pos_rho_0_y + 1;
        end
    end
    
    % -- pick the one with maximal set of feature points --
    if pos_rho_0_x >= pos_rho_0_y
        rho_0 = rho_0_from_x;
    else
        rho_0 = rho_0_from_y;
    end

    %covis_tracks(:,3, 1) = rho_0;
end
function self_drawColoredFeatureTracks_5_tracks_with_img(projectedImgNames, covisible_frames, start_fr, covisible_feature_tracks, reprojected_feature_tracks, plot_track_num)
    

    % Determine the number of tracks to plot
    num_tracks = size(covisible_feature_tracks, 1);
    
    % Compute distances between the first and last feature points of each track
    distances = zeros(num_tracks, 1);
    for i = 1:num_tracks
        start_point = covisible_feature_tracks(i, :, 1);
        end_point = covisible_feature_tracks(i, :, covisible_frames);
        distances(i) = sqrt((end_point(1) - start_point(1))^2 + (end_point(2) - start_point(2))^2);
    end
    
    % Sort tracks by distance and select the top 10
    [~, sorted_indices] = sort(distances, 'descend');
    if num_tracks > 10
        selected_indices = sorted_indices(1:10);
        covisible_feature_tracks = covisible_feature_tracks(selected_indices, :, :);
        reprojected_feature_tracks = reprojected_feature_tracks(selected_indices, :, :);
        num_tracks = 10;
    else
        covisible_feature_tracks = covisible_feature_tracks(sorted_indices, :, :);
        reprojected_feature_tracks = reprojected_feature_tracks(sorted_indices, :, :);
    end
    
    % Create a colormap - MATLAB automatically handles up to about 60 distinct colors well
    colors = lines(num_tracks);  % 'lines' is just one option, 'jet', 'hsv', 'cool', etc., can also be used

    
    % -- compute u_vec and v_vec of feature tracks --
    u_vec_feature = zeros(size(covisible_feature_tracks, 1), covisible_frames);
    v_vec_feature = zeros(size(covisible_feature_tracks, 1), covisible_frames);
    for j = 2:covisible_frames
        u_vec_feature(:, j) = covisible_feature_tracks(:, 1, j) - covisible_feature_tracks(:, 1, j-1);
        v_vec_feature(:, j) = covisible_feature_tracks(:, 2, j) - covisible_feature_tracks(:, 2, j-1);
    end

     % -- compute u_vec and v_vec of reprojected feature tracks --
    u_vec_feature_reprojected = zeros(size(reprojected_feature_tracks, 1), covisible_frames);
    v_vec_feature_reprojected = zeros(size(reprojected_feature_tracks, 1), covisible_frames);
    for j = 2:covisible_frames
        u_vec_feature_reprojected(:, j) = reprojected_feature_tracks(:, 1, j) - reprojected_feature_tracks(:, 1, j-1);
        v_vec_feature_reprojected(:, j) = reprojected_feature_tracks(:, 2, j) - reprojected_feature_tracks(:, 2, j-1);
    end


    % Plotting feature tracks with different colors
    for i = 1:covisible_frames
        figure;
        imshow(projectedImgNames(start_fr+i));
        hold on;  % Move hold on outside loop for better performance
        
        for k = 1:num_tracks

            % Check if track should be plotted (not all zeros)
            if any(covisible_feature_tracks(k, :, i) ~= 0)
                plot(covisible_feature_tracks(k, 1, i), covisible_feature_tracks(k, 2, i), 'g+', 'Color', 'yellow');
                if i >= 2
                    quiver(covisible_feature_tracks(k, 1, i-1), covisible_feature_tracks(k, 2, i-1), u_vec_feature(k, i), v_vec_feature(k, i), 0, 'Color', colors(k, :));
                end
            end
            % Check if track should be plotted (not all zeros)
            if any(reprojected_feature_tracks(k, :, i) ~= 0)
                plot(reprojected_feature_tracks(k, 1, i), reprojected_feature_tracks(k, 2, i), 'o-', 'Color', 'yellow');
                if i >= 2
                    quiver(reprojected_feature_tracks(k, 1, i-1), reprojected_feature_tracks(k, 2, i-1), u_vec_feature_reprojected(k, i), v_vec_feature_reprojected(k, i), 0, 'Color', colors(k, :));
                end
            end

            pause(0.1);

        end
    end
    
   
  

end

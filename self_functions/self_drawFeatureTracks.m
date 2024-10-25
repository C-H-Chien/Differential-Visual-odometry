% -- draw and connect features that projected to the first frame to form feature tracks --
function self_drawFeatureTracks(projectedImgName, covisible_frames, covisible_feature_tracks, covisible_gamma_tracks, plot_track_num)
    figure;
    imshow(projectedImgName);
    
    % -- compute u_vec and v_vec of feature tracks --
    u_vec_feature = zeros(size(covisible_feature_tracks,1), covisible_frames);
    v_vec_feature = zeros(size(covisible_feature_tracks,1), covisible_frames);
    for j = 2:covisible_frames
        u_vec_feature(:,j) = covisible_feature_tracks(:,1,j) - covisible_feature_tracks(:,1,j-1);
        v_vec_feature(:,j) = covisible_feature_tracks(:,2,j) - covisible_feature_tracks(:,2,j-1);
    end
    
%     % -- compute u_vec and v_vec of gamma tracks --
%     u_vec_gamma = zeros(size(covisible_feature_tracks,1), covisible_frames);
%     v_vec_gamma = zeros(size(covisible_feature_tracks,1), covisible_frames);
%     for j = 2:covisible_frames
%         u_vec_gamma(:,j) = covisible_gamma_tracks(:,1,j) - covisible_gamma_tracks(:,1,j-1);
%         v_vec_gamma(:,j) = covisible_gamma_tracks(:,2,j) - covisible_gamma_tracks(:,2,j-1);
%     end
    
    for i = 1:covisible_frames
        if plot_track_num == -1

            % -- plote feature tracks --
            hold on;
            plot(covisible_feature_tracks(:,1,i), covisible_feature_tracks(:,2,i), 'g+');
            if i >= 2
                hold on;
                quiver(covisible_feature_tracks(:,1,i-1), covisible_feature_tracks(:,2,i-1), u_vec_feature(:,i), v_vec_feature(:,i), 0, 'g');
            end

            % -- plot gamma(t) --
%             hold on;
%             plot(covisible_gamma_tracks(:,1,i), covisible_gamma_tracks(:,2,i), 'ro');
%             if i >= 2
%                 hold on;
%                 quiver(covisible_gamma_tracks(:,1,i-1), covisible_gamma_tracks(:,2,i-1), u_vec_gamma(:,i), v_vec_gamma(:,i), 0, 'r');
%             end
        else        
           
            % -- plot feature tracks --
            hold on;
            plot(covisible_feature_tracks(1:plot_track_num,1,i), covisible_feature_tracks(1:plot_track_num,2,i), 'g+');
            if i >= 2
                hold on;
                quiver(covisible_feature_tracks(1:plot_track_num,1,i-1), covisible_feature_tracks(1:plot_track_num,2,i-1), u_vec_feature(1:plot_track_num,i), v_vec_feature(1:plot_track_num,i), 0, 'g');
            end

            % -- plot gamma(t) --
            hold on;
            plot(covisible_gamma_tracks(1:plot_track_num,1,i), covisible_gamma_tracks(1:plot_track_num,2,i), 'ro');
            if i >= 2
                hold on;
                quiver(covisible_gamma_tracks(1:plot_track_num,1,i-1), covisible_gamma_tracks(1:plot_track_num,2,i-1), u_vec_gamma(1:plot_track_num,i), v_vec_gamma(1:plot_track_num,i), 0, 'r');
            end

        end
    end
end
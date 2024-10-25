function plot_geometry_dynamics_params(speed, acceleration, curvature, torsion, jolt, start_fr, end_fr)

    % -- derivative of curvature --

    % ----------- PLOTTINGS -------------
    % -- plot Gaussian convolved speed --
    figure;
    subplot(1,5,1);
    plot(start_fr:end_fr, speed(3:end, 1), 'b-');
    xlim([start_fr, end_fr]);
    title({'speed'});
    xlabel("points");
    ylabel("speed (m/s)");

    % -- plot Gaussian convolved acceleration --
    subplot(1,5,2);
    plot(start_fr:end_fr, acceleration(3:end, 1), 'b-');
    xlim([start_fr, end_fr]);
    %ylim([0, 3.5]);
    title('acceleration');
    xlabel("points");
    ylabel("acceleration (m/s^2)");
    
    % -- plot thrid derivative of curve --
    subplot(1,5,3);
    plot(start_fr:end_fr, jolt(3:end, 1), 'b-');
%     upper_bound = round(max(jolt(3:end, 1)), 1);
%     lower_bound = round(min(jolt(3:end, 1)), 1);
%     if upper_bound == lower_bound
%         
%     end
    ylim([-1, 2]);
    xlim([start_fr, end_fr]);
    title('jolt');
    xlabel("points");
    ylabel("jolt (m/s^3)");

    % -- plot Gaussian convolved curvature --
    subplot(1,5,4);
    plot(start_fr:end_fr, curvature(3:end, 1), 'b-');
    xlim([start_fr, end_fr]);
    %ylim([-0.15, 0.15]);
    title('curvature');
    xlabel("points");
    ylabel("curvature (1/m)");

    % -- plot Gaussian convolved torsion --
    subplot(1,5,5);
    plot(start_fr:end_fr, torsion(3:end, 1), 'b-');
    xlim([start_fr, end_fr]);
    ylim([-1, 1]);
    title('torsion');
    xlabel("points");
    ylabel("torsion");

    
    axis equal;
    set(gcf,'color','w');
end


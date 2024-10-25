function self_drawGeometryModelFit(t, GT_curvature, GT_torsion, curvature_prime, torsion_prime, opt_curvature, opt_torsion)

    % ----------- PLOTTINGS -------------
    % -- curvature --
    opt_curvature_copies = opt_curvature*ones(size(t(3:end,1), 1), 1);
    subplot(1,4,1);
    plot(t(3:end,1), GT_curvature(:,1), 'b-');
    hold on;
    plot(t(3:end,1), opt_curvature_copies, 'r-');
    xlim([t(3,1), t(end,1)]);
    title({'curvature'});
    xlabel("time");
    ylabel("curvature");
    
    % -- torsion --
    opt_torsion_copies = opt_torsion*ones(size(t(3:end,1), 1), 1);
    subplot(1,4,2);
    plot(t(3:end,1), GT_torsion(:,1), 'b-');
    hold on;
    plot(t(3:end,1), opt_torsion_copies, 'r-'); 
    xlim([t(3,1), t(end,1)]);
    title('torsion');
    xlabel("time");
    ylabel("torsion");

    % -- curvature prime --
    subplot(1,4,3);
    plot(t(3:end,1), curvature_prime(:, 1), 'b-');
    hold on;
    plot(t(3:end,1), zeros(size(t(3:end,1), 1), 1), 'r-');
    xlim([t(3,1), t(end,1)]);
    title({'curvature','derivative', '        '});
    xlabel("time");
    ylabel("curvature derivative");
    
    % -- torsion prime --
    subplot(1,4,4);
    plot(t(3:end,1), torsion_prime(:, 1), 'b-', 'DisplayName', 'ground-truth');
    hold on;
    plot(t(3:end,1), zeros(size(t(3:end,1), 1), 1), 'r-', 'DisplayName', 'best-fit');
    xlim([t(3,1), t(end,1)]);
    title({'torsion','derivative', '     '});
    xlabel("time");
    ylabel("torsion derivative");
    
    %legend;
end
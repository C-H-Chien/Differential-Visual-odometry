function drawDynamicModelFit(t, arcLength, speed, acceleration, jolt, colorLabel)

    % ----------- PLOTTINGS -------------
    % -- s(t) --
    subplot(1,4,1);
    plot(t(3:end,1), arcLength(:, 1), 'bo-');
    hold on;
    plot(t(3:end,1), arcLength(:, 2), 'ro-');
    hold on;
    plot(t(3:end,1), arcLength(:, 3), 'go-');
    xlim([t(3,1), t(end,1)]);
    title({'s(t)'});
    xlabel("time");
    ylabel("arc-length (m)");
    
    % -- speed --
    subplot(1,4,2);
    plot(t(3:end,1), speed(:, 1), colorLabel(1));
    hold on;
    plot(t(3:end,1), speed(:, 2), colorLabel(2));
    hold on;
    plot(t(3:end,1), speed(:, 3), colorLabel(3));
    xlim([t(3,1), t(end,1)]);
    title('v(t)');
    xlabel("time");
    ylabel("speed (m/s)");
    axis equal;

    % -- acceleration --
    subplot(1,4,3);
    plot(t(3:end,1), acceleration(:, 1), colorLabel(1));
    hold on;
    plot(t(3:end,1), acceleration(:, 2), colorLabel(2));
    hold on;
    plot(t(3:end,1), acceleration(:, 3), colorLabel(3));
    xlim([t(3,1), t(end,1)]);
    %ylim([0, 3.5]);
    title('a(t)');
    xlabel("time");
    ylabel("acceleration (m/s^2)");
    %axis equal;
    
    % -- jolt --
    subplot(1,4,4);
    plot(t(3:end,1), jolt(:, 1), colorLabel(1), 'DisplayName', 'ground truth');
    hold on;
    plot(t(3:end,1), jolt(:, 2), colorLabel(2), 'DisplayName', '2nd-order best fit');
    hold on;
    plot(t(3:end,1), jolt(:, 3), colorLabel(3), 'DisplayName', '3rd-order best fit');
    %ylim([-1, 2]);
    xlim([t(3,1), t(end,1)]);
    title('j(t)');
    xlabel("time");
    ylabel("jolt (m/s^3)");
    %axis equal;
    
    %legend;
    set(gcf,'color','w');
end
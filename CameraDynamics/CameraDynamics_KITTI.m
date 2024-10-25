

clc; clear all; close all;

%> define parameters
datasetName = 'KITTI';
sequenceName = '10';
category_left = '/image_0/';

%> read R, T, and time stamp of each frame
[all_T, all_R, time_per_fr, ~] = readDataset(datasetName, sequenceName, category_left);

start_fr = 2;
end_fr = 12;
seq_end_fr = size(all_T,2)-50;
window_size = 5;

run_segment_wise = 1;
debugging = 1;

%> dynamic model params and geometry model params
dynamicParams = [];
geometryParams = [];

%> drawings and printings
print_segment_wise = run_segment_wise;
drawErrChangeAlongSequence = ~run_segment_wise;

enable_drawDynamicModelFit = run_segment_wise;
drawTrajectory = run_segment_wise;
drawAvgErrHistogram = 1;

drawHistogram = 0;  %> temporarily not used
drawArcLengthCompare = 0;
piecewise_constant = 1;
slide_length = 0;

if drawTrajectory
    %> plot trajectory
    figure;
    p_ground_truth = plot3(all_T(1,start_fr:end_fr), all_T(2,start_fr:end_fr), all_T(3,start_fr:end_fr), 'bo--');
    hold on;
    plot3(all_T(1,start_fr), all_T(2,start_fr), all_T(3,start_fr), 'ko');
    text(all_T(1,start_fr), all_T(2,start_fr), all_T(3,start_fr), 'start', 'FontSize', 10, 'Color', 'k');
    hold on;
    plot3(all_T(1,end_fr), all_T(2,end_fr), all_T(3,end_fr), 'ko');
    text(all_T(1,end_fr), all_T(2,end_fr), all_T(3,end_fr), 'end', 'FontSize', 10, 'Color', 'k');
    xlabel("x (m)");
    ylabel("y (m)");
    zlabel("z (m)");
    axis equal;
    set(gcf,'color','w');
end

if run_segment_wise
    end_fr = start_fr;
else
    end_fr = seq_end_fr;
end

%> record all end-point errors
collection_of_end_point_errs = zeros(size(start_fr:seq_end_fr,2), 2);
collection_of_arcLength_diff = zeros(size(start_fr:seq_end_fr,2), 1);
collection_of_err_percentage = zeros(size(start_fr:seq_end_fr,2), 2);

%> record all average errors
collection_of_avg_errs = zeros(size(start_fr:seq_end_fr,2), 2);
record_err_idx = 1;
for fr_idx = start_fr:end_fr

    [dynamicParams, ~, r, t] = self_dynamics_geometry_params_extractor(fr_idx, fr_idx+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);
    speed = dynamicParams.speed;
    a0 = dynamicParams.a0;
    jolt = dynamicParams.jolt;

    %> compute the arc length of the curve using second or third order dynamics
    s_GT = zeros(size(speed,1), 1);
    s_second_order = zeros(size(speed,1)-2, 1);
    s_third_order = s_second_order;

    %> start the arc-length from 0 for the ground truths
    cumulative_arc_length = zeros(size(s_GT, 1), 1);
    for i = 4:size(speed,1)
        pose_dist = sqrt((r(1,i)-r(1,i-1))^2 + (r(2,i)-r(2,i-1))^2 + (r(3,i)-r(3,i-1))^2);
        cumulative_arc_length(i,1) = cumulative_arc_length(i-1,1) + pose_dist;
        s_GT(i,1) = cumulative_arc_length(i,1);
    end
    arcLength_total_diff = s_GT(end,1);

    %> not used
    for i = 4:size(speed,1)
        %> dynamic models
        s_second_order(i,1) = speed(3,1)*t(1,i) + 0.5*a0(3,1)*t(1,i)^2;
        s_third_order(i,1) = speed(3,1)*t(1,i) + 0.5*a0(3,1)*t(1,i)^2 + (jolt(3,1)*t(1,i)^3)/6;
    end

    %> fit the s_GT(t) using polyfit 
    %> s(t) = v_0*t + 0.5*a_0*t^2
    second_order_coeffs = polyfit(t(1,3:end), s_GT(3:end,1), 2);
    third_order_coeffs = polyfit(t(1,3:end), s_GT(3:end,1), 3);

    %> print out the best fit velocity, acceleration, and jolt
    if print_segment_wise
        fprintf('second-order best fit:\nv0=%f,\ta0=%f\n', second_order_coeffs(2), second_order_coeffs(1)*2);
        fprintf('third-order best fit:\nv0=%f,\ta0=%f,\tj0=%f\n', third_order_coeffs(3), third_order_coeffs(2)*2, 6*third_order_coeffs(1));
    end

    s_second_order_fit = zeros(size(s_GT,1)-2, 1);
    s_third_order_fit = zeros(size(s_GT,1)-2, 1);
    predicted_second_order_vel = zeros(size(s_GT,1)-2, 1);
    predicted_second_order_acc = zeros(size(s_GT,1)-2, 1);
    predicted_third_order_vel = zeros(size(s_GT,1)-2, 1);
    predicted_third_order_acc = zeros(size(s_GT,1)-2, 1);
    predicted_second_order_jolt = zeros(size(s_GT,1)-2, 1);
    predicted_third_order_jolt = zeros(size(s_GT,1)-2, 1);
    fit_idx = 1;
    for p = 3:size(s_GT, 1)
        %> second order
        s_second_order_fit(fit_idx,1) = second_order_coeffs(1)*t(1,p)^2 + second_order_coeffs(2)*t(1,p) + second_order_coeffs(3);
        s_second_order_fit(fit_idx,1) = s_second_order_fit(fit_idx,1) - second_order_coeffs(3);
        predicted_second_order_vel(fit_idx,1) = second_order_coeffs(2) + 2*second_order_coeffs(1)*t(1,p);
        predicted_second_order_acc(fit_idx,1) = 2*second_order_coeffs(1);

        %> third order
        s_third_order_fit(fit_idx,1) = third_order_coeffs(1)*t(1,p)^3 + third_order_coeffs(2)*t(1,p)^2 + third_order_coeffs(3)*t(1,p) + third_order_coeffs(4);
        s_third_order_fit(fit_idx,1) = s_third_order_fit(fit_idx,1) - third_order_coeffs(4);
        predicted_third_order_vel(fit_idx,1) = 3*third_order_coeffs(1)*t(1,p)^2 + 2*third_order_coeffs(2)*t(1,p) + third_order_coeffs(3);
        predicted_third_order_acc(fit_idx,1) = 6*third_order_coeffs(1)*t(1,p) + 2*third_order_coeffs(2);
        predicted_third_order_jolt(fit_idx,1) = 6*third_order_coeffs(1);

        fit_idx = fit_idx + 1;
    end

    %> compute end-point error
    endpoint_err_second_order = abs(s_GT(end,1) - s_second_order_fit(end,1));
    endpoint_err_third_order = abs(s_GT(end,1) - s_third_order_fit(end,1));

    %> compute average point-pairs error
    sum_over_pt_pair_err_second_order = 0;
    sum_over_pt_pair_err_third_order = 0;
    for p_idx = 1:size(s_second_order_fit, 1)
        sum_over_pt_pair_err_second_order = sum_over_pt_pair_err_second_order + (s_GT(p_idx+2,1) - s_second_order_fit(p_idx,1))^2;
        sum_over_pt_pair_err_third_order = sum_over_pt_pair_err_third_order + (s_GT(p_idx+2,1) - s_third_order_fit(p_idx,1))^2;
    end
    collection_of_avg_errs(record_err_idx, 1) = sqrt(sum_over_pt_pair_err_second_order) / size(s_second_order_fit, 1);
    collection_of_avg_errs(record_err_idx, 2) = sqrt(sum_over_pt_pair_err_third_order) / size(s_second_order_fit, 1);

    collection_of_end_point_errs(record_err_idx, 1) = endpoint_err_second_order;
    collection_of_end_point_errs(record_err_idx, 2) = endpoint_err_third_order;
    collection_of_arcLength_diff(record_err_idx, 1) = arcLength_total_diff;
    collection_of_err_percentage(record_err_idx, 1) = (endpoint_err_second_order / arcLength_total_diff)*100;
    collection_of_err_percentage(record_err_idx, 2) = (endpoint_err_third_order / arcLength_total_diff)*100;
    record_err_idx = record_err_idx + 1;

    if mod(record_err_idx, 20)==0
        fprintf('. ');
    elseif mod(record_err_idx, 100) == 0
        fprintf('\n');
    end
end

if print_segment_wise
    fprintf('2nd-order end point error: %f\n', endpoint_err_second_order);
    fprintf('3rd-order end point error: %f\n', endpoint_err_third_order);
end

%> drawings
if enable_drawDynamicModelFit    
    %> stacking
    arcLength_stack = [s_GT(3:end,1), s_second_order_fit, s_third_order_fit];
    speed_stack = [speed(3:end,1), predicted_second_order_vel, predicted_third_order_vel];
    acceleration_stack = [a0(3:end,1), predicted_second_order_acc, predicted_third_order_acc];
    jolt_stack = [jolt(3:end,1), predicted_second_order_jolt, predicted_third_order_jolt];
    colors_stack = ["b-", "r-", "g-"];

    %> draw out
    figure;
    drawDynamicModelFit(t', arcLength_stack, speed_stack, acceleration_stack, jolt_stack, colors_stack);    
    set(gcf,'color','w');
end

if drawHistogram
    %[binN,~] = histcounts(avg_err_second_order);
    h1 = histogram(avg_err_second_order, 'DisplayName', 'second-order');
    h1.BinWidth = 0.02;
    hold on;
    h2 = histogram(avg_err_third_order, 'DisplayName', 'third-order');
    h2.BinWidth = 0.02;
    xlabel('error (m)');
    ylabel('number of cases');
    xlim([-0.05, 0.8]);
    str_title = strcat(datasetName, {' '}, 'sequence', {' '}, sequenceName);
    title(str_title);
    legend;
    set(gcf,'color','w');
end

if drawArcLengthCompare
    %fprintf('2nd-order dynamic avg error: %f\n', avg_err_second_order);
    %fprintf('3rd-order dynamic avg error: %f\n', avg_err_third_order);

    %> plot the results
    figure;
    frame_index = start_fr:start_fr+N;
    plot(t(1,:), s_GT(1:end,1), 'bo-', 'DisplayName', 'ground truth');
    hold on;
    plot(t(1,:), s_second_order(1:end,1), 'ro-', 'DisplayName', 'second order dynamic');
    hold on;
    plot(t(1,:), s_third_order(1:end,1), 'go-', 'DisplayName', 'third order dynamic');
    legend;
    ylabel('arc-length');
    xlabel('frame index');
    set(gcf,'color','w');
end

if drawErrChangeAlongSequence
    figure;
    subplot(3,1,1);
    plot(start_fr:seq_end_fr, collection_of_arcLength_diff(:,1), 'b-');
    xlabel('frame index');
    ylabel({'end-point ground','truth arcLength','(m)'});
    xlim([start_fr, seq_end_fr]);

    subplot(3,1,2);
    plot(start_fr:seq_end_fr, collection_of_end_point_errs(:,1), 'r-', 'DisplayName', '2nd-order best fit');
    hold on;
    plot(start_fr:seq_end_fr, collection_of_end_point_errs(:,2), 'g-', 'DisplayName', '3rd-order best fit');
    xlabel('frame index');
    ylabel({'end-point','error (m)'});
    xlim([start_fr, seq_end_fr]);
    legend;

    subplot(3,1,3)
    plot(start_fr:seq_end_fr, collection_of_err_percentage(:,1), 'r-', 'DisplayName', '2nd-order best fit');
    hold on;
    plot(start_fr:seq_end_fr, collection_of_err_percentage(:,2), 'g-', 'DisplayName', '3rd-order best fit');
    xlabel('frame index');
    ylabel({'error','percentage','(%)'});
    xlim([start_fr, seq_end_fr]);
    legend;

    set(gcf,'color','w');
end

%> draw a histogram of all average errors for both dynamic models
if drawAvgErrHistogram && ~run_segment_wise
    title_string = strcat(datasetName, {' '}, sequenceName);
    nBins = 40;
    figure;
    h1 = histogram(collection_of_avg_errs(:,1), nBins);
    hold on;
    h2 = histogram(collection_of_avg_errs(:,2), nBins);
    h2.BinWidth = h1.BinWidth;
    xlabel("average error in a 10-frame local segment");
    ylabel("number of 10-frame local segments");
    title(title_string);
    legend('2nd-order','3rd-order');
    set(gcf,'color','w');
end



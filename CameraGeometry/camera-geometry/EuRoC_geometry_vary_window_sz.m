%% Use gradient descent optimization to fit a EuRoC dataset trajectory from a helix model
% -- varying window size --

clc; clear all; close all;

% -- define parameters --
datasetName = 'EuRoC';
sequenceName = 'MH_03_medium/mav0/';

% -- read all R, T, and time stamp of each frame from ground truth dataset --
[all_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, '');

start_fr = 1905;
end_fr = size(all_T,2)-50;
window_size = 5;

% -- dynamic model params and geometry model params --
dynamicParams = [];
geometryParams = [];

% ==================== CURVE GEOEMETRY MODEL FITTING =================
[dynamicParams, geometryParams, r, t] = self_dynamics_geometry_params_extractor(start_fr, start_fr+window_size, all_T, time_per_fr, '3d', dynamicParams, geometryParams);

vel = dynamicParams.vel_vec(3,:)';
acc = dynamicParams.acc_vec(3,:)';

% -- calculate the initial frenet frame --
vec_acc_cross_prod = cross(vel, acc);
init_binormal = vec_acc_cross_prod / norm(vec_acc_cross_prod);
init_tangent = vel / norm(vel);
tangent_binormal_cross_prod = -cross(init_tangent, init_binormal);
init_normal = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);

init_FrenetFrame.T = init_tangent;
init_FrenetFrame.N = init_normal;
init_FrenetFrame.B = init_binormal;

% -- valid trajectory points and scale up --
trajectory_pts = r(:,3:end);

% -- compute arcLength --
cumulative_arcLength = 0;
arcLength = zeros(size(trajectory_pts, 2),1);
for i = 2:size(trajectory_pts, 2)
    cumulative_arcLength = cumulative_arcLength + norm(trajectory_pts(:,i)-trajectory_pts(:,i-1));
    arcLength(i,1) = cumulative_arcLength;
end

% -- initial guesses --
init_guess.T = init_FrenetFrame.T;
init_guess.N = init_FrenetFrame.N;
init_guess.B = init_FrenetFrame.B;
init_guess.curvature = 20;
init_guess.torsion = -1;

% -- Numerical Optimization: GRADIENT DESCENT!!!! --
% -- [Approach 1] optimize only curvature, torsion, and normal vector --
GD_params.lr = 0.01;
GD_params.iter_nums = 10000;
GD_params.first_success_counts = 200;
GD_params.second_success_counts = 300;

[updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_curve_fit(init_guess, GD_params, trajectory_pts, arcLength);

opt_params.k = updated_geometry.k;
opt_params.tau = updated_geometry.tau;
opt_params.T = updated_geometry.T;
opt_params.N = updated_geometry.N;
opt_params.B = updated_geometry.B;

fprintf("\n");
fprintf("==== [Approach 1] optimize curvature and torsion ====\n");
fprintf("optimal curvature: %f\n", opt_params.k);
fprintf("optimal torsion: %f\n", opt_params.tau);
fprintf("average point error: %f\n", collect_err(end,1));
fprintf("error percentage: %f\n", err_pctge);
fprintf("\n");

% -- create the helix curve using optimized geoemtry parameters --
[curveFromModel, ~] = self_generateHelixFromModel(opt_params.k, opt_params.tau, opt_params.T, opt_params.N, opt_params.B, arcLength', trajectory_pts(:,1));

% -- store necessary data for drawings --
draw_A1.k = opt_params.k;
draw_A1.tau = opt_params.tau;
draw_A1.T = opt_params.T;
draw_A1.N = opt_params.N;
draw_A1.B = opt_params.B;
draw_A1.curve = curveFromModel;

%% ------------------------------------------------------------------------
% -- [Approach 2] Optimize curvature, torsion, and initial Frenet frame --
% -- i) convert initial frenet frame to unit sphere angles --
[thetaT, phiT, theta_mkortho, phi_mkortho] = self_cvt_frenetframe_to_unit_sphere_angle(init_guess.T, init_guess.N);
% -- ii) initial guesses --
init_guess.theta1 = thetaT;
init_guess.phi1 = phiT;
init_guess.theta2 = theta_mkortho;
init_guess.phi2 = phi_mkortho;
% -- iii) gradient descent parameters --
GD_params.lr = 0.05;
GD_params.lr_frenetframe = 0.01;
GD_params.iter_nums = 10000;
GD_params.first_success_counts = 100;
GD_params.second_success_counts = 300;
% -- iv) do gradient descent --
fprintf('Doing Gradient Descent ...\n');
[updated_geometry, collect_lr, collect_err, err_pctge] = do_gradient_EuRoC_local_geometry_frenetframe_curve_fit(init_guess, GD_params, trajectory_pts, arcLength);
% -- v) get output updates --
opt_params.k = updated_geometry.k;
opt_params.tau = updated_geometry.tau;
opt_thetaT = updated_geometry.theta1;
opt_phiT = updated_geometry.phi1;
opt_theta_mkortho = updated_geometry.theta2;
opt_phi_mkortho = updated_geometry.phi2;
% -- vi) convert unit sphere angles to frenet frame --
[prop_T, prop_N, prop_B] = self_cvt_unit_sphere_angle_to_frenetframe(opt_thetaT, opt_phiT, opt_theta_mkortho, opt_phi_mkortho);

opt_params.T = prop_T;
opt_params.N = prop_N;
opt_params.B = prop_B;

fprintf("\n");
fprintf("==== [Approach 2] optimize curvature, torsion, and initial Frenet frame ====\n");
fprintf("optimal curvature: %f\n", opt_params.k);
fprintf("optimal torsion: %f\n", opt_params.tau);
fprintf("average point error: %f\n", collect_err(end,1));
fprintf("error percentage: %f\n", err_pctge);
fprintf("\n");

% -- create the helix curve using optimized geoemtry parameters --
[curveFromModel, ~] = self_generateHelixFromModel(opt_params.k, opt_params.tau, opt_params.T, opt_params.N, opt_params.B, arcLength', trajectory_pts(:,1));

% -- store necessary data for drawings --
draw_A2.k = opt_params.k;
draw_A2.tau = opt_params.tau;
draw_A2.T = opt_params.T;
draw_A2.N = opt_params.N;
draw_A2.B = opt_params.B;
draw_A2.curve = curveFromModel;
    


% =============== DRAWINGS ===========
% -- PLOTTINGS --
% -- deisplay curves --
r1 = trajectory_pts';
r2 = draw_A1.curve';
r3 = draw_A2.curve';
figure;
plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'ground truth');
hold on;
plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', '[Approach 1]');
hold on;
plot3(r3(:,1), r3(:,2), r3(:,3), 'go-', 'DisplayName', '[Approach 2]');
hold on;

scale_factor = 0.005;

% -- draw original initial frenet frame at the initial point --
drawTangent = (draw_A1.T).*scale_factor;
drawNormal = (draw_A1.N).*scale_factor;
drawBinormal = (draw_A1.B).*scale_factor;
endpoint_T = r1(1,:) + drawTangent';
endpoint_N = r1(1,:) + drawNormal';
endpoint_B = r1(1,:) + drawBinormal';

line([r1(1,1), endpoint_T(1,1)],[r1(1,2), endpoint_T(1,2)],[r1(1,3), endpoint_T(1,3)],'Color','k','LineStyle','-', 'DisplayName', 'original T');
hold on;
line([r1(1,1), endpoint_N(1,1)],[r1(1,2), endpoint_N(1,2)],[r1(1,3), endpoint_N(1,3)],'Color','c','LineStyle','-', 'DisplayName', 'original N');
hold on;
line([r1(1,1), endpoint_B(1,1)],[r1(1,2), endpoint_B(1,2)],[r1(1,3), endpoint_B(1,3)],'Color','m','LineStyle','-', 'DisplayName', 'original B');
hold on;

% -- draw optimized initial frenet frame at the initial point --
drawOptTangent = (draw_A2.T).*scale_factor;
drawOptNormal = (draw_A2.N).*scale_factor;
drawOptBinormal = (draw_A2.B).*scale_factor;
endpoint_OptT = r1(1,:) + drawOptTangent';
endpoint_OptN = r1(1,:) + drawOptNormal';
endpoint_OptB = r1(1,:) + drawOptBinormal';

line([r1(1,1), endpoint_OptT(1,1)],[r1(1,2), endpoint_OptT(1,2)],[r1(1,3), endpoint_OptT(1,3)],'Color','k','LineStyle','--', 'DisplayName', 'optimized T');
hold on;
line([r1(1,1), endpoint_OptN(1,1)],[r1(1,2), endpoint_OptN(1,2)],[r1(1,3), endpoint_OptN(1,3)],'Color','c','LineStyle','--', 'DisplayName', 'optimized N');
hold on;
line([r1(1,1), endpoint_OptB(1,1)],[r1(1,2), endpoint_OptB(1,2)],[r1(1,3), endpoint_OptB(1,3)],'Color','m','LineStyle','--', 'DisplayName', 'optimized B');
hold on;

legend;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
set(gcf,'color','w');

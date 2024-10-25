%% Use optimization to fit a curve by bi-helix, given a curve as two connected helix
% i.e. Create two connected helix, and try to fit this curve by bi-helix

clc; clear all; close all;

% =============================================
% -- This part creates an explicit circular helix used for verify the
% model-generated helix. Should be unused for generating double connected 
% helix from model --
% % -- create circular helix in terms of arc length s --
% s = 0:3:15;
% % -- define parameters --
% a1 = 10;
% b1 = 7;
% add_noise = 0;
% [r1, dynamicParams, geometryParams] = self_create_explicit_helix(s1, a1, b1, add_noise, 1);
% 
% curvature = geometryParams.curvature(1,1);
% torsion = geometryParams.torsion(1,1);
% 
% % -- Calculate initial Frenet Frame --
% vel_vector = dynamicParams.vel_vec(1,:)';
% acc_vector = dynamicParams.acc_vec(1,:)';
% vec_acc_cross_prod = cross(vel_vector, acc_vector);
% binormal = vec_acc_cross_prod / norm(vec_acc_cross_prod);
% tangent = vel_vector / norm(vel_vector);
% tangent_binormal_cross_prod = -cross(tangent, binormal);
% normal = tangent_binormal_cross_prod / norm(tangent_binormal_cross_prod);
% % -- initial helix curve point --
% init_C = r1(1,:)';
% ===============================================
s1 = 0:0.3:1.5;
curvature1 = 0.6;
torsion1 = 0.01;
tangent = [1;0;0];
normal = [0;1;0];
binormal = [0;0;1];

% -- initial Frenet frame and start curve point --
init_T = tangent;
init_N = normal;
init_B = binormal;
init_C = [0;0;0];

% -- create circular helix using peicewise constant model --
[C1, last_frenet_frame] = self_create_propogated_Frenet_helix(s1, curvature1, torsion1, init_T, init_N, init_B, init_C);

% -- for the continual second model --
s2 = 0:0.2:1;
curvature2 = 3;
torsion2 = 0.3;
init_T = last_frenet_frame.tangent;
init_N = last_frenet_frame.normal;
init_B = last_frenet_frame.binormal;
init_C = C1(:,end);
[C2, ~] = self_create_propogated_Frenet_helix(s2, curvature2, torsion2, init_T, init_N, init_B, init_C);

% r1 = C1';
% r2 = C2';
% figure;
% plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'first helix');
% hold on;
% plot3(r2(:,1), r2(:,2), r2(:,3), 'go-', 'DisplayName', 'second helix');
% legend;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis equal;
% set(gcf,'color','w');

% -- concatenate two circular helix --
full_C = [C1, C2(:,2:end)];

% -- approximate the first helix by looping over all intermediate points as
% junction points, and pick the one with the smallest average point
% pairwise error --

num_of_pts = size(s1,2) + size(s2,2) - 1;

% -- TEST ON FRIDAY, THE FULL BIHELIX CURVE MODEL!!!! (UPDATE: PASSED!)--
junction_pt_idx = 9;
% [fullCurveModel, arcLength] = self_generateBihelixFromModel(curvature1, torsion1, curvature2, torsion2, tangent, normal, binormal, full_C, junction_pt_idx);

% -- apply gradient descent to find optimal curvature and torsion --
% 1) define parameters
gdParams.lr = 0.001;
gdParams.num_of_iter = 5000;
gdParams.num_of_first_success = 500;
gdParams.num_of_second_success = 600;

initParams.init_k0 = 0.25;
initParams.init_tau0 = 0.01;
initParams.init_k1 = 0.9;
initParams.init_tau1 = 0.01;

initFrenetFrame.T = tangent;
initFrenetFrame.N = normal;
initFrenetFrame.B = binormal;
[opt_res, collect_gd_errors, collect_learning_rate] = self_GD_min_curve_err(gdParams, initParams, initFrenetFrame, full_C, junction_pt_idx);

fprintf("optimal curvature #0: %f\n", opt_res.k0);
fprintf("optimal torsion #0: %f\n", opt_res.tau0);
fprintf("optimal curvature #1: %f\n", opt_res.k1);
fprintf("optimal torsion #1: %f\n", opt_res.tau1);
fprintf("average point error: %f\n", collect_gd_errors(end,1));

% -- generate the bihelix curve --
[fullCurveModel, arcLength] = self_generateBihelixFromModel(opt_res.k0, opt_res.tau0, opt_res.k1, opt_res.tau1, tangent, normal, binormal, full_C, junction_pt_idx);


% [curveFromModel_first, ~] = self_generateCurveFromModel(opt_res.curvature, opt_res.torsion, opt_res.T, opt_res.N, opt_res.B, trajectory_pts, trajectory_pts(:,1), 'peicewise_constant_model');
% [curveFromModel_second, ~] = self_generateCurveFromModel(opt_res_second.curvature, opt_res_second.torsion, opt_res_second.T, opt_res_second.N, opt_res_second.B, trajectory_pts_second, trajectory_pts_second(:,1), 'peicewise_constant_model');

r1 = C1';
r2 = C2';
r3 = fullCurveModel';
figure;
plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'first helix');
hold on;
plot3(r2(:,1), r2(:,2), r2(:,3), 'go-', 'DisplayName', 'second helix');
hold on;
plot3(r3(:,1), r3(:,2), r3(:,3), 'ro-', 'DisplayName', 'bihelical model');
legend;
title_str = strcat("Connecting at point ", string(junction_pt_idx));
title(title_str);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
set(gcf,'color','w');

%% -- take a test! --
% total_s = [s1 s1(1,end)+s2(2:end)];
% pt_idx = 7;
% test_teajectory = full_C(:,1:pt_idx);
% x = full_C(1,pt_idx);
% y = full_C(2,pt_idx);
% z = full_C(3,pt_idx);
% L = total_s(1,pt_idx);
% alpha1 = z / (L-x);
% R = (y^2*(1+alpha1^2) + (x-alpha1*z)^2)/(2*y*(1+alpha1^2));
% test_k = 1/(R*(1+alpha1^2));
% test_tau = alpha1 / (R*(1+alpha1^2));
% fprintf('\n');
% fprintf("test curvature: %f\n", test_k);
% fprintf("test torsion: %f\n", test_tau);
% [test_curveFromModel_second, ~] = self_generateCurveFromModel(test_k, test_tau, tangent, normal, binormal, test_teajectory, test_teajectory(:,1), 'peicewise_constant_model');

%% -- take a test for analytic solution? --
% syms x y z
% total_s = [s1 s1(1,end)+s2(2:end)];
% pt_idx = 6;
% test_teajectory = full_C(:,1:pt_idx);
% L = total_s(1,pt_idx);
% alpha1 = z / (L-x);
% R = (y^2*(1+alpha1^2) + (x-alpha1*z)^2)/(2*y*(1+alpha1^2));
% k = 1/(R*(1+alpha1^2));
% tau = alpha1 / (R*(1+alpha1^2));
% omega = sqrt(k^2+tau^2);
% 
% E1 = ((alpha1^2+1)*x - alpha1^2*L)*omega == sin(omega*L);
% E2 = (alpha1*L - z*(alpha1^2+1))*omega == sin(omega*L)*alpha1;
% E3 = (x-alpha1*z)^2-y^2*(1+alpha1^2) == cos(omega*L)*((x+alpha1*z)^2+y^2*(1+alpha1^2));
% % [solx,soly,solz] = vpasolve([E1 E2 E3], [x y z], [1.0; 0.5; 0]);

%% -- take another test!! --
% MH_01_easy, frame 185~195
test_trajectory2 = [-0.2444, -0.2486, -0.2494, -0.2471, -0.2417, -0.2332, -0.2216, -0.2076, -0.1915, -0.1738; ...
   -0.0129, -0.0134, -0.0138, -0.0143, -0.0147, -0.0151, -0.0160, -0.0172, -0.0187, -0.0201; ...
   -0.1047, -0.1060, -0.1073, -0.1088, -0.1103, -0.1120, -0.1134, -0.1151, -0.1165, -0.1176];
test_trajectory2 = test_trajectory2*10;
init_T = [-0.9781; -0.0359; -0.2051];
init_N = [0.2081; -0.1932; -0.9588];
init_B = [-0.0052; -0.9805; 0.1964];

% MH_03_medium, frame 1905~1915
% test_trajectory2 = [0.8998, 0.8986, 0.8970, 0.8953, 0.8933, 0.8908, 0.8877, 0.8840, 0.8796, 0.8745; ...
%     0.1178, 0.1233, 0.1271, 0.1293, 0.1300, 0.1295, 0.1275, 0.1244, 0.1200, 0.1142; ...
%    -1.0822, -1.0751, -1.0693, -1.0649, -1.0619, -1.0600, -1.0596, -1.0608, -1.0635, -1.0671];
% test_trajectory2 = test_trajectory2*10;
% init_T = [-0.1044; 0.5914; 0.7996];
% init_N = [-0.7438; -0.5801; 0.3320];
% init_B = [0.6602; -0.5600; 0.5004];


gdParams.lr = 0.01;
gdParams.num_of_iter = 20000;
gdParams.num_of_first_success = 150;
gdParams.num_of_second_success = 300;

initParams.init_k0 = 35;
initParams.init_tau0 = 3;
initParams.init_k1 = 2;
initParams.init_tau1 = 1;

initFrenetFrame.T = init_T;
initFrenetFrame.N = init_N;
initFrenetFrame.B = init_B;
j_pts = 4;
[opt_res, collect_gd_errors, collect_learning_rate, err_percentage] = self_GD_min_curve_err(gdParams, initParams, initFrenetFrame, test_trajectory2, j_pts);

% error percentage
%err_percentage = (collect_gd_errors(end,1) / arcLength(end,1))*100;

fprintf('\n');
fprintf("optimal curvature #0: %f\n", opt_res.k0);
fprintf("optimal torsion #0: %f\n", opt_res.tau0);
fprintf("optimal curvature #1: %f\n", opt_res.k1);
fprintf("optimal torsion #1: %f\n", opt_res.tau1);
fprintf("average point error: %f\n", collect_gd_errors(end,1));
fprintf("error percentage: %f\n", err_percentage);


[fullCurveModel, arcLength] = self_generateBihelixFromModel(opt_res.k0, opt_res.tau0, opt_res.k1, opt_res.tau1, opt_res.T, opt_res.N, opt_res.B, test_trajectory2, j_pts);
[initCurveModel, ~] = self_generateBihelixFromModel(initParams.init_k0, initParams.init_tau0, initParams.init_k1, initParams.init_tau1, init_T, init_N, init_B, test_trajectory2, j_pts);

r1 = test_trajectory2';
r2 = fullCurveModel';
r3 = initCurveModel';
figure;
plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'ground truth');
hold on;
plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', 'optimal curve');
%hold on;
%plot3(r3(:,1), r3(:,2), r3(:,3), 'ko-', 'DisplayName', 'initial guess');
legend;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
set(gcf,'color','w');

%% -- plot approximated curve and true double helix --
% figure;
% r = full_C';
% plot3(r(:,1), r(:,2), r(:,3), 'bo-', 'DisplayName', 'true double helix');
% hold on;
% r2 = curveFromModel_first';
% plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', 'approximated helix 1');
% hold on;
% % r3 = curveFromModel_second';
% % plot3(r3(:,1), r3(:,2), r3(:,3), 'go-', 'DisplayName', 'approximated helix 2');
% % hold on;
% r4 = test_curveFromModel_second';
% plot3(r4(:,1), r4(:,2), r4(:,3), 'mo-', 'DisplayName', 'test approximate helix');
% legend;
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis equal;
% set(gcf,'color','w');

%% -- plot the circular helix trajectory --
r1 = full_C';
r2 = fullCurveModel';
figure;
plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'true');
hold on;
plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', 'model');
legend;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
set(gcf,'color','w');

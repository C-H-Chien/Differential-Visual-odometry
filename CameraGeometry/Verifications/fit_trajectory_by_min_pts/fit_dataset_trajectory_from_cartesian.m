%% Use gradient descent optimization to fit a real-dataset trajectory from a caetesian parametrized circular helix curve

clc; clear all; close all;

% MH_01_easy, frame 185~195
% full_C = [-0.2444, -0.2486, -0.2494, -0.2471, -0.2417, -0.2332, -0.2216, -0.2076, -0.1915, -0.1738; ...
%           -0.0129, -0.0134, -0.0138, -0.0143, -0.0147, -0.0151, -0.0160, -0.0172, -0.0187, -0.0201; ...
%           -0.1047, -0.1060, -0.1073, -0.1088, -0.1103, -0.1120, -0.1134, -0.1151, -0.1165, -0.1176];
% full_C = full_C*10;
% init_T = [-0.9781; -0.0359; -0.2051];
% init_N = [0.2081; -0.1932; -0.9588];
% init_B = [-0.0052; -0.9805; 0.1964];

% MH_03_medium, frame 1905~1915
full_C = [0.8998, 0.8986, 0.8970, 0.8953, 0.8933, 0.8908, 0.8877, 0.8840, 0.8796, 0.8745; ...
          0.1178, 0.1233, 0.1271, 0.1293, 0.1300, 0.1295, 0.1275, 0.1244, 0.1200, 0.1142; ...
          -1.0822, -1.0751, -1.0693, -1.0649, -1.0619, -1.0600, -1.0596, -1.0608, -1.0635, -1.0671];
full_C = full_C*10;
init_T = [-0.1044; 0.5914; 0.7996];
init_N = [-0.7438; -0.5801; 0.3320];
init_B = [0.6602; -0.5600; 0.5004];
start_fr = 1905;
seq_end_fr = 1915;

% compute arcLength
cumulative_arcLength = 0;
arcLength = zeros(1, size(full_C, 2));
for i = 2:size(full_C, 2)
    cumulative_arcLength = cumulative_arcLength + norm(full_C(:,i)-full_C(:,i-1));
    arcLength(1,i) = cumulative_arcLength;
end


%% -- do gradient descent, point by point --
tangent = [1;0;0];
normal = [0;1;0];
binormal = [0;0;1];

% gradient descent parameters
collect_opt_curve_pts = zeros(size(arcLength(1,:), 2)-1, 3);
collect_opt_geometry = zeros(size(arcLength(1,:), 2)-1, 2);
init_pt = full_C(:,1);
last_pt = init_pt;
init_Frenetframe.t = tangent;
init_Frenetframe.n = normal;
init_Frenetframe.b = binormal;
propogated_Frenetframe.t = init_T;
propogated_Frenetframe.n = init_N;
propogated_Frenetframe.b = init_B;

collect_tangent = zeros(size(full_C, 2), 3);
collect_normal = zeros(size(full_C, 2), 3);
collect_binormal = zeros(size(full_C, 2), 3);
collect_tangent(1,:) = init_T';
collect_normal(1,:) = init_N';
collect_binormal(1,:) = init_B';
collect_opt_curve_pts(1,:) = full_C(:,1);
for i = 2:size(full_C, 2)
%for i = 2:4
    % arc length
    sc = arcLength(1,i)-arcLength(1,i-1);
    
    % gradient descent parameters
    lr = 0.5;
    gd_iter_num = 20000;
    first_success_count = 100;
    success_count = 0;
    collect_err = zeros(gd_iter_num,1);
    collect_lr = zeros(gd_iter_num,1);
    
    % move the start curve point to origin and Frenet frame to cartesian basis
    transl = last_pt - [0;0;0];
    rotm = [propogated_Frenetframe.t'; propogated_Frenetframe.n'; propogated_Frenetframe.b' ];
    
    % transform the ground truth to the frenet frame coordinate
    transform_GT = rotm*(full_C(:,i) - transl);
    
    k_guess = 3;
    tau_guess = 0.1;
    
    % generate the helix and create initial guess of curve point
    [helixFromModel, ~] = self_generateHelixFromModel(k_guess, tau_guess, init_Frenetframe.t, init_Frenetframe.n, init_Frenetframe.b, sc, [0;0;0]);
    cand_x = helixFromModel(1);
    cand_y = helixFromModel(2);
    cand_z = helixFromModel(3);
    
    for gd_iter = 1:gd_iter_num
        
        % calculate error before updating
        before_err = sqrt((cand_x - transform_GT(1,1))^2+(cand_y - transform_GT(2,1))^2+(cand_z - transform_GT(3,1))^2);
        
        % -- get gradients of curvature and torsion --
        Df = self_getGradientsOfCurvePts(cand_x,cand_y,cand_z,sc);
%         gradx = 2*(cand_x - transform_GT(1,1))*Df.f1x + 2*(cand_y - transform_GT(2,1))*Df.f2x + 2*(cand_z - transform_GT(3,1))*Df.f3x;
%         grady = 2*(cand_x - transform_GT(1,1))*Df.f1y + 2*(cand_y - transform_GT(2,1))*Df.f2y + 2*(cand_z - transform_GT(3,1))*Df.f3y;
%         gradz = 2*(cand_x - transform_GT(1,1))*Df.f1z + 2*(cand_y - transform_GT(2,1))*Df.f2z + 2*(cand_z - transform_GT(3,1))*Df.f3z;
        
        gradx = 2*(cand_x - transform_GT(1,1))*Df.f1x;
        grady = 2*(cand_y - transform_GT(2,1))*Df.f2y;
        gradz = 2*(cand_z - transform_GT(3,1))*Df.f3z;
        
        % -- update the curve point --
        upd_x = cand_x - lr * gradx;
        upd_y = cand_y - lr * grady;
        upd_z = cand_z - lr * gradz;
        
        % calculate error
        err = sqrt((upd_x - transform_GT(1,1))^2+(upd_y - transform_GT(2,1))^2+(upd_z - transform_GT(3,1))^2);
        collect_err(gd_iter, 1) = err;
        collect_lr(gd_iter, 1) = lr;
        
        % stop the gradient descent if the error is too small --
        if err < 0.0005
            % get the updated curve point (x,y,z)
            cand_x = upd_x;
            cand_y = upd_y;
            cand_z = upd_z;
            break;
        end
        
        % adaptivelt adjust the learning rate
        if (err > before_err)
            lr = lr * 0.5;
            success_count = 0;
            continue;
        else
            % update the candidate curve pt (x,y,z)
            cand_x = upd_x;
            cand_y = upd_y;
            cand_z = upd_z;
            % increase the learning rate if the err is still large and there are many successes
            if (success_count >= first_success_count)
                if err > 0.5
                    lr = lr * 3.0;
                elseif err < 0.005
                    lr = lr * 1.5;
                else
                    lr = lr * 2.0;
                end
                success_count = 0;
            end
            success_count = success_count + 1;
        end
    end
    
    % record the optimal geometry parameters of helix
    alpha1 = cand_z / (sc-cand_x);
    R = (cand_y^2*(1+alpha1^2) + (cand_x-alpha1*cand_z)^2)/(2*cand_y*(1+alpha1^2));
    k = 1/(R*(1+alpha1^2));
    tau = alpha1 / (R*(1+alpha1^2));
    collect_opt_geometry(i-1,1) = k;
    collect_opt_geometry(i-1,2) = tau;
    opt_pt = [cand_x; cand_y; cand_z];
    
    % transform back to the original coordinate
    collect_opt_curve_pts(i,:) = inv(rotm)*opt_pt + transl;
    
    % propogate the frenet frame
    [~, FrenetFrame] = self_generateHelixFromModel(k, tau, propogated_Frenetframe.t, propogated_Frenetframe.n, propogated_Frenetframe.b, sc, last_pt);
    propogated_Frenetframe.t = FrenetFrame.T;
    propogated_Frenetframe.n = FrenetFrame.N;
    propogated_Frenetframe.b = FrenetFrame.B;
    
    collect_tangent(i,:) = propogated_Frenetframe.t';
    collect_normal(i,:) = propogated_Frenetframe.n';
    collect_binormal(i,:) = propogated_Frenetframe.b';
    
    % record the optimal curve point
    %collect_opt_curve_pts(i,:) = endCurvePt;

    last_pt = collect_opt_curve_pts(i,:)';
end

% compute average pairwise point error and err percentage
avg_err = 0;
for i = 2:size(full_C, 2)
    avg_err = avg_err + norm(full_C(:,i)-collect_opt_curve_pts(i,:)');
end
avg_err = avg_err / (size(full_C,2)-1);
err_percentage = (avg_err / arcLength(1,end))*100;

%% -- plot the circular helix trajectory --
%pt_idx = 1;
r1 = full_C';
%r2 = collect_opt_curve_pts(1:pt_idx,:);
r2 = collect_opt_curve_pts(1:end,:);
r1 = r1*0.1;
r2 = r2*0.1;
figure;
plot3(r1(:,1), r1(:,2), r1(:,3), 'bo-', 'DisplayName', 'ground truth');
hold on;
plot3(r2(:,1), r2(:,2), r2(:,3), 'ro-', 'DisplayName', 'optimal curve');
hold on;

% draw frenet frame
% -- scale down --
for pt_idx = 1:size(collect_tangent, 1)
%for pt_idx = 1:1
    drawTangent = collect_tangent(pt_idx,:);
    drawNormal = collect_normal(pt_idx,:);
    drawBinormal = collect_binormal(pt_idx,:);
    scale_factor = 0.003;
    drawTangent = drawTangent.*scale_factor;
    drawNormal = drawNormal.*scale_factor;
    drawBinormal = drawBinormal.*scale_factor;

    % -- draw frenet frame at the initial point --
    endpoint_T = r1(pt_idx,:) + drawTangent;
    endpoint_N = r1(pt_idx,:) + drawNormal;
    endpoint_B = r1(pt_idx,:) + drawBinormal;

    line([r1(pt_idx,1), endpoint_T(1,1)],[r1(pt_idx,2), endpoint_T(1,2)],[r1(pt_idx,3), endpoint_T(1,3)],'Color','g','LineStyle','-', 'DisplayName', 'tangent');
    hold on;
    line([r1(pt_idx,1), endpoint_N(1,1)],[r1(pt_idx,2), endpoint_N(1,2)],[r1(pt_idx,3), endpoint_N(1,3)],'Color','c','LineStyle','-', 'DisplayName', 'normal');
    hold on;
    line([r1(pt_idx,1), endpoint_B(1,1)],[r1(pt_idx,2), endpoint_B(1,2)],[r1(pt_idx,3), endpoint_B(1,3)],'Color','m','LineStyle','-', 'DisplayName', 'binormal');
    hold on;
end
legend;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
set(gcf,'color','w');

%% plot curvature and torsion of each segment
% figure;
% subplot(1,2,1);
% plot(start_fr:seq_end_fr, collect_opt_geometry(:,1), 'b-');
% xlabel('frame index');
% ylabel({'total','arc length','(m)'});
% xlim([start_fr, seq_end_fr]);
% set(gcf,'color','w');

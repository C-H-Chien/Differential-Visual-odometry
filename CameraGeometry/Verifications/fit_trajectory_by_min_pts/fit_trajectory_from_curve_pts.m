%% Use optimization to fit a curve by bi-helix, given a curve as two connected helix
% i.e. Create two connected helix, and try to fit this curve by bi-helix

clc; clear all; close all;

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
curvature2 = 1.2;
torsion2 = 0.01;
init_T = last_frenet_frame.tangent;
init_N = last_frenet_frame.normal;
init_B = last_frenet_frame.binormal;
init_C = C1(:,end);
[C2, ~] = self_create_propogated_Frenet_helix(s2, curvature2, torsion2, init_T, init_N, init_B, init_C);

% -- concatenate two circular helix --
full_C = [C1, C2(:,2:end)];
arcLength = [s1, s1(1,end)+s2(1,2:end)];


%% -- do gradient descent, point by point --
% initial guess of curvature and torsion
k_guess = 0.1;
tau_guess = 0.1;

% gradient descent parameters
collect_opt_curve_pts = zeros(size(arcLength(1,:), 2)-1, 3);
collect_opt_geometry = zeros(size(arcLength(1,:), 2)-1, 2);
init_pt = [0;0;0];
last_pt = init_pt;
init_Frenetframe.t = tangent;
init_Frenetframe.n = normal;
init_Frenetframe.b = binormal;
propogated_Frenetframe.t = tangent;
propogated_Frenetframe.n = normal;
propogated_Frenetframe.b = binormal;
for i = 2:size(full_C, 2)
%for i = 2:8
    % gradient descent parameters
    lr = 0.5;
    gd_iter_num = 3000;
    first_success_count = 50;
    success_count = 0;
    collect_err = zeros(gd_iter_num,1);
    collect_lr = zeros(gd_iter_num,1);
    
    % move the start curve point to origin and Frenet frame to cartesian basis
    transl = full_C(:,i-1) - [0;0;0];
    rotm = [propogated_Frenetframe.t'; propogated_Frenetframe.n'; propogated_Frenetframe.b' ];
    
    % transform the ground truth to the frenet frame coordinate
    transform_GT = rotm*(full_C(:,i) - transl);
    
    % generate the helix and create initial guess of curve point
    sc = arcLength(1,i)-arcLength(1,i-1);
    [helixFromModel, ~] = self_generateHelixFromModel(k_guess, tau_guess, init_Frenetframe.t, init_Frenetframe.n, init_Frenetframe.b, sc, [0;0;0]);
    cand_x = helixFromModel(1);
    cand_y = helixFromModel(2);
    cand_z = helixFromModel(3);
    
    for gd_iter = 1:gd_iter_num
        
        % calculate error before updating
        before_err = sqrt((cand_x - transform_GT(1,1))^2+(cand_y - transform_GT(2,1))^2+(cand_z - transform_GT(3,1))^2);
        
        % -- get gradients of curvature and torsion --
        Df = self_getGradientsOfCurvePts(cand_x,cand_y,cand_z,sc);
        gradx = 2*(cand_x - transform_GT(1,1))*Df.x;
        grady = 2*(cand_y - transform_GT(2,1))*Df.y;
        gradz = 2*(cand_z - transform_GT(3,1))*Df.z;
        
        % -- update the curve point --
        upd_x = cand_x - lr * gradx;
        upd_y = cand_y - lr * grady;
        upd_z = cand_z - lr * gradz;
        
        % calculate error
        err = sqrt((upd_x - transform_GT(1,1))^2+(upd_y - transform_GT(2,1))^2+(upd_z - transform_GT(3,1))^2);
        collect_err(gd_iter, 1) = err;
        collect_lr(gd_iter, 1) = lr;
        
        % stop the gradient descent if the error is too small --
        if err < 0.001
            % get the updated curve point (x,y,z)
            cand_x = upd_x;
            cand_y = upd_y;
            cand_z = upd_z;
            break;
        end
        
        % adaptivelt adjust the learning rate
        if (err > before_err)
            lr = lr * 0.01;
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
                    lr = lr * 2.0;
                elseif err < 0.1
                    lr = lr * 1.2;
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
    
    % record the optimal curve point
    %collect_opt_curve_pts(i,:) = endCurvePt;

    last_pt = collect_opt_curve_pts(i,:)';
end


%% -- plot the circular helix trajectory --
r1 = full_C';
r2 = collect_opt_curve_pts;
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

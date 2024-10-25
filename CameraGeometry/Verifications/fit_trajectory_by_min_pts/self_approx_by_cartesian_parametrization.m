%% make fit_dataset_trajectory_from_cartesian a function
function [collect_opt_curve_pts, collect_opt_geometry, collectFrenetFrame, avg_err] = self_approx_by_cartesian_parametrization(full_C, init_FrenetFrame, dataset)

if strcmp(dataset, 'EuRoC')
    full_C = full_C*10;
end

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
propogated_Frenetframe.t = init_FrenetFrame.T;
propogated_Frenetframe.n = init_FrenetFrame.N;
propogated_Frenetframe.b = init_FrenetFrame.B;

collect_tangent = zeros(size(full_C, 2), 3);
collect_normal = zeros(size(full_C, 2), 3);
collect_binormal = zeros(size(full_C, 2), 3);
collect_tangent(1,:) = init_FrenetFrame.T';
collect_normal(1,:) = init_FrenetFrame.N';
collect_binormal(1,:) = init_FrenetFrame.B';
collect_opt_curve_pts(1,:) = full_C(:,1);
for i = 2:size(full_C, 2)
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

    k_guess = 1;
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

    last_pt = collect_opt_curve_pts(i,:)';
end

collectFrenetFrame.collect_T = collect_tangent;
collectFrenetFrame.collect_N = collect_normal;
collectFrenetFrame.collect_B = collect_binormal;

% compute average pairwise point error
avg_err = 0;
for i = 2:size(full_C, 2)
    avg_err = avg_err + norm(full_C(:,i)-collect_opt_curve_pts(i,:)');
end
avg_err = avg_err / (size(full_C,2)-1);

if strcmp(dataset, 'EuRoC')
    collect_opt_curve_pts = collect_opt_curve_pts*0.1;
    %collect_opt_geometry(:,1) = collect_opt_geometry(:,1)*10;
end

end
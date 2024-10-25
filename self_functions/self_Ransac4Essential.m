function [finalE, inlierIndx] = self_Ransac4Essential(RANSACIterations, reprojectionErrorThreshold, matchImg1, matchImg2, K, debug_show)
    % Initially set maximal number of inliers to 0
    inlierNumMax = 0;
    
    % iterate a fixed number of times for RANSAC
    for i = 1 : RANSACIterations
        gamma1 = zeros(2, 5);
        gamma2 = zeros(2, 5);

        % Select 5 random matches
        idx = randperm(length(matchImg1), 5);
        for j = 1 : 5
            gamma1(1, j) = matchImg1(1, idx(j));
            gamma1(2, j) = matchImg1(2, idx(j));
            gamma2(1, j) = matchImg2(1, idx(j));
            gamma2(2, j) = matchImg2(2, idx(j));
        end

        % Construct the essential matrix E based on the 5 randomly selected points
        E = ComputeEssentialMatrix(gamma1, gamma2, K);

        % Take inverse operation of the intrinsic matrix
        K_inv = inv(K);
        
        % Compute coefficients of a line equation
        A = zeros(size(E, 3), length(matchImg1));
        B = zeros(size(E, 3), length(matchImg1));
        C = zeros(size(E, 3), length(matchImg1));
        for j = 1 : size(E, 3)
            calE = K_inv' * E{j} * K_inv;
            A(j, :) = calE(1, :) * matchImg1;
            B(j, :) = calE(2, :) * matchImg1;
            C(j, :) = calE(3, :) * matchImg1;
        end
        
        % Compute the distance from a point to a line for all matches
        dist = zeros(size(E, 3), length(matchImg1));
        denomOfDist = zeros(size(E, 3), length(matchImg1));
        numerOfDist = zeros(size(E, 3), length(matchImg1));
        A_ep = zeros(size(E, 3), length(matchImg1));
        B_it = zeros(size(E, 3), length(matchImg1));
        for k = 1 : length(matchImg1)
            A_ep(:,k) = A(:,k).*matchImg2(1, k);
            B_it(:,k) = B(:,k).*matchImg2(2, k);
        end
        numerOfDist = abs(A_ep + B_it + C);
        denomOfDist = A.^2 + B.^2;
        denomOfDist = sqrt(denomOfDist);
        dist = numerOfDist./denomOfDist;

        % Use Algebraic error for objective measurement
%         err = [];
%         for m = 1 : size(E,3)
%             e = diag(matchImg2' * K_inv' * E{m} * K_inv * matchImg1);
%             err(m,:) = e';
%         end
        
        % Investigate the maximal inliers among the essential matrix
        % candidates
        inlierIndx4all = [];
        for j = 1 : size(E, 3)
            inlierIndx4all = find(dist(j,:) < reprojectionErrorThreshold);
            %inlierIndx4all = find(abs(err(j,:)) < params.reprojectionErrorThreshold);
            NumOfInliers = size(inlierIndx4all, 2);
            if (NumOfInliers > inlierNumMax)
                inlierNumMax = NumOfInliers;
                finalE = E{j};
                inlierIndx = inlierIndx4all;
            end
            inlierIndx4all = [];
        end
    end
    
    if debug_show
        fprintf("Originally there are %d matches.\n", length(matchImg1));
        fprintf("Found %d inliers\n",inlierNumMax);
        disp("Essential matrix after RANSAC:");
        disp(finalE);
    end
end
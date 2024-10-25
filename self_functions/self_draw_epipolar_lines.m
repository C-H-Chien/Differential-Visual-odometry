function self_draw_epipolar_lines(im1, im2, matched_f1, matched_f2, F, inlier_Thresh, ...
                                  disp_Img_SideBySide, disp_Match_Lines, ...
                                  disp_Epipolar_Lines, disp_Num_Inliers, disp_With_Text)
    figure; clf;
    if disp_Img_SideBySide == 1
        imshow([im1, im2]);
    else
        imshow([im1; im2]);
    end
    hold on;
    axis image off;

    xa = matched_f1(1,:);
    xb = matched_f2(1,:);
    ya = matched_f1(2,:);
    yb = matched_f2(2,:);
    gamma1 = [xa; ya; ones(1,length(xa))];
    
    %> Compute epipolar line coefficients in pixels
    Apixel = F(1,:) * gamma1;
    Bpixel = F(2,:) * gamma1;
    Cpixel = F(3,:) * gamma1;

    %> Compute the distance from the features to the epipolar line
    A_xi  = Apixel.*xb;
    B_eta = Bpixel.*yb;
    numerOfDist = abs(A_xi + B_eta + Cpixel);
    denomOfDist = Apixel.^2 + Bpixel.^2;
    denomOfDist = sqrt(denomOfDist);
    dist21 = numerOfDist./denomOfDist;
    inlier_indices = find(dist21 > 0);
    % inlier_indices = find(dist21 <= inlier_Thresh);
    % inlier_indices = find(dist21 > inlier_Thresh);
    
    %> Pick inliers for visualization
    % indx = randperm(length(inlier_indices), disp_Num_Inliers);
    % disp_Inlier_Index = inlier_indices(indx);
    disp_Inlier_Index = inlier_indices;
%     disp_dist_Epipolar_Line = dist21(indx);

    sz = [length(disp_Inlier_Index) 3];
    disp_RGB_colors = unifrnd(0, 1, sz);
    
    %> Plot Inliers
    %> 1) Points
    if disp_Img_SideBySide == 1
        for i = 1:length(disp_Inlier_Index)
            ii = disp_Inlier_Index(i);
            if disp_With_Text == 1
                text(xa(ii), ya(ii), string(dist21(ii)), 'Color', 'red', 'BackgroundColor', 'w');
                hold on;
            end
            plot(xa(ii), ya(ii), 's', 'Color', disp_RGB_colors(i,:), 'LineWidth', 2);
            hold on;
            plot(xb(ii)+size(im1,2), yb(ii), 's', 'Color', disp_RGB_colors(i,:), 'LineWidth', 2);
            hold on;
        end
    else
        for i = 1:length(disp_Inlier_Index)
            ii = disp_Inlier_Index(i);
            if disp_With_Text == 1
                text(xa(ii), ya(ii), string(dist21(ii)), 'Color', 'red', 'BackgroundColor', 'w');
                hold on;
            end
            plot(xa(ii), ya(ii), 's', 'Color', disp_RGB_colors(i,:), 'LineWidth', 2);
            hold on;
            plot(xb(ii), yb(ii)+size(im1,1), 's', 'Color', disp_RGB_colors(i,:), 'LineWidth', 2);
            hold on;
        end
    end
    
    %> 2) Matched Lines
    if disp_Match_Lines
        for i = 1:length(disp_Inlier_Index)
            ii = disp_Inlier_Index(i);
            if disp_Img_SideBySide == 1
                h1 = line([xa(ii) ; xb(ii) + size(im1,2) ], ...
                          [ya(ii) ; yb(ii)]);
                hold on;
                set(h1,'linewidth', 1, 'Color', disp_RGB_colors(i,:));
            else
                h1 = line([xa(ii) ; xb(ii) ], ...
                          [ya(ii) ; yb(ii) + size(im1,1)]);
                hold on;
                set(h1,'linewidth', 1, 'Color', disp_RGB_colors(i,:));
            end
        end
    end
    
    %> 3) Epipolar Lines
    if disp_Epipolar_Lines
        yMin = -Cpixel./Bpixel;
        yMax = (-Cpixel - Apixel*size(im1,2)) ./ Bpixel;
        disp_yMin = yMin(disp_Inlier_Index);
        disp_yMax = yMax(disp_Inlier_Index);
        if disp_Img_SideBySide == 1
            for i = 1:length(disp_yMin)
                h2 = line([size(im1,2), 2*size(im1,2)], [round(disp_yMin(i)), round(disp_yMax(i))]);
                set(h2, 'Color', disp_RGB_colors(i,:), 'LineWidth', 1);
                hold on;
            end
        else
            for i = 1:length(disp_yMin)
                h2 = line([1, size(im1,2)], [round(disp_yMin(i)) + size(im1,1), round(disp_yMax(i)) + size(im1,1)]);
                set(h2, 'Color', disp_RGB_colors(i,:), 'LineWidth', 1);
                hold on;
            end
        end
    end
    hold off;
    
    if disp_Img_SideBySide == 1
        xlim([1, size(im1,2) + size(im1,2)]);
        ylim([1, size(im1,1)]);
    else
        xlim([1, size(im1,2)]);
        ylim([1, size(im1,1) + size(im1,1)]);
    end
end
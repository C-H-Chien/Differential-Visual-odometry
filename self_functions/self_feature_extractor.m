function [featureDescriptions, featurePoints, rows, cols] = self_feature_extractor(imgName, featureType, i, start_fr, datasetName, K_val)

    % -- load image and convert it to gray scale --
    img_load = imread(imgName);
    [rows, cols, channels] = size(img_load);
    if channels > 1
        grayScaleImg = im2double(rgb2gray(img_load));
    else
        grayScaleImg = im2double(img_load);
    end
    [height, width, ~] = size(grayScaleImg);
    
    %imshow(grayScaleImg);
    if strcmp(datasetName, 'EuRoC')
        %%%%%%%%%%%%%%%%% camera Parameters for EuRoC from sensor.yaml %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        radialDistortion = [-0.28340811, 0.07395907];
        tangentialDistortion = [0.00019359, 1.76187114e-05]';
        imageSize = [height, width];
        cameraParams = cameraParameters("K",K_val,"RadialDistortion",radialDistortion, "TangentialDistortion",tangentialDistortion, "ImageSize", imageSize);

        grayScaleImg = undistortImage(grayScaleImg, cameraParams); 
    end

    %imshow(grayScaleImg);
    
    % -- extract Features --
    if strcmp(featureType, 'ORB')
        keyPointsObj = detectORBFeatures(grayScaleImg);
        [binaryDescriptions, featurePoints] = extractFeatures(grayScaleImg, keyPointsObj);
        featureDescriptions = binaryDescriptions.Features;
    elseif strcmp(featureType, 'SURF')
        keyPointsObj_prev = detectSURFFeatures(grayScaleImg);
        [featureDescriptions, featurePoints] = extractFeatures(grayScaleImg, keyPointsObj_prev);
    elseif strcmp(featureType, 'vlfeat-SIFT')
        grayScaleImg_single = single(grayScaleImg);
        %[featurePoints, featureDescriptions] = vl_sift(grayScaleImg_single, 'PeakThresh', vlfeat_params.SiftThreshPeak, 'edgethresh', vlfeat_params.SiftThreshEdge);
        [fpts, fdes] = vl_sift(grayScaleImg_single);
        featureDescriptions = fdes';
        featurePoints = fpts';
    elseif strcmp(featureType, 'TAPIR')
        %> TAPIR read feature points
        TAPIR_file_path = ['/Users/qiwuzhang/Desktop/Brown/LEMS/code/feature_points/TAPIR_feature_points_', num2str(start_fr-1), '.mat'];
        TAPIR_file_id = fopen(TAPIR_file_path);
        if TAPIR_file_id == -1
            error('TAPIR file is not opened successfully');
        end
        TAPIR = load(TAPIR_file_path).matrix_3d;
        TAPIR = TAPIR(:, i-start_fr+1, :);
        TAPIR_feature_points = squeeze(TAPIR(:, 1, :));
        featureDescriptions = [];  
        featurePoints = TAPIR_feature_points;  
     elseif strcmp(featureType, 'TEEATO')
         TEEATO_folder_path = '/Users/qiwuzhang/Desktop/Brown/LEMS/code/TEEATO_feature_points';
         TEEATO_file = [TEEATO_folder_path, '/id_', num2str(start_fr), '.mat'];
         disp(TEEATO_file)
         TEEATO_file_id = fopen(TEEATO_file);
         if TEEATO_file_id == -1
            error('TEEATO file is not opened successfully');
         end

         TEEATO_feature_points = load(TEEATO_file).filtered_data;
         TEEATO = TEEATO_feature_points(i-start_fr+1, :, :);
         TEEATO_feature_points = squeeze(TEEATO(1, :, :));
         featureDescriptions = [];  
         featurePoints = TEEATO_feature_points; 

    end

end
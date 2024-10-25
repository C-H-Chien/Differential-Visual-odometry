function [rgb_filename, depth_filename] = ...
         get_ICL_NUIM_rgb_depth_lists(PARAMS, is_noisy)
    
    if is_noisy == 1
        sequence_path = strcat(PARAMS.DATASET_PATH_NOISY, PARAMS.SEQUENCE_NAME);
    else
        sequence_path = strcat(PARAMS.DATASET_PATH, PARAMS.SEQUENCE_NAME);
    end

    images_Path = strcat(sequence_path, 'rgb/');
    depths_Path = strcat(sequence_path, 'depth/');
    
    imageList = dir(images_Path);
    imageList(1:2, :) = [];
    Num_Of_Images = length(imageList)-1;

    rgb_filename     = strings(Num_Of_Images, 1);
    depth_filename   = strings(Num_Of_Images, 1);
    for i = 1 : Num_Of_Images
        imgStringName   = strcat(images_Path, string(i), ".png");
        depthStringName = strcat(depths_Path, string(i), ".png");
        rgb_filename(i,1)   = imgStringName;
        depth_filename(i,1) = depthStringName;
    end
end






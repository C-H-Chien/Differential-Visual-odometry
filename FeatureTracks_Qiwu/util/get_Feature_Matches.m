function [matched_f1_location, matched_f2_location] = get_Feature_Matches(PARAMS, Gray_img1, Gray_img2)
   
    if strcmp(PARAMS.FEATURE_TYPE, "SIFT")
        [f1, d1] = vl_sift(single(Gray_img1));
        [f2, d2] = vl_sift(single(Gray_img2));
        [matches, scores] = vl_ubcmatch(d1, d2, 1.6);
    elseif strcmp(PARAMS.FEATURE_TYPE, "SURF")
        keyPointsObj_img1 = detectSURFFeatures(Gray_img1);
        keyPointsObj_img2 = detectSURFFeatures(Gray_img2);
        [d1, f1] = extractFeatures(Gray_img1, keyPointsObj_img1);
        [d2, f2] = extractFeatures(Gray_img2, keyPointsObj_img2);
        f1 = f1.Location';
        f2 = f2.Location';
        matches = matchFeatures(d1, d2)';
    elseif strcmp(PARAMS.FEATURE_TYPE, "ORB")
        keyPointsObj_img1 = detectORBFeatures(Gray_img1);
        keyPointsObj_img2 = detectORBFeatures(Gray_img2);
        [d1, f1] = extractFeatures(Gray_img1, keyPointsObj_img1);
        [d2, f2] = extractFeatures(Gray_img2, keyPointsObj_img2);
        f1 = f1.Location';
        f2 = f2.Location';
        matches = matchFeatures(d1, d2)';
    else
        error("Invalid feature type!");
    end
    

    %> Get the matched keypoint locations
    matched_f1_location = double(f1(1:2, matches(1, :))');
    matched_f2_location = double(f2(1:2, matches(2, :))');
end
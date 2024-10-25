function clothoid = self_generateClothoidFromModel(params, arcLength)

    %> start pt of the clothoid
    init_pt   = params.curve_start_pt;

    %> collection of colothoid points, based on the number of arcLength
    clothoid = [];
    clothoid = [clothoid, init_pt];

    %> loop over all arc lengths and collect clothoid points
    for i = 1:size(arcLength, 1)
        params.arcLength = arcLength(i,1);
        curve_pt = evaluate_clothoid3d_model(params);
        clothoid = [clothoid, [curve_pt.x; curve_pt.y; curve_pt.z]];
    end

end
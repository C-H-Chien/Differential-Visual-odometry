function [R21, T21] = get_Rel_Pose(PARAMS, img_id1, img_id2, K)
    
    sequence_path = strcat(PARAMS.DATASET_PATH, PARAMS.SEQUENCE_NAME);
    
    %> Get Ground Truth
    ft_filePath = strcat(sequence_path, 'groundtruth.txt');
    gt_file = importdata(ft_filePath);
    time_stamps = gt_file(:,1);
    T_gts = gt_file(:,2:4);
    R_gts = gt_file(:,5:8);
    pos1 = img_id1;
    pos2 = img_id2;
    C1 = T_gts(pos1,:);
    C2 = T_gts(pos2,:);
    R1 = R_gts(pos1,:);
    R2 = R_gts(pos2,:);
    R1 = transpose(quat2rotm(R1([4 1 2 3])));
    R2 = transpose(quat2rotm(R2([4 1 2 3])));
    R21 = R2 * inv(R1);
    T21 = R2 * (C1' - C2');
    % T12x = [0 -T12(3) T12(2);
    %         T12(3) 0  -T12(1);
    %         -T12(2) T12(1) 0];

    %> Compute the Fundamental Matrix
%     E12 = T12x * R12;
    %T12=C2'-R12*C1';
%     F12 = transpose(inv(K)) * E12 * inv(K);
end

% -- function of reading dataset --
% -- input: dataset name, sequence name --
% -- output: all sequence of T, R, timestamp --
function [all_T, all_R, time_per_fr, imgFileName] = readDataset(datasetName, sequenceName, imgCategory)
    if strcmp(datasetName, 'KITTI')
        sequence = sequenceName;
        % -- image data directories --
        imgFolder = '/gpfs/data/bkimia/zqiwu/KITTI_gray/sequences/';
        imgDir = strcat(imgFolder, sequence);
        
        % -- read real poses from KITTI ground truths --
        poseFolder = '/gpfs/data/bkimia/zqiwu/';
        dataset = 'KITTI-gt/poses/';
        postfix = '.txt';
        
        % -- read R, T of each frame --
        pose_FileName = strcat(poseFolder, dataset, sequence, postfix);
        [all_R, all_T] = read_groundtruths_kitti(pose_FileName);
        
        % -- read time stamp of each pose --
        timeFile = '/times.txt';
        times_FileName = fullfile(imgDir, timeFile);
        disp(times_FileName);
        timesFileRd = fopen(times_FileName, 'r');
        ldata = textscan(timesFileRd, '%s', 'CollectOutput', true);
        line = string(ldata{1});
        time_per_fr = str2double(line);
        time_per_fr = time_per_fr';
        
        % -- read all image directories as a set of strings (imgFileName) --
        imgFileName = strings(size(all_T,2), 1);
        for i = 0:size(all_T, 2)-1
            imgName = num2str(i,'%06.f.png');
            imgFileName(i+1, 1) = fullfile(imgFolder, sequence, imgCategory, imgName);
        end
        
    elseif strcmp(datasetName, 'EuRoC')
        sequence = sequenceName;
        
        % -- read associated files --
        dataFolder = '/oscar/data/bkimia/zqiwu/EuRoC/';
        img_gt_align_files = 'associate_img_gt.txt';
        align_full_dir = fullfile(dataFolder, sequence, img_gt_align_files);

        % -- extract timestamps, body frame ground truths, and image file names 
        gt_img = fopen(align_full_dir, 'r');

        % Read_Data = importdata(gt_img);
        ldata = textscan(gt_img, '%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f', 'CollectOutput', true );
        Time_Stamps = ldata{1,1};
        Img_File_Names = string(ldata{1,2});
        GT_Poses = ldata{1,3};

        %> Transformation of a camera with respect to the body frame
        %  In EuRoC dataset, I found this the same across all sequences.
        %  SO let's make it constant here.
        R_BC = [0.0148655429818, -0.999880929698, 0.00414029679422;
                0.999557249008, 0.0149672133247, 0.025715529948;
               -0.0257744366974, 0.00375618835797, 0.999660727178];
        T_BC = [-0.0216401454975; -0.064676986768; 0.00981073058949];
        Tr_BC = [R_BC, T_BC; 0, 0, 0, 1];
        inv_R_BC = inv(Tr_BC);

        %> Time stamps
        sample_time = (Time_Stamps(2,1) - Time_Stamps(1,1))*1e-9;   % -- seconds --
        time_per_fr = zeros(1, size(Time_Stamps,1)); 
        for i = 1:size(Time_Stamps,1), time_per_fr(1,i) = (i-1)*sample_time; end
        Num_Of_Poses = size(GT_Poses, 1);
        Tr_cam = zeros(4, 4, Num_Of_Poses);
        for p = 1:Num_Of_Poses
            T_g = GT_Poses(p, 2:4)';
            Q_g = GT_Poses(p, 5:8);
            R_g = quat2rotm(Q_g);
            Tr_gt = [R_g, T_g; 0, 0, 0, 1];
            Tr_cam(:,:,p) = inv(inv_R_BC*inv(Tr_gt));
        end

        all_T(:,:) = Tr_cam(1:3,4,:);
        all_R = Tr_cam(1:3,1:3,:);

        %> A list of image names
        imgFileName = strings(Num_Of_Poses, 1);
        for i = 1:size(all_T, 2)
            imgName = string(Img_File_Names(i));
            imgFileName(i, 1) = fullfile(dataFolder, sequence, '/cam0/data/', imgName);
        end
        
    elseif strcmp(datasetName, 'BPOD')
        sequence = sequenceName;
        % -- image data directories --
        imgFolder = '/home/chchien/datasets/BPOD/bpod_pedestrian_dataset/gt_415/';
        imgDir = strcat(imgFolder, sequence);
        
        % -- read real poses from KITTI ground truths --
        poseFolder = '/home/chchien/datasets/BPOD/bpod_pedestrian_dataset/';
        dataset = 'gt_415/';
        postfix = '.txt';
        
        % -- read R, T of each frame --
        pose_FileName = strcat(poseFolder, dataset, sequence, postfix); 
        gt_bpod = fopen(pose_FileName, 'r');
        ldata = textscan(gt_bpod, '%f\t%f\t%f', 'CollectOutput', true );
        time_pose_gt = ldata{1,1};
        
        time_per_fr = time_pose_gt(:,1);
        time_per_fr = time_per_fr.*0.000001;
        all_T = time_pose_gt(:,2:3);
        
        % -- TODO: requires BPOD images and rotation data --
        all_R = zeros(size(all_T));
        imgFileName = strings(size(time_per_fr));
        
%         img_filename = ldata{1,2};
%         
%         % -- read time stamp of each pose --
%         timeFile = '/times.txt';
%         times_FileName = fullfile(imgDir, timeFile);
%         timesFileRd = fopen(times_FileName, 'r');
%         ldata = textscan(timesFileRd, '%s', 'CollectOutput', true);
%         line = string(ldata{1});
%         time_per_fr = str2double(line);
%         time_per_fr = time_per_fr';
%         
%         % -- read all image directories as a set of strings (imgFileName) --
%         imgFileName = strings(size(all_T,2), 1);
%         for i = 0:size(all_T, 2)-1
%             imgName = num2str(i,'%06.f.png');
%             imgFileName(i+1, 1) = fullfile(imgFolder, sequence, imgCategory, imgName);
%         end
    end
end


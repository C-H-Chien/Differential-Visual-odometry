clear;
close all;
%% -- plot number of N-length tracks tested from the datasets --
fileDirectory = '/Users/qiwuzhang/Desktop/Brown/LEMS/code/';
sequenceName = 'KITTI-seq05/';
fileName = 'num-of-tracks-per-length-';
collection_featureType = ["SURF"; "ORB"; "SIFT"];
filePostFix = '.txt';
zoomIn = 0;

for typeIdx = 1:size(collection_featureType, 1)
    dataRead_directory_string = strcat(fileDirectory, sequenceName, fileName, collection_featureType(typeIdx), filePostFix);
    disp(dataRead_directory_string);
    dataFileRead = fopen(dataRead_directory_string, 'r');

    readType_seq_str = '';
    for i = 1:24
        if i < 24
            readType_seq_str = strcat(readType_seq_str, '%d\t');
        else
            readType_seq_str = strcat(readType_seq_str, '%d');
        end
    end

    ldata = textscan(dataFileRead, readType_seq_str, 'CollectOutput', true);
    line = string(ldata{1});

    seq_numOfTracks_in_Lengths = str2double(line);
    %fclose(dataRead_directory_string);
    %max_min_numOfTracks_per_Length = zeros(size(seq_numOfTracks_in_Lengths, 2), 2);
    for i = 1:size(seq_numOfTracks_in_Lengths, 2)
        avg_numOfTracks_per_Length(i,typeIdx) = sum(seq_numOfTracks_in_Lengths(:,i));
        avg_numOfTracks_per_Length(i,typeIdx) = avg_numOfTracks_per_Length(i,typeIdx) / size(seq_numOfTracks_in_Lengths, 1);
        %max_min_numOfTracks_per_Length(i,1) = max(seq_numOfTracks_in_Lengths(:,i));
        %max_min_numOfTracks_per_Length(i,2) = min(seq_numOfTracks_in_Lengths(:,i));
    end
end



% -- drawings --
color_label = ["bo-"; "ro-"; "go-"];
figure;
for typeIdx = 1:size(collection_featureType, 1)
    plot(2:25, avg_numOfTracks_per_Length(:,typeIdx), color_label(typeIdx,1), 'DisplayName', collection_featureType(typeIdx,1));
    hold on;
end
xlim([2, 25]);
xlabel('track length N');
ylabel('average number of tracks');
legend;
set(gcf,'color','w');

if zoomIn
    figure;
    for typeIdx = 1:size(collection_featureType, 1)
        plot(2:25, avg_numOfTracks_per_Length(:,typeIdx), color_label(typeIdx,1), 'DisplayName', collection_featureType(typeIdx,1));
        hold on;
        for x = 19:24
            y = avg_numOfTracks_per_Length(x,typeIdx);
            %text(x+1,y+2,[num2str(round(y,2))])

            hold on;
        end
        hold on;
    end
    xlim([20, 25]);
    ylim([0, 50]);
    set(gcf,'color','w');
end

%% -- draw inlier ratio over different track lengths for SIFT, SURF, and ORB --
sequenceName = 'KITTI-seq05/';
fileName = 'inlier_Ratio_';
collection_featureType = ["SURF"; "ORB"; "SIFT"];
filePostFix = '.txt';

[numOfFrames, window_length] = size(seq_numOfTracks_in_Lengths);
seq_inlier_ratio_in_Lengths = zeros(numOfFrames, 24);
avg_inlier_ratio_in_Lengths = zeros(window_length, 3);

for typeIdx = 1:size(collection_featureType, 1)
    dataRead_directory_string = strcat(fileDirectory, sequenceName, fileName, collection_featureType(typeIdx), filePostFix);
    seq_inlier_ratio_in_Lengths = dlmread(dataRead_directory_string);
    
    % -- compute the average over number of frames --

    for i = 1:window_length
        avg_inlier_ratio_in_Lengths(i,typeIdx) = sum(seq_inlier_ratio_in_Lengths(:,i))/numOfFrames;
        avg_inlier_ratio_in_Lengths(i,typeIdx) = avg_inlier_ratio_in_Lengths(i,typeIdx).*0.87;
    end
end

% -- drawings --
color_label = ["bo-"; "ro-"; "go-"];
figure;
for typeIdx = 1:size(collection_featureType, 1)
    plot(2:25, avg_inlier_ratio_in_Lengths(:,typeIdx), color_label(typeIdx,1), 'DisplayName', collection_featureType(typeIdx,1));
    hold on;
end
xlim([2, 25]);
xlabel('track length N');
ylabel('average inlier ratio');
ylim([0, 1]);
legend;
set(gcf,'color','w');
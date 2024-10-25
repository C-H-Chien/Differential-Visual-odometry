function [matR, T] = read_groundtruths_kitti(p0_FileName)

    poseFileRd = fopen(p0_FileName, 'r');

    ldata = textscan(poseFileRd, '%s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput', true);
    line = string(ldata{1});

    T = zeros(3, size(line, 1));
    matR = zeros(3,3,size(line, 1));
    for i = 1:size(line, 1)
        for t = 1:3
            T(t,i) = line(i,t+3*t);
        end

        for k = 1:3
            for j = 1:3
                matR(k,j,i) = line(i, j+4*(k-1));
            end
        end
    end
    
    fclose(poseFileRd);
end

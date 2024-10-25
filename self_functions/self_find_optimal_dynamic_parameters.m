function w = self_find_optimal_dynamic_parameters(time, arcLength, order)

    % -- 2) build a system of linear equations --
    phi = zeros(order+1, 1);
    matA = zeros(order+1, order+1);
    vecZ = zeros(order+1, 1);
    for indx = 1:size(time, 2)
        phi(1) = 1;
        for i = 1:order
            phi(i+1) = time(1,indx)^i;
        end
        matA = matA + phi * phi';
        vecZ = vecZ + arcLength(indx,1) * phi;
    end

    % -- 3) solve the linear system
    w = matA \ vecZ;

end

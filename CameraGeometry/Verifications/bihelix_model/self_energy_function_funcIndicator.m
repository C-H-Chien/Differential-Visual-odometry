function energy_value = self_energy_function_funcIndicator(x, tangent0, normal0, binormal0, C_GT, junction_pt_idx)

    k0 = x(1);
    tau0 = x(2);
    k1 = x(3);
    tau1 = x(4);
    energy_value = 0;
    
    [C0, ~] = self_generateBihelixFromModel(k0, tau0, k1, tau1, tangent0, normal0, binormal0, C_GT, junction_pt_idx);

    for i = 1:size(C_GT, 2)
        energy_value = energy_value + (C_GT(1,i)-C0(1,i))^2;
    end
    
    %energy_value = energy_value / size(C_GT, 2);
end
function  power_alloc = Heuristic_Power_Allocation(Pmax,bkl,L,K)
    % input parameters:
    % L = number of AP
    % K = number of UE
    % Pmax = maximum downlink power

    % output parameters:
    % power_alloc = K*L matrix, downlink power allocation between UE to AP
    power_alloc = zeros(K,L);
    v = 0.6; % scaling parameter
    power_alloc = sqrt(Pmax) .* (bkl.^v)./sum(bkl.^v,1);
    
end

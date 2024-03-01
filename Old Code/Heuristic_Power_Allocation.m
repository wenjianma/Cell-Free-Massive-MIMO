function  rho_pp = Heuristic_Power_Allocation(Pmax,bkl,L,K)
    % input parameters:
    % L = number of AP
    % K = number of UE
    % Pmax = maximum downlink power

    % output parameters:
    % rho_pp = K*L matrix, downlink power allocation between UE to AP
    %rho_pp = zeros(K,L);
    v = 0.8; % scaling parameter
    rho_pp = Pmax .* ((bkl.^v)./sum(bkl.^v,1));
%     rho_pp = zeros(K, L);
%     v = 0.6;

%     for k = 1:K
%         for l = 1:L
%             denominator = sum(bkl(k, :).^v);
%             rho_pp(k, l) = sqrt(Pmax) * (bkl(k, l)^v) / denominator;
%         end
%     end
  
end

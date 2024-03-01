%% General Description
% Function to allocate how much power each AP is using to serve its UEs
% The allocation policy is based on Large Scale Fading (bkl)

%% Function Description

% Input parameters:
% L = number of AP
% K = number of UE
% Pmax = maximum downlink power

% Output parameters:
% rho_pp = matrix of downlink power allocation between UE to AP K*L
%% Function Part
function  rho_pp = Heuristic_Power_Allocation(Pmax,bkl,L,K, sortedMaxIndices)
    v = 0.6; % scaling parameter
    rho_pp = zeros(K,L);

    for l = 1:L
        [~,served_UEs] = find(sortedMaxIndices==l);
        if isempty(served_UEs) == 0
            rho_pp(served_UEs,l) = Pmax.* ( bkl(served_UEs,l).^v ) ./ sum((bkl(served_UEs,l).^v),1);
        end
        % for kk = 1:size(served_UEs)
        %     if served_UEs(kk) == 0
        %         break;
        %     end
        % 
        % 
        % end

    end


    % for k = 1:K
    %     for l = 1:3 
    %         rho_pp(k,sortedMaxIndices(l,k)) = Pmax .* ((bkl(k,sortedMaxIndices(l,k)).^v)./sum(bkl(k,sortedMaxIndices(l,k)).^v)); 
    %     end 
    % end 

    % v = 0.6; % scaling parameter
    % rho_pp = Pmax .* ((bkl.^v)./sum(bkl.^v,1)); 
end
%% This function is used for spectral efficiency calculation for all UEs
% In this function, the spectral efficiency of each user is calculated
% first and then plot the CDF
% This function is based on Corollary 6.3 of Reference 1

%% Input of this function
% prelogFactor = tau_d / tau_c
% w = normalized precoding vector
% power_alloc = K*L size matrix, downlink transmission power of UE k to AP l under heuristic power allocation policy
% R = N*N*K*L size matrix, R(:,:,k,l) is the spatial correlation between UE k and AP l
% D = K*L size matrix, D(k,l) = 1 means UE k is served by AP l
% H = N*Number_of_channelrealizations*K*L, true channel condition
% Hhat = N*Number_of_channelrealizations*K*L, estimated channel
% numOfRealization = number of the channel realization
% K = number of UE
% L = number of AP
%% Important variables
% numerator = size K vector 
% denominator = size K vector

%notes: this version hasn't add noise 
%%
function  SE_MR = Spectral_Efficiency_Calculation(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L)
    
    numerator = zeros(K,1);
    denominator = zeros(K,1);
    SE_MR = zeros(K,1);
    % Go through all the UE
    for k =1:K
        
        %calculate numerator
        
        % For each UE, go through all the AP (all ue is served by all ap)
        val0 = zeros(L,1);
        for l = 1:L
             %Go through all channel realization
             val1 = 0;
             for n = 1:numOfRealization
                val1 = val1 + (conj(H(:,n,k,l)))'*(D(:,:,k,l)*w_MR(:,k,l)); % w is normalized so we need to add power allocation in the next line
             end
             %val1 = sum(power_alloc(k,:))*val1/numOfRealization; % calculate the expectation
             val1 = val1/numOfRealization;
             val0(l,1) = val1;
        end
        numerator(k,1) = abs(sum(val0))^2;
        
        %calculate the first part of denominator
        val2 = zeros(K,1);
        for k2 = 1:K
            for n = 1:numOfRealization
                val3 = 0;
                for l = 1:L
                   val3 = val3 + (conj(H(:,n,k,l)))'*(D(:,:,k2,l)*w_MR(:,k2,l));
                end
                %val3 = abs(sum(power_alloc(k2,:))*val3)^2;
                val3 = abs(val3)^2;
                val2(k2,1) = val2(k2,1) + val3;
            end
            val2(k2,1) = val2(k2,1)/numOfRealization;
        end
        denominator(k,1) = sum(val2);

        %SE calculation
        SE_MR(k,1) = prelogFactor*real(log2(1+numerator(k,1)/(denominator(k,1)-numerator(k,1)+1)));
    end

    

end
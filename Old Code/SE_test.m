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
function  SE_MR = SE_test(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L)
    
    SINR = zeros(K,1);
    numerator = zeros(K,1);
    denominator = zeros(K,1);
    SE_MR = zeros(K,1);
    for k = 1:K
        all_mean = 0;
        for l = 1:L
            hdw = 0;
            for no = 1:numOfRealization
                t1 = power_alloc(k,l);
                t2 = H(:, no, k, l)';
                t3 = D(:, :, k, l);
                t4 = w_MR(:, no,k, l);
                hdw = hdw + sqrt(power_alloc(k,l))* H(:, no, k, l)' * D(:, :, k, l) * w_MR(:, no,k, l)/numOfRealization;
            end
            all_mean = all_mean + hdw;
        end
        % numerator k
        numerator(k,1) = abs(all_mean)^2;

        interference = 0;     
        for i = 1:K
            hdw3 = 0;
            for no = 1: numOfRealization
                hdw2 = 0;
                for l = 1:L
                    hdw2 = hdw2 + sqrt(power_alloc(i,l))*H(:, no, k, l)' * D(:, :, i, l) * w_MR(:,no, i, l); 
                end
                hdw3 = hdw3 + abs(hdw2).^2/numOfRealization;
            end
            interference = interference + hdw3;
        end
        denominator(k,1) = interference - numerator(k,1) + 1; % think about noise

        SINR(k,1) = numerator(k,1)/denominator(k,1);
        SE_MR(k,1) = prelogFactor*real(log2(1+SINR(k,1)));
    end 

end

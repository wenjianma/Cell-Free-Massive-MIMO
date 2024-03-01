%% General Description
% Function to calculate the precoding vector using Maximum Ratio (MR) and Regulated Zero Forcing (RZF)

%% Function Description

% Input parameters:
% Hhat               =  Estimated channel matrix
% p                  =  Uplink transmit power per UE
% numRealizations    =  Number of channel realizations
% N                  =  Number of antennas per Access Point
% K                  =  Number of User Equipments per cell
% L                  =  Number of Access Points
% sortedMaxIndices   =  CxK matrix where C is the number of APs that are serving one UE
% C                  =  Cluster size.

% Output parameters:
% w_MR               =  Normalized Precoded MR vector  N*numofRealization*K*L
% w_RZF              =  Normalized Precoded RZF vector N*numofRealization*K*L

%% Function Part
function [w_MR,w_RZF] = MR_RZF_Precoder(Hhat, p, numRealizations, N, K, L, sortedMaxIndices, C)

    I_N = eye(N);

    % Initialize the MR and RZF vectors
    v_MR = zeros(N,numRealizations,K,L);
    v_RZF = zeros(N,numRealizations,K,L);

    for k = 1:K
        for l_idx = 1:C % Loop over the 3 APs serving UE k
            if ~(sortedMaxIndices(l_idx, k) == 0)
                l = sortedMaxIndices(l_idx, k); % AP index from sortedMaxIndices

                % MR Vector Calculation
                v_MR(:, :, k, l) = Hhat(:, :, k, l);

                % RZF Vector Calculation
                for n = 1:numRealizations
                    sum_term = zeros(N, N);
                    for i = 1:K
                        if ismember(l, sortedMaxIndices(:, i)) % Check if AP l serves UE i
                            sum_term = sum_term + p * Hhat(:, n, i, l) * Hhat(:, n, i, l)';
                        end
                    end
                    inv_term = inv(sum_term + I_N);
                    v_RZF(:, n, k, l) = inv_term * p * Hhat(:, n, k, l);
                end
            end
        end
    end

    % Normalization of Precoded Vectors
    w_MR = v_MR ./ vecnorm(v_MR);
    w_RZF = v_RZF ./ vecnorm(v_RZF);

    % Fix NaN values by converting them to zeros
    a = find(isnan(w_MR));
    b = find(isnan(w_RZF));
    w_MR(a) = 0;
    w_RZF(b) = 0;

end

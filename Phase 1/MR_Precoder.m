% Function to Implement MR Precoder 

% Input parameters:
% p    = the transmit power allocated to UE i
% Hhat = Estimated channel matrix.

% Output parameters:
% w = precoder vector that AP l assigns to UE i. 

function  [w] = MR_Precoder(Hhat,p,L,K)

% MR Precoding for each UE and AP combination
for l = 1:L
    for i = 1:K
        % Calculate the MR precoding weight for the i-th UE
        w = sqrt(p) * Hhat(:, i) / sqrt(mean(abs(Hhat(:, i)).^2));
    end
end
end

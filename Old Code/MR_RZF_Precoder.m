% Function to calculate the precoder using Maximum Ratio (MR) and 
% Regulated Zero Forcing (RZF)

% Input parameters:
% Hhat = Estimated channel matrix.
%power              =  Uplink transmit power per UE
%numRealizations    =  Number of channel realizations
%numAntennas        =  Number of antennas per Access Point
%numUEs             =  Number of User Equipments per cell
%numAPs             =  Number of Access Points

% Output parameters:
% w_MR = Precoded MR vector. 
% w_RZF = Precoded RZF vector. 

function [w_MR,w_RZF] = MR_RZF_Precoder(Hhat,p,numRealizations,N,K,L)

    % Store identity matrix of size numAntennas x numAntennas
    I_N = eye(N);


        %Calculation of the Vector
        v_MR = Hhat;

        v_RZF = zeros(N,numRealizations,K,L); 
        for l = 1:L 
            for k = 1:K
                for n = 1:numRealizations
                    sum_term = zeros(N,N);
                    for i = 1:K
                        sum_term = sum_term + p*Hhat(:,n,i,l)*Hhat(:,n,i,l)';
                    end 
                    inv_term = inv(sum_term + I_N);
                    v_RZF(:,n,k,l) = inv_term*p * Hhat(:,n,k,l);
                end 
            end 
        end 

       
% Calculating Normalizer Precoded Vector 
        w_MR = v_MR ./ pagenorm(v_MR);
        w_RZF = v_RZF ./ pagenorm(v_RZF);

        %end
    %end
end

% % Function to calculate the precoder using Maximum Ratio (MR) and 
% % Regulated Zero Forcing (RZF)
% 
% % Input parameters:
% % Hhat = Estimated channel matrix.
% %power              =  Uplink transmit power per UE
% %numRealizations    =  Number of channel realizations
% %numAntennas        =  Number of antennas per Access Point
% %numUEs             =  Number of User Equipments per cell
% %numAPs             =  Number of Access Points
% 
% % Output parameters:
% % w_MR = Precoded MR vector. 
% % w_RZF = Precoded RZF vector. 
% 
function [w_MR,w_RZF] = MR_ZRF_Precoder(Hhat,p,numRealizations,N,K,L)

    % Prepare to store results through initilization
    %w_MR = zeros(numAntennas,numRealizations,numUEs,numAPs);
    %w_RZF = zeros(numAntennas,numRealizations,numUEs,numAPs);

    % Store identity matrix of size numAntennas x numAntennas
    I_N = eye(N);

    % Normalizing the precoder
    %for nr = 1:numRealizations

        %for ap = 1:numAPs

        %Calculation of the Vector
        v_MR = Hhat;
        %v_RZF = (((v_MR*pagectranspose(v_MR))+identityMatrix)*(-1))*v_MR;
        %v3 = pagectranspose(v_MR);
        %v2 = pagemtimes(v_MR, v3);
        %v1 = (power*v2);
        %v4 = (v1+identityMatrix);
        %v5 = v4.^(-1);
        %v_RZF = pagemtimes(v5, v_MR);
%         v_RZF = pagemtimes(((power*(pagemtimes(v_MR,pagectranspose(v_MR))))+identityMatrix).^(-1), v_MR);

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

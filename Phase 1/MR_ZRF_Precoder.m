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

function [w_MR,w_RZF] = MR_ZRF_Precoder(Hhat,power,numRealizations,numAntennas,numUEs,numAPs)

    % Prepare to store results through initilization
    % w_MR = zeros(numAntennas,numUEs,numAPs);
    % w_RZF = zeros(numAntennas,numUEs,numAPs);

    % Store identity matrix of size numAntennas x numAntennas
    identityMatrix = eye(numAntennas);
    v_MR = Hhat; 
    w_MR = v_MR./pagenorm(v_MR);

    % Normalizing the precoder
%     for nr = 1:numRealizations
%         for ap = 1:numAPs
%             v_MR = Hhat(:,nr,:,ap);
%             w_MR(:,:,ap) = v_MR./pagenorm(v_MR);
%         end 
%     end 

    w_RZF = 0;


end
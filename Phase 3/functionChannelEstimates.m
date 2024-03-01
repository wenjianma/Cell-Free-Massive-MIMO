%% General Description
% This function generate channel estimation for cell free massive MIMO using MMSE estimator and also we assumed pilot are orthogonal to each other
%% Function Description

% Input parameters:
% cMat               =  Correlation matrix between each Access Point (AP) and User Equipment (UE). 
% numRealizations    =  Number of channel realizations
% numAPs             =  Number of Access Points
% numUEs             =  Number of User Equipments per cell
% numAntennas        =  Number of antennas per Access Point
% numOPilots         =  Number of orthogonal pilots
% pilotIndex         =  Vector containing the pilot assigned to each UE
% power              =  Uplink transmit power per UE

% Output parameters:
% H = the actual channel matrix numAntennasxnumRealizationxnumUEsxnumAPs (NxnumOfRealizationxKxL)
% Hhat = the estimated channel matrix numAntennasxnumRealizationxnumUEsxnumAPs (NxnumOfRealizationxKxL)

%% Function Part
function [H, Hhat] = functionChannelEstimates(cMat, numRealizations, numAPs, numUEs, numAntennas, numOPilots, pilotIndex, power)

    % Generate uncorrelated Rayleigh fading channel realizations
    CH = sqrt(0.5)*(randn(numAntennas,numRealizations,numUEs,numAPs)+1i*randn(numAntennas,numRealizations,numUEs,numAPs));
    H = zeros(numAntennas,numRealizations,numUEs,numAPs);

    % Apply the spatial corelation matrices to all channels
    for ap = 1:numAPs
        for ue = 1:numUEs
            % Apply correlation to the uncorrelated channel realizations
            cSqrt = sqrtm(cMat(:,:,ue,ap));
            H(:,:,ue,ap) = cSqrt*CH(:,:,ue,ap);
        end
    end

    % Perform channel estimation
    % Store identity matrix of size numAntennas x numAntennas
    identityMatrix = eye(numAntennas);

    % Generate realizations of normalized noise
    NormalizedNoise = sqrt(0.5)*(randn(numAntennas,numRealizations,numAPs,numOPilots) + 1i*randn(numAntennas,numRealizations,numAPs,numOPilots));

    % Prepare to store results
    Hhat = zeros(numAntennas,numRealizations,numUEs,numAPs);

    % Go through all APs and pilots
    for ap = 1:numAPs
    
        for tk = 1:numOPilots
    
            % Compute processed pilot signal for all UEs that use pilot
            yp = sqrt(power*numOPilots)*sum(H(:,:,pilotIndex == tk,ap),3) + NormalizedNoise(:,:,ap,tk);
            
            % Compute the matrix that is inverted in the MMSE estimator 
            iMat = (power*numOPilots*sum(cMat(:,:,pilotIndex == tk,ap),3) + identityMatrix); % we just add identityMatrix here because we have normalized the signal and the noise 
            
            % Go through all UEs that use pilot and compute the MMSE Hhat
            UE = find(pilotIndex == tk);
    
            for ue = 1:size(UE)
                ciMat = cMat(:,:,UE(ue),ap)*iMat^(-1);
                Hhat(:,:,UE(ue),ap) = sqrt(power*numOPilots)*ciMat*yp;
            end
            
        end
    
    end

end
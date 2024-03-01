%% General Description
% This function calculates the spatial correlation matrix under correlated rayleigh fading

%% Function Description
% Input parameters:
% N                 =   num of antennas
% theta_deg         =   angle between AP and UE
% ASD_deg           =   angular standard deviation
% antennaSpacing    =   distance between antennas

% Output parameters:
% R                 =   spatial correlation matrix N*N

%% Function Part
function R = calculateR(N, theta_deg, ASD_deg,antennaSpacing)
    
    ASD = deg2rad(ASD_deg);
    theta = deg2rad(theta_deg);
    
    % Initialize the first row of the correlation matrix
    firstRow = zeros(N, 1);

    for l = 1:N
        % Compute the approximated integral as in (2.24)
        % we use l-1 for l-m because we are calculating the distance from the first antenna ;
        firstRow(l) = exp(1i*2*pi*antennaSpacing*sin(deg2rad(theta))*(l - 1)) * exp(-ASD^2 / 2 * (2 * pi * antennaSpacing * cos(deg2rad(theta_deg)) * (l - 1))^2); 
    end

    % Construct the Toeplitz matrix using the first row
    R = toeplitz(firstRow);
end
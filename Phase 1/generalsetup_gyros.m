%% General Description
% Function to generate the general setup of our network simulator. 
% Here we generate the AP locations based on a grid format with randomness
% factor. The UE locations are then generated to be randomly and uniformly distributed over the area. 
%
% In this function we also define our channel model (noise, pathloss,
% shadow fading) based on --------. 
%% Function Description 
% INPUTS: 
% K = number of UEs
% L = number of APs
% N = number of antennas per AP
% pilotLength = length of orthogonal pilots 
% numOfSim = number of Monte-Carlo Simulations
% 
% OUTPUTS: 
% R = NxNxKxLxnumOfsim array that corresponds to the spatial correlation
% matrix between UEs and APs. 
%
% bklT = KxLxnumOfSim array that corresponds to the channel gain or bkl
% coefficients between UEs and APs 
%
% pilotIndex = KxnumOfSim array that shows the assigned pilot for each UE. 
%% Begin Code
function [R,bklT,pilotIndex] = generalsetup_gyros(K, L, N,pilotLength,flag)

% Define the side length of our area
areaSizeX = 2000; 

% Distance between two antennas = half-wavelength
antenna_spacing = 1/2; 

% AP antenna height based on 3GPP Urban Microcell Model 
h_bs = 10; % 10m 

% UE antenna height based on 3GPP Urban microcell model 
h_ue = 1; 
%% Propagation Model

%Total channel bandwidth in Hz
W = 20e6; 

%Noise figure in dB
noise_fig = 7;

%Calculate total receiver noise power
%thermal noise = boltzmann constant*room temperature = -174 dBm
thermalNoise = -174;
%Total receiver noise power. 
noiseVardBm = thermalNoise + 10*log10(W) + noise_fig; %in dBm
%noiseVardbW = noiseVardBm - 30; % in dBW
% Pathloss exponent
a = 3.76; 

% Standard deviation of shadow fading, based on 3GPP Urban Microcell
% pathloss model.
sigma_sf = 4; 

%Average channel gain in dB at a reference distance of 1 meter. Note that
% -35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;

% Angular Standard Deviation for the angle needed in R (in degrees)
ASDdeg = 10;

%% Generate AP and UE locations
% During simulation, the locations of APs should be changed. The
% locations of UEs can be changed but we need to make sure two
% algorithms have the same UE distribution. We can use "the same seed"
% to control it.

% The distribution of AP can be grid distribution (in this case, sqrt(number of AP) should be integer)

% number of APs per row/column
% note on this: example for L = 3, then M = 4 but we don't have 4 APs. 
M = round(sqrt(L)); 

% distance between two neighbor APs
D = areaSizeX/(M+1); 



stream1 = RandStream('mt19937ar', 'Seed', 123);
[x,y] = meshgrid(D:D:areaSizeX - D,D:D:areaSizeX - D);
% add randomness to the grid distribution, maximized randomness equals to D/2
APpositions = (x + rand(stream1,M,M)*D/2) + 1i*(y + rand(stream1,M,M)*D/2);
%convert M*M matrix AP to (M^2)*1 vector
APpositions = APpositions(:); 
%APpositions_wrapped = APpositions;

% Random UE locations with uniform distribution
UEpositions = (rand(stream1,K,1) + 1i*rand(stream1,K,1))*areaSizeX;

%% Gryos Method
% Plot AP and UE locations

% figure(1)
% scatter(real(APpositions), imag(APpositions), 'DisplayName', 'APs')
% hold on;
% scatter(real(UEpositions), imag(UEpositions), "filled", 'DisplayName','UEs')
% legend
% xlabel('Area - X (m)')
% ylabel('Area - Y(m)')

% Use wrap around method to compute the distances between APs-UEs.
% If calculated distance is larger than half of the side length then
% correct_distance = sidelength - distance. Otherwise, correct distance
% = distance.

distances_old = abs(APpositions - UEpositions.');
% distances is a variable only for the current simulation. 
distances = abs(APpositions - UEpositions.');
% Take into consideration APs antenna height.
distances = sqrt((h_bs-h_ue)^2 + distances.^2);

 % we will also store the angle between each UE-AP and then modify it to
% the correct angle.
%angles = atan2(imag(APpositions)-imag(UEpositions.'),real(APpositions)-real(UEpositions.'));

%Calculate the virtual AP positions if they need to be wrapped
%So here the logic is as follows:
% 1. find the indices where we need to apply wrap around 
% 2. convert to row-column format so that we can have the AP-UE mindset
% 3. Now for each case where we need to apply the wrap around distance
% if the real part of the AP position is larger than half of the side
% length then the new real part is sidelength - the previous real part.
% And then if the imaginary part is larger than half of the area size
% then the same logic applies. Otherwise we keep the same values. 
% Similar logic for when the real part is smaller than half of the side
% length. 

% inside this loop we change the angles that need to be calculated
% based on the wrap around positions.

% wrapAroundIndices = find(distances > areaSizeX/2);
% for idx = 1:length(wrapAroundIndices)
%     [x, y] = ind2sub(size(APpositions), wrapAroundIndices(idx));
% 
%     if real(APpositions(x,1)) >= areaSizeX/2 
%         APpositions_wrapped(x,1) = areaSizeX - real(APpositions(x,1)); 
%         if imag(APpositions(x,1)) >= areaSizeX/2
%             APpositions_wrapped(x,1) = areaSizeX*1i - imag(APpositions(x,1)); 
%         else 
%             APpositions_wrapped(x,1) = imag(APpositions(x,1));
%         end 
%     else 
%         APpositions_wrapped(x,1)= real(APpositions(x,1));
%         if imag(APpositions(x,1)) >= areaSizeX/2 
%             APpositions_wrapped(x,1) = areaSizeX*1i - imag(APpositions(x,1));
%         else 
%             APpositions_wrapped(x,1) = imag(APpositions(x,1));
%         end 
%     end 
% 
%     % angles(x,y) = atan2(imag(APpositions_wrapped(x,1)) - imag(UEpositions(y,1)),real(APpositions_wrapped(x,1))-real(UEpositions(x,1)));
% end
% 
% 
distances(distances > areaSizeX/2) = distances(distances > areaSizeX/2) - areaSizeX/2;
% anglesDeg = rad2deg(angles);
     

%% Compute channel gain and channel gain over noise
% 
channelgaindB = constantTerm - a*10*log10(distances) + sigma_sf*randn([L,1]);
channelgainOverNoisedB = channelgaindB - noiseVardBm; 

%% Calculate βkl 
bkl = 10.^(channelgainOverNoisedB/10);

%% Calculate R matrix for uncorrelated Rayleigh fading 
%Identity matrix of size N, i.e., number of antennas per AP.
I = eye(N);

% need to transpose bkl to match dimensions
% use page-..- to do operations on N-D arrays.
bklT = bkl';
% anglesDegT = anglesDeg';
%find the correlation matrix 
%find the correlation matrix 
Runcorr = zeros(N, N, K, L);
Rcorr = zeros(N,N,K,L);
Rcorr_final = zeros(N,N,K,L);
for k = 1:K
    for l = 1:L
        %using uncorrelated Rayleigh fading 
        Runcorr(:,:,k,l) = bklT(k,l).*I;
        %using correlated Rayleigh fading
        % Rcorr(:,:,k,l) = calculateR(N,anglesDegT(k,l),ASDdeg,antenna_spacing);
        % Rcorr_final(:,:,k,l) = bklT(k,l)*Rcorr(:,:,k,l);
    end 
end 

if strcmp("uncorrelatedRayleigh",flag) == 1
    R = Runcorr;
elseif strcmp("correlatedRayleigh",flag) == 1
    R = Rcorr_final; 
else 
    fprintf("wrong flag indication");
end 

%% Pilot Assignment 
% Do random pilot assignment to each UE 
pilotIndex = randi([1 pilotLength],[K 1]);

    
end 
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
function [R,bklT,pilotIndex,sortedMinIndices] = generalsetup(K, L, N,pilotLength,flag)

% Define the side length of our area
areaSizeX = 1000; 

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
constantTerm = -30.5;

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


stream1 = RandStream('mt19937ar', 'Seed', randi(10000));
[x,y] = meshgrid(D:D:areaSizeX - D,D:D:areaSizeX - D);
%add randomness to the grid distribution, maximized randomness equals to D/2
APpositions = (x + rand(stream1,M,M)*D/2) + 1i*(y + rand(stream1,M,M)*D/2);
%convert M*M matrix AP to (M^2)*1 vector
APpositions = APpositions(:); 
%APpositions_wrapped = APpositions;
% 
% % Random UE locations with uniform distribution
% UEpositions = (rand(stream1,K,1) + 1i*rand(stream1,K,1))*areaSizeX;
%Random AP locations with uniform distribution
% APpositions = (rand(stream1, L,1) + 1i*rand(stream1, L,1)) * areaSizeX;
% 
%Random UE locations with uniform distribution
UEpositions = (rand(K,1) + 1i*rand(stream1, K,1)) * areaSizeX;

%Compute alternative AP locations by using wrap around 
wrapHorizontal = repmat([-areaSizeX 0 areaSizeX],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

% plot the locations of the APs & UEs
 scatter(real(APpositionsWrapped), imag(APpositionsWrapped))
 figure(2)
 scatter(real(APpositionsWrapped), imag(APpositionsWrapped),'DisplayName','APs')
 hold on;
 scatter(real(UEpositions),imag(UEpositions),"filled",'DisplayName','UEs')
 legend
 xlabel('Area - X (m)')
 ylabel('Area - Y (m)')

distances = zeros(L,K);
Runcorr = zeros(N, N, K, L);
Rcorr = zeros(N,N,K,L);
Rcorr_final = zeros(N,N,K,L);

for current_UE = 1:K
        %replicate the position of the UE under study so that the
        %APpositionsWrapped matrix and the UEPositionWrapped matrix are of
        %the same size. 
        %UEPositionWrapped = repmat(UEpositions(current_UE),size(APpositionsWrapped)); 

        %calculate the minimum absolute distance between the current UE and
        %all of the APs (in wrap around area).
        %min(..., [], 2) calculates the minimum value along the second dimension (columns). 
        [distancetoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(current_UE),size(APpositionsWrapped))),[],2);
        %disp('distancetoUE:');
        %disp(distancetoUE);
        distances(:,current_UE) = sqrt(9^2+distancetoUE.^2);
        %disp('distances:');
        %disp(distances(:, current_UE));
        %disp(distances);

        %save the number of the "cell" (=expanded area, for wrap around method) where the minimum distance from UE to AP is found.
       % whichpos_tot(:,current_UE) = whichpos; 

        %now we need to determine the channel gain based on the large scale
        %parameters in dB and in linear scale
        %channel gain is parameter βkl in book
        
        sf = sigma_sf*randn(size(distances(:,current_UE)));
        channel_gainOverNoisedB(:,current_UE) = constantTerm - a*10*log10(distances(:,current_UE)) + sf - noiseVardBm;
        channel_gainOverNoise(:,current_UE) = db2pow(channel_gainOverNoisedB(:,current_UE));
        %make the channel gain in W not mW 
        
        %calculate small-scale parameters, i.e., normalized R matrix (spatial
        %correlation matrix) for all of the APs! 
        %normalized R matrix = channel_gain*R 
        % we are calculating R between AP l (current AP) and UE k (current
        % UE

        % for uncorrelated rayleigh fading R = identity matrix
        for current_AP = 1:L
            Runcorr(:,:,current_UE,current_AP) = channel_gainOverNoise(current_AP,current_UE)*eye(N);
            thetaAPtoUE = angle(UEpositions(current_UE)-APpositionsWrapped(current_AP,whichpos(current_AP)));
            Rcorr(:,:,current_UE,current_AP) = calculateR(N,thetaAPtoUE,ASDdeg,antenna_spacing);
            Rcorr_final(:,:,current_UE,current_AP) = channel_gainOverNoise(current_AP,current_UE)'*Rcorr(:,:,current_UE, current_AP);
        end 

end
disp(distances);

% Added code changes for Clustering pairing up APs that are close to a
% specific UE
% Initialize arrays to store two minimum distances and their indices for each UE

minDistances = zeros(2, K);
minIndices = zeros(2, K);

for i = 1:K
    % Find the two minimum distances and their indices for the current UE (column-wise)
    [sortedDistances, sortedIndices] = sort(distances(:, i));
    minDistances(:, i) = sortedDistances(1:2);
    minIndices(:, i) = sortedIndices(1:2);
    sortedMinIndices = sort(minIndices);
end
disp(minDistances);
disp(sortedMinIndices);


%% Compute channel gain and channel gain over noise
bklT = channel_gainOverNoise;
bklT = bklT';

if strcmp("uncorrelatedRayleigh",flag) == 1
    R = Runcorr;
elseif strcmp("correlatedRayleigh",flag) == 1
    R = Rcorr_final; 
else 
    fprintf("wrong flag indication");
end 

%% Pilot Assignment 
% Do random pilot assignment to each UE
%pilotIndex = ones(K,1);
%pilotIndex = randi([1 pilotLength],[K 1]);
pilotIndex = simple_algorithm(R,pilotLength,K,L,sortedMinIndices);
disp(pilotIndex);
    
end 

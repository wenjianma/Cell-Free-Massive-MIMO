%% General Description
% Function to generate the general setup of our network simulator. 
% Here we generate the AP locations based on a grid format with randomness factor. 
% The UE locations are then generated to be randomly and uniformly distributed over the area. 
% In this function we also define our channel model (noise, pathloss, shadow fading) based on 3GPP Urban Microcell Model. 
%% Function Description

% INPUTS: 
% K = number of UEs  scalar
% L = number of APs  scalar
% N = number of antennas per AP  scalar
% pilotLength = length of orthogonal pilots  scalar
% flag_ch = flag of channel type  scalar
% flag_alg = flag of algorithm type  scalar

% OUTPUTS: 
% R = array that corresponds to the spatial correlation matrix between UEs and APs  size: NxNxKxL 
% bklT = array that corresponds to the channel gain coefficients between UEs and APs  size: KxL
% pilotIndex = vector that shows the assigned pilot for each UE  size: K*1

%% Function part
function [R,bklT,pilotIndex,sortedMaxIndices] = generalsetup(K, L, N,pilotLength,flag_ch, flag_alg, flag_serv, C)

%% Parameter Initilization & Propagation Model

% Define the side length of our area (unit: meter)
areaSizeX = 1000; 

% Distance between two antennas (1/2 means half-wavelength)
antenna_spacing = 1/2; 

% AP antenna height based on 3GPP Urban Microcell Model (unit: meter) 
h_ap = 10;  

% UE antenna height based on 3GPP Urban microcell model (unit: meter) 
h_ue = 1; 

% Total channel bandwidth in Hz
W = 20e6; 

% Noise figure in dB
noise_fig = 7;

% Calculate total receiver noise power
% thermal noise = boltzmann constant*room temperature = -174 dBm
thermalNoise = -174;

% Total receiver noise power in dBm
noiseVardBm = thermalNoise + 10*log10(W) + noise_fig; 

% Pathloss exponent
a = 3.76; 

% Standard deviation of shadow fading, based on 3GPP Urban Microcell pathloss model.
sigma_sf = 4; 

%Average channel gain in dB at a reference distance of 1 meter. Note that -30.5 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -30.5;

% Angular Standard Deviation for the angle needed in R (in degrees)
ASDdeg = 10;

%% Generate AP and UE locations & Calculate LSF & Spatial Correlation Matrix

% The distribution of AP can be grid distribution (in this case, sqrt(number of AP) should be integer)

% number of APs per row/column (note on this: example for L = 3, then M = 4 but we don't have 4 APs.) 
M = round(sqrt(L)); 

% distance between two neighbor APs
D = areaSizeX/(M+1); 

stream1 = RandStream('mt19937ar', 'Seed', randi(10000));
[x,y] = meshgrid(D:D:areaSizeX - D,D:D:areaSizeX - D);

% add randomness to the grid distribution, maximized randomness equals to D/2
APpositions = (x + rand(stream1,M,M)*D/2) + 1i*(y + rand(stream1,M,M)*D/2);

%convert M*M matrix AP to (M^2)*1 vector
APpositions = APpositions(:); 

% Random UE locations with uniform distribution
UEpositions = (rand(stream1,K,1) + 1i*rand(stream1, K,1)) * areaSizeX;

% Compute alternative AP locations by using wrap around 
wrapHorizontal = repmat([-areaSizeX 0 areaSizeX],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

% plot the locations of the APs & UEs (for test)
%  scatter(real(APpositionsWrapped), imag(APpositionsWrapped))
%  figure(2)
%  scatter(real(APpositionsWrapped), imag(APpositionsWrapped),'DisplayName','APs')
%  hold on;
%  scatter(real(UEpositions),imag(UEpositions),"filled",'DisplayName','UEs')
%  legend
%  xlabel('Area - X (m)')
%  ylabel('Area - Y (m)')

distances = zeros(L,K);
Runcorr = zeros(N, N, K, L);
Rcorr = zeros(N,N,K,L);
Rcorr_final = zeros(N,N,K,L);

for current_UE = 1:K
    % calculate the minimum absolute distance between the current UE and all of the APs (in wrap around area).
    % min(..., [], 2) calculates the minimum value along the second dimension (columns). 
    [distancetoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEpositions(current_UE),size(APpositionsWrapped))),[],2);
    distances(:,current_UE) = sqrt((h_ap - h_ue)^2+distancetoUE.^2);

    %now we need to determine the channel gain based on the large scale parameters in dB and in linear scale (channel gain is parameter βkl in book)
    sf = sigma_sf*randn(size(distances(:,current_UE)));
    channel_gainOverNoisedB(:,current_UE) = constantTerm - a*10*log10(distances(:,current_UE)) + sf - noiseVardBm;
    channel_gainOverNoise(:,current_UE) = db2pow(channel_gainOverNoisedB(:,current_UE));
    
    %calculate small-scale parameters, i.e., normalized R matrix (spatial correlation matrix) for all of the APs normalized R matrix = channel_gain*R 
    % we are calculating R between AP l (current AP) and UE k (current UE) 
    % for uncorrelated rayleigh fading R = identity matrix
    % for correlated Rayleigh, calculate R using calculateR function
    for current_AP = 1:L
        Runcorr(:,:,current_UE,current_AP) = channel_gainOverNoise(current_AP,current_UE)*eye(N);
        thetaAPtoUE = angle(UEpositions(current_UE)-APpositionsWrapped(current_AP,whichpos(current_AP)));
        Rcorr(:,:,current_UE,current_AP) = calculateR(N,thetaAPtoUE,ASDdeg,antenna_spacing);
        Rcorr_final(:,:,current_UE,current_AP) = channel_gainOverNoise(current_AP,current_UE)'*Rcorr(:,:,current_UE, current_AP);
    end 

end

bklT = channel_gainOverNoise;
bklT = bklT';
%% Cluster Formation

% Added code changes for Clustering pairing up APs that are close to a specific UE
% Initialize arrays to store two minimum distances and their indices for each UE

maxIndices = zeros(C, K);

% for i = 1:K
%     % Find the two minimum distances and their indices for the current UE (column-wise)
%     [sortedDistances, sortedIndices] = sort(distances(:, i));
%     %minDistances(:, i) = sortedDistances(1:2);
%     maxIndices(:, i) = sortedIndices(1:C);
%     sortedMaxIndices = sort(maxIndices);
% end

% clustering formation algorithm 
if strcmp(flag_serv, "cluster") == 1 
    for i = 1:K
        % Find the two minimum distances and their indices for the current UE (column-wise)
        [sortedGains, sortedIndices] = sort(bklT(i, :));
        % minDistances(:, i) = sortedGains(1:2);
        maxIndices(:, i) = sortedIndices(L-C+1:L);
        sumchannelgain = sum(sortedGains(L-C+1:L));
        for cc = 1:C
           prec = sortedGains(L-C+cc)/sumchannelgain;
           if  prec < 0.05 % 0.1 is our threshold
                maxIndices(cc,i) = 0;
           end
        end
        sortedMaxIndices = sort(maxIndices);

    end 
elseif strcmp(flag_serv, "all") == 1 % no clustering
    sortedMaxIndices = repmat((1:L)', 1, K);
elseif strcmp(flag_serv, "predetermined") == 1 % pre-determined clustering
    sortedMaxIndices = predetermined_clustering(K,L,bklT);
else 
    fprintf("wrong flag indication");
end 


%Form the clusters using K-Mean algorithm  
% UE coordinates
 UE_X = real(UEpositions);
 UE_Y = imag(UEpositions);


% disp(minDistances);
% disp(sortedMinIndices);

%% Compute channel gain and channel gain over noise

if strcmp("uncorrelatedRayleigh",flag_ch) == 1
    R = Runcorr;
elseif strcmp("correlatedRayleigh",flag_ch) == 1
    R = Rcorr_final; 
else 
    fprintf("wrong flag indication");
end 

%% Pilot Assignment 

if strcmp("random",flag_alg) == 1
    % Do random pilot assignment to each UE
    % pilotIndex = randi([1 pilotLength],[K 1]);
    for i = 1:K
        if  mod(i,pilotLength) == 0
            pilotIndex(i,1) = pilotLength;
        else
            pilotIndex(i,1) = mod(i,pilotLength);
        end
    end
elseif strcmp("simple",flag_alg) == 1
    % Do pilot assignment to each UE using simple algorithm
    cluster_matrix = zeros(2, 2);
    pilotIndex = simple_algorithm(R,pilotLength,K,L,sortedMaxIndices,C, cluster_matrix, 0);
elseif strcmp("kmeans", flag_alg) == 1
    % Do pilot assignment to each UE using k-means algorithm
    pilotIndex = Kmeansfinal(UE_X, UE_Y, K, pilotLength,R,L,sortedMaxIndices,C);
else
    fprintf("wrong flag indication");
end 

end 
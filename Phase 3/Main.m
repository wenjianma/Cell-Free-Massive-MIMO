%% IK2200 Communication System Design 
% Project Title: Pilot Assignment and Cluster Formation in Cell-Free Massive MIMO Networks
% Team Gyros 

%% Parameter Initilization
% Use this file to generate desirable plots and call all functions. 

% Start measuring time 
tic;
% Number of UEs 
K = 20; 

% Number of APs 
L = 36; 

% Number of APs/Cluster
C = 7; 
% C = 36 (= L) is the case for no clustering.  (no clustering case)
% C = 7 Choose maximum 7 APs in one cluster under simple clustering formation algorithm. But we will remove the AP that can't provide good service.(simple clustering case)
% C = sqrt(L) if we choose predetermined clustering formation (predetermined clustering case)

% Number of antennas per AP
N = 4; 

% prelogFactor
% in our scenario, we have tp = 10, tc = 200, td = 200-10 = 190;so prelogfactor should be 0.95
prelogFactor = 0.95; 

% Length of pilot sequence which means the number of orthogonal pilots we have 
pilotLength = 10;

%number of channel realization per Monte-Carlo simulation
numOfRealization = 10;

%number of Monte-Carlo simulations
numOfSim = 10; 

%Uplink transmission power from UE to AP in mW
p = 100;  

%Downlink max transmission power from AP to UE in mW
Pmax = 1000;

% Select the channel type
flag_channel = ["correlatedRayleigh", "uncorrelatedRayleigh"];
flag_ch = flag_channel(1);

% Select the algorithm
flag_algorithm = ["random","simple","kmeans"];
flag_alg = flag_algorithm(2);

% Select AP-UE association 
flag_serve = ["all","predetermined","cluster"];
flag_serv = flag_serve(3);
%%
% initilization of SE results
SE_MR_all = zeros(K,1,numOfSim);
SE_RZF_all = zeros(K,1,numOfSim);
SINR_MR_all = zeros(K,1,numOfSim);
SINR_RZF_all = zeros(K,1,numOfSim);
%% Simulation Part
for n = 1:numOfSim
    %% General Setup Generation
    % generate the general setup with AP and UE locations 
    % Output: spatial correlation and pilot index
    % Inside generalsetup, we have pilot assignment & cluster formation block
    [R,bkl,pilotIndex,sortedMaxIndices]= generalsetup(K,L,N,pilotLength,flag_ch, flag_alg, flag_serv ,C);

    %% Channel Estimation
    % Output: estimated channel matrix 
    [H,Hhat] = functionChannelEstimates(R,numOfRealization,L,K,N,pilotLength,pilotIndex,p); 

    %% Precoder block   
    % Output: precoder vector that AP l assigns to UE i
    [w_MR, w_RZF] = MR_RZF_Precoder(Hhat,p,numOfRealization,N,K,L,sortedMaxIndices,C);
   
    %% Heuristic Power Allocation
    % Output: power allocation in the downlink
    power_alloc = Heuristic_Power_Allocation(Pmax,bkl,L,K, sortedMaxIndices);
   
    %% Spectral Efficiency calculation
    % Output: Spectral Efficiency using Maximum Ratio (MR) Precoding
    [SE_MR, SINR_MR] = SE_calculation(prelogFactor,w_MR,power_alloc,H,numOfRealization,K,L,N, sortedMaxIndices, C);
    SE_MR_all(:,:,n) = SE_MR;
    SINR_MR_all(:,:,n) = SINR_MR;
    
    % Output: Spectral Efficiency using Regularized Zero Forcing (RZF) Precoding
    [SE_RZF, SINR_RZF] = SE_calculation(prelogFactor,w_RZF,power_alloc,H,numOfRealization,K,L, N, sortedMaxIndices, C);
    SE_RZF_all(:,:,n) = SE_RZF;
    SINR_RZF_all(:,:,n) = SINR_RZF;

end 

%% Save elapsed time
elapsedTime = toc; 
RunTimefilename = sprintf('runtime_%dUES%dAPs_%s_%s_%s.mat', K, L, flag_alg, flag_ch, flag_serv);
save(RunTimefilename, "elapsedTime");
movefile(RunTimefilename,"Results/");


%% Save results to plot
% Create a temporary structure for SE_MR_all and SINR_MR_all
tempStructMR = struct('SE_MR', SE_MR_all, 'SINR_MR', SINR_MR_all);
tempStructRZF = struct('SE_RZF', SE_RZF_all, 'SINR_RZF', SINR_RZF_all);

% Save the struct for SE_MR_all and SINR_MR_all
save(sprintf("MR_all_%dUES%dAPs_%s_%s_%s.mat", K, L, flag_alg, flag_ch, flag_serv), '-struct', 'tempStructMR')
save(sprintf("RZF_all_%dUES%dAPs_%s_%s_%s.mat", K, L, flag_alg, flag_ch, flag_serv), '-struct', 'tempStructRZF')


%% Generate Plots 
run_plots(K, L, numOfSim, flag_alg, flag_ch, flag_serv);

%% Move result files to Results folder
movefile(sprintf("MR_all_%dUES%dAPs_%s_%s_%s.mat",K,L,flag_alg,flag_ch, flag_serv),"Results/");
movefile(sprintf("RZF_all_%dUES%dAPs_%s_%s_%s.mat",K,L,flag_alg,flag_ch, flag_serv),"Results/");



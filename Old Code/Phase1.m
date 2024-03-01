%% IK2200 Communication System Design 
% Project Title: Pilot Assignment and Cluster Formation in Cell-Free Massive MIMO Networks
% Team Gyros 

%% Main File 
% Use this file to generate desirable plots and call all functions. 

% Numbe of UEs 
K = 20; 

% Number of APs 
L = 36; 

% Number of antennas per AP
N = 4; 

% prelogFactor
prelogFactor = 0.95; 
% in our scenario, we have tp = 10, tc = 200, td = 200-10 = 190 
% so prelogfactor should be 0.95

% Length of pilot sequence 
pilotLength = 10;

%number of channel realization
numOfRealization = 10;

%number of Monte-Carlo Simulations
numOfSim = 10; 

%Uplink transmission power from UE to AP in mW
p = 100;  % correct value = 0.1, using large value for testing

%Downlink max transmission power from AP to UE in mW
Pmax = 1000;

% D matrix -> 4 dimensional
%D = repmat(eye(N), [1, 1, K, L]);
D = ones(N,N,K,L);

%spectral efficiency calculation result
SE_MR = zeros(K,1);

% generate the general setup with AP and UE locations and get the spatial
% outputs: correlation matrix, LSF coefficient,pilotindex
flag = ["correlatedRayleigh", "uncorrelatedRayleigh"];
%use flag(1) for correlatedRayleigh 
%use flag(2) for uncorrelatedRayleigh

SE_MR_all = zeros(K,1,numOfSim);
SE_RZF_all = zeros(K,1,numOfSim);
for n = 1:numOfSim
    [R,bkl,pilotIndex,sortedMinIndices]= generalsetup(K,L,N,pilotLength,flag(2)); % note: check bkl -> LSF
    %[R,bkl,pilotIndex]= generalsetup_gyros(K,L,N,pilotLength,flag(1));
    %R = permute(R, [1, 2, 4, 3]);
    
    %% Channel Estimation
    % % channel estimation block
   
    [H,Hhat] = functionChannelEstimates(R,numOfRealization,L,K,N,pilotLength,pilotIndex,p); % compare H and Hhat 

    %Hhat_all(:,:,:,n) = Hhat;
    % % Precoder block   % Output parameters -> w = precoder vector that AP l assigns to UE i
    [w_MR, w_RZF] = MR_ZRF_Precoder(Hhat,p,numOfRealization,N,K,L);
    %[w_MR, w_RZF] = Precoder(Hhat,p,numOfRealization,N,K,L);
   
    %% Heuristic Power Allocation 
    % Heuristic_Power_Allocation
    power_alloc = Heuristic_Power_Allocation(Pmax,bkl,L,K);
    %power_alloc = Pmax/K*ones(K,L);
    
    %% Spectral Efficiency calculation
    %SE_MR = Spectral_Efficiency_Calculation(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L);
    SE_MR = SE_test(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L);
    %[nMR, dMR, SE_MR]= SE_lina2(prelogFactor,w_MR,power_alloc,D,H);
    SE_MR_all(:,:,n) = SE_MR;

    SE_RZF = SE_test(prelogFactor,w_RZF,power_alloc,D,H,numOfRealization,K,L);
    %[nRZF, dRZF, SE_RZF] = SE_lina2(prelogFactor,w_RZF,power_alloc,D,H);
    SE_RZF_all(:,:,n) = SE_RZF;

end 
%% Save results to plot whenever 
%when running random assignmet 
% save("SE_MR_all_20UES64APs_random_Uncor.mat","SE_MR_all")
% save("SE_RZF_all_20UES64APs_random.mat","SE_RZF_all")
% %when running the simple algorithm 
% save("SE_MR_all_20UES64APs_alg_Uncor.mat","SE_MR_all")
% save("SE_RZF_all_20UES64APs_alg.mat","SE_RZF_all")
% %%
% fig = figure();
% cdfplot(SE_MR_all)
% xlabel('Spectral Efficiency per UE [bits/s/Hz]')
% ylabel('CDF')
% title("Phase 1 simulator")
% saveas(fig, "SE_20UEs.png")

%%
% plot(sort(SE_MR_all(:)),linspace(0,1,K*numOfSim),'r-.','LineWidth',2);
% hold on
% % plot(sort(SE_RZF_all(:)),linspace(0,1,K*numOfSim),'b-.','LineWidth',2);
% % xlabel(['Spectral effici' ...
% %     '.ency [bit/s/Hz]'],'Interpreter','Latex');
% % ylabel('CDF','Interpreter','Latex');

%% Debugging notes (Wenjian)

% 1. check generalsetup uncorrelated & correlated (the difference)
% 2. check the theory of RZF precoder





%% IK2200 Communication System Design 
% Project Title: Pilot Assignment and Cluster Formation in Cell-Free Massive MIMO Networks
% Team Gyros 

%% Main File 
% Use this file to generate desirable plots and call all functions. 

% Numbe of UEs 
K = 20; 

% Number of APs 
L = 16; 

% Number of antennas per AP
N = 4; 

% prelogFactor
prelogFactor = 0.9;

% Length of pilot sequence 
pilotLength = 10;

%number of channel realization
numOfRealization = 100;

%number of Monte-Carlo Simulations
numOfSim = 1000; 

%Uplink transmission power from UE to AP in mW
p = 100;  % correct value = 0.1, using large value for testing

%Downlink max transmission power from AP to UE in mW
Pmax = 1000;

% D matrix -> 4 dimensional
D = ones(N,N,K,L);

%spectral efficiency calculation result
SE_MR = zeros(K,1);

% generate the general setup with AP and UE locations and get the spatial
% outputs: correlation matrix, LSF coefficient,pilotindex
flag = ["correlatedRayleigh", "uncorrelatedRayleigh"];
%use flag(1) for correlatedRayleigh 
%use flag(2) for uncorrelatedRayleigh

SE_MR_all = zeros(K,1,numOfSim);
for n = 1:numOfSim
    [R,bkl,pilotIndex]= generalsetup(K,L,N,pilotLength,flag(2)); % note: check bkl -> LSF
    %[R,bkl,pilotIndex]= generalsetup_gyros(K,L,N,pilotLength,flag(2));
    R = permute(R, [1, 2, 4, 3]);
    bkl = bkl';
    %% Channel Estimation
    % % channel estimation block
    [H,Hhat] = functionChannelEstimates(R,numOfRealization,L,K,N,pilotLength,pilotIndex,p); % compare H and Hhat 
    %something wrong with this function!
    %Hhat_all(:,:,:,n) = Hhat;
    % % Precoder block   % Output parameters -> w = precoder vector that AP l assigns to UE i
    [w_MR, w_RZF] = MR_ZRF_Precoder(Hhat,p,numOfRealization,N,K,L);
   
    %% Heuristic Power Allocation 
    % Heuristic_Power_Allocation
    power_alloc = Heuristic_Power_Allocation(Pmax,bkl,L,K);
    
    %% Spectral Efficiency calculation
    %SE_MR = Spectral_Efficiency_Calculation(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L);
    SE_MR = SE_test(prelogFactor,w_MR,power_alloc,D,H,numOfRealization,K,L);
    SE_MR_all(:,:,n) = SE_MR;

end 
%%
figure();
cdfplot(SE_MR_all)
%% 

% flat_data = reshape(SE_MR_all, [K, numOfSim]);
% 
% % Initialize the figure
% figure;
% 
% for user = 1:K
%     user_data = flat_data(user, :);
%     sorted_data = sort(user_data);
%     cdf_data = (1:length(sorted_data)) / length(sorted_data);
% 
%     % Plot the CDF for each user
%     plot(sorted_data, cdf_data, 'DisplayName', ['User ', num2str(user)]);
%     hold on;
% end
% 
% % Add labels and legend
% xlabel('Spectral Efficiency');
% ylabel('CDF');
% title('CDF of Spectral Efficiency per User');
% legend('Location', 'Best');
% grid on;
% 
% % Show the plot
% hold off;







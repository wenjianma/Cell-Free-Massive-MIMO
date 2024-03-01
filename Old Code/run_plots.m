SE_MR_random_Uncor = load("SE_MR_all_20UES36APs_random_Uncor.mat");
SE_MR_random_Cor = load("SE_MR_all_20UES36APs_random_Cor.mat");

% SE_RZF_random = load("SE_RZF_all_20UES64APs_random.mat");
SE_MR_alg_Uncor = load("SE_MR_all_20UES36APs_alg_Uncor.mat");
SE_MR_alg_Cor = load("SE_MR_all_20UES36APs_alg_Cor.mat");

% SE_RZF_alg = load("SE_RZF_all_20UES64APs_alg.mat");

%plot(sort(SE_MR_random_Uncor.SE_MR_all(:)),linspace(0,1,K*numOfSim),'r:','LineWidth',3,'DisplayName','MR Random');
hold on;
plot(sort(SE_MR_random_Cor.SE_MR_all(:)),linspace(0,1,K*numOfSim),'r:','LineWidth',3,'DisplayName','MR Random');
hold on;
%plot(sort(SE_MR_alg_Uncor.SE_MR_all(:)),linspace(0,1,K*numOfSim),'b--','LineWidth',3,'DisplayName','MR Algorithm');
plot(sort(SE_MR_alg_Cor.SE_MR_all(:)),linspace(0,1,K*numOfSim),'b--','LineWidth',3,'DisplayName','MR Algorithm');
legend('Position', [0.7, 0.7, 0.1, 0.1],'FontSize',14)
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
title('36 APs 20 UEs - Correlated Rayleigh','FontSize',14)

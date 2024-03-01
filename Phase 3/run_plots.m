%% General Description
% This function is used to generate the required plots automatically.

%% Function Description
% Input parameters:
% K          =  number of UEs
% L          =  number of APs
% numOfSim   =  number of simulations
% flag_alg   =  flag that decides which algorithm is used
% flag_ch    =  flag that decided which channel is used
% flag_serv  =  flag that decides the UE-AP association

% Output parameters:
% No output parameters. Just generated files and plots that are automatically saved.

%% Function Part

function run_plots(K, L, numOfSim, flag_alg, flag_ch, flag_serv)
    % load data
    files = dir("*.mat");
    for i = 1:length(files)
        fileData = load(files(i).name);
        fieldName = fieldnames(fileData);
        allData.(fieldName{1}) = fileData.(fieldName{1});
    end 
    
    % get variable names from struc
    SE = fieldnames(allData);
    fig = figure();
    plot(sort(allData.(SE{1})(:)),linspace(0,1,K*numOfSim),'r:','LineWidth',2,'DisplayName','MR');
    hold on;
    plot(sort(allData.(SE{2})(:)),linspace(0,1,K*numOfSim),'b--','LineWidth',2,'DisplayName','RZF');
    legend('Position', [0.7, 0.7, 0.1, 0.1],'FontSize',10);
    xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
    ylabel('CDF','Interpreter','Latex');
    title(sprintf('%d APs %d UEs - %s %s %s',L,K,flag_alg,flag_ch, flag_serv),'FontSize',14)
    saveas(fig, sprintf('Results/%d APs %d UEs - %s %s %s',L,K,flag_alg,flag_ch, flag_serv),'jpg')
end
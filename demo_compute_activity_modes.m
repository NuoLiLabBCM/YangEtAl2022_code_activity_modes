% This demo script computes the activity modes using 1/2 of the data and
% use independent trials for activity projections

clc
clear all
close all

addpath('./func')

%% load data: 
% this demo uses dataset from Yang et al 2022, which can be downloaded from
% Zenodo.org <http://dx.doi.org/10.5281/zenodo.6846161>
%     see description of dataset structure and variables in the documentation included on Zenoodo



% this loads the dataset, depending on where you store the file on your local computer, the path name will differ.
% load the spike times information specifically
load '..\YangEtAl_2022\Data_CompileData1_YangEtAl22.mat' neuron_spike_times

spk_times_no_correct_ALL = neuron_spike_times(:,1);         % instructed lick left trial, correct trial
spk_times_yes_correct_ALL = neuron_spike_times(:,2);        % instructed lick right trial, correct trial
spk_times_no_error_ALL = neuron_spike_times(:,3);           % instructed lick left trial, error trial
spk_times_yes_error_ALL = neuron_spike_times(:,4);          % instructed lick right trial, error trial

disp('File loaded. This will run for a long time, ~ 1 hour...')

%% define time points
time_epochs = [-2.6 -1.3 0];        % start of sample delay response
start_t = time_epochs(1)-.4;        % trials start (0.4s prior to sample start)


%% Sort trials into training and testing data, compute PSTH
PSTH_yes_correct_training = [];
PSTH_no_correct_training = [];
PSTH_yes_error_training = [];
PSTH_no_error_training = [];

PSTH_yes_correct_test = [];
PSTH_no_correct_test = [];
PSTH_yes_error_test = [];
PSTH_no_error_test = [];

N_trial_train = [];
N_trial_test = [];
for i_cell = 1:size(spk_times_yes_correct_ALL,1)

    if rem(i_cell,1000)==0
        disp(['Computing PSTHs, processing cell ',num2str(i_cell)]);
    end
    
    % yes correct trials
    spk_times_tmp = spk_times_yes_correct_ALL{i_cell,1};
    i_sort = 1:length(spk_times_tmp);
    iTraining = randsample(i_sort,round(length(spk_times_tmp)/2));
    iTest = i_sort(~ismember(i_sort,iTraining));

    [psth0 t] = func_getPSTH(spk_times_tmp(iTraining),-3.5,2);
    [psth1 t] = func_getPSTH(spk_times_tmp(iTest),-3.5,2);
    
    PSTH_yes_correct_training(i_cell,:) = psth0;
    PSTH_yes_correct_test(i_cell,:) = psth1;
    
    n_trials_yes_correct_train = length(iTraining);
    n_trials_yes_correct_test = length(iTest);

    
    % no correct trials
    spk_times_tmp = spk_times_no_correct_ALL{i_cell,1};
    i_sort = 1:length(spk_times_tmp);
    iTraining = randsample(i_sort,round(length(spk_times_tmp)/2));
    iTest = i_sort(~ismember(i_sort,iTraining));

    [psth0 t] = func_getPSTH(spk_times_tmp(iTraining),-3.5,2);
    [psth1 t] = func_getPSTH(spk_times_tmp(iTest),-3.5,2);
    
    PSTH_no_correct_training(i_cell,:) = psth0;
    PSTH_no_correct_test(i_cell,:) = psth1;

    n_trials_no_correct_train = length(iTraining);
    n_trials_no_correct_test = length(iTest);

    
    % yes error trials
    spk_times_tmp = spk_times_yes_error_ALL{i_cell,1};
    if ~isempty(spk_times_tmp)
        i_sort = 1:length(spk_times_tmp);
        iTraining = randsample(i_sort,round(length(spk_times_tmp)/2));
        iTest = i_sort(~ismember(i_sort,iTraining));
        
        [psth0 t] = func_getPSTH(spk_times_tmp(iTraining),-3.5,2);
        [psth1 t] = func_getPSTH(spk_times_tmp(iTest),-3.5,2);
        
        n_trials_yes_error_train = length(iTraining);
        n_trials_yes_error_test = length(iTest);
    else
        psth0 = nan(1,5101);
        psth1 = nan(1,5101);

        n_trials_yes_error_train = 0;
        n_trials_yes_error_test = 0;
    end
    
    PSTH_yes_error_training(i_cell,:) = psth0;
    PSTH_yes_error_test(i_cell,:) = psth1;
    

    
    % no error trials
    spk_times_tmp = spk_times_no_error_ALL{i_cell,1};
    if ~isempty(spk_times_tmp)
        i_sort = 1:length(spk_times_tmp);
        iTraining = randsample(i_sort,round(length(spk_times_tmp)/2));
        iTest = i_sort(~ismember(i_sort,iTraining));
        
        [psth0 t] = func_getPSTH(spk_times_tmp(iTraining),-3.5,2);
        [psth1 t] = func_getPSTH(spk_times_tmp(iTest),-3.5,2);

        n_trials_no_error_train = length(iTraining);
        n_trials_no_error_test = length(iTest);
        
    else
        psth0 = nan(1,5101);
        psth1 = nan(1,5101);
        
        n_trials_no_error_train = 0;
        n_trials_no_error_test = 0;
    end
    PSTH_no_error_training(i_cell,:) = psth0;
    PSTH_no_error_test(i_cell,:) = psth1;

    
    
    N_trial_train = [N_trial_train;...
        n_trials_yes_correct_train,...
        n_trials_no_correct_train,...
        n_trials_yes_error_train,...
        n_trials_no_error_train,...
        ];
    
    
    N_trial_test = [N_trial_test;...
        n_trials_yes_correct_test,...
        n_trials_no_correct_test,...
        n_trials_yes_error_test,...
        n_trials_no_error_test,...
        ];

       
end




%% compute activity modes
i_pts = find(t>start_t);     % cut off the trace before -3, some dataset has a clipping artifact
i_sel = min(N_trial_test,[],2)>=2;

PSTH_yes_correct_sel = PSTH_yes_correct_training(i_sel,i_pts);
PSTH_no_correct_sel = PSTH_no_correct_training(i_sel,i_pts);
PSTH_yes_error_sel = PSTH_yes_error_training(i_sel,i_pts);
PSTH_no_error_sel = PSTH_no_error_training(i_sel,i_pts);

PSTH_yes_correct_sel_test = PSTH_yes_correct_test(i_sel,i_pts);
PSTH_no_correct_sel_test = PSTH_no_correct_test(i_sel,i_pts);
PSTH_yes_error_sel_test = PSTH_yes_error_test(i_sel,i_pts);
PSTH_no_error_sel_test = PSTH_no_error_test(i_sel,i_pts);

T_cue_aligned_sel = t(i_pts);


[orthonormal_basis var_allDim] = func_compute_activity_modes_DRT(PSTH_yes_correct_sel, PSTH_no_correct_sel, PSTH_yes_error_sel, PSTH_no_error_sel, T_cue_aligned_sel, time_epochs);


save Test_Dataset_Yang_et_al_2021_activity_modes PSTH* T_cue_aligned_sel orthonormal_basis var_allDim i_sel
load Test_Dataset_Yang_et_al_2021_activity_modes PSTH* T_cue_aligned_sel orthonormal_basis var_allDim i_sel


%% plot activity modes
% correct trials, independent trials
activityRL_train = [PSTH_yes_correct_sel PSTH_no_correct_sel PSTH_yes_error_sel PSTH_no_error_sel];
activityRL_test = [PSTH_yes_correct_sel_test PSTH_no_correct_sel_test];
activityRL_test = activityRL_test-repmat(mean(activityRL_train,2),1,size(activityRL_test,2));        % remove mean (from correct trials)
proj_allDim = activityRL_test'*orthonormal_basis;

% error trials, independent trials
activityRL_train = [PSTH_yes_correct_sel PSTH_no_correct_sel PSTH_yes_error_sel PSTH_no_error_sel];
activityRL_err = [PSTH_yes_error_sel_test PSTH_no_error_sel_test];
activityRL_err = activityRL_err-repmat(mean(activityRL_train,2),1,size(activityRL_test,2));        % remove mean (from correct trials)
proj_allDim_err = activityRL_err'*orthonormal_basis;



figure
bar(var_allDim(1:30));
xlabel('Activity modes')
ylabel('Frac var.');

figure
for i_pc = 1:16
    subplot(4,4,i_pc); hold on
    plot(T_cue_aligned_sel(1,:),proj_allDim(1:size(T_cue_aligned_sel,2),i_pc),'b')
    plot(T_cue_aligned_sel(1,:),proj_allDim(size(T_cue_aligned_sel,2)+1:end,i_pc),'r')

    plot(T_cue_aligned_sel(1,:),proj_allDim_err(1:size(T_cue_aligned_sel,2),i_pc),'color',[.7 .7 1])
    plot(T_cue_aligned_sel(1,:),proj_allDim_err(size(T_cue_aligned_sel,2)+1:end,i_pc),'color',[1 .7 .7])
    
    %ylim([-500 500])
    title(['mode ',num2str(i_pc)]);
end
subplot(4,4,1); ylabel('Activity proj.'); xlabel('Time')





%% plot top behavior-relevant activity modes
mode_ID = [1 2 6 3 7 8 9];
mode_name = {'stimulus', 'choice', 'action', 'outcome', 'ramping', 'go', 'response'};

proj_allDim = activityRL_test'*orthonormal_basis;
var_allDim = sum(proj_allDim.^2);

var_allDim = var_allDim/sum(var_allDim);    

figure
bar(var_allDim(mode_ID));
xlabel('Activity modes')
ylabel('Frac var.');
title(['Total Cross Validated Var Explained: ',num2str(sum(var_allDim(mode_ID)))]);


n_plot = 0;
figure
for i_mode = mode_ID
    n_plot = n_plot+1;
    disp(['ploting mode ',num2str(n_plot)]);
    
    proj_iPC_allBtstrp = [];
    projErr_iPC_allBtstrp = [];
    for i_btstrp = 1:20
        i_sample = randsample(size(activityRL_test,1),size(activityRL_test,1),'true');
        proj_iPC_allBtstrp(i_btstrp,:) = activityRL_test(i_sample,:)'*orthonormal_basis(i_sample,i_mode);
        projErr_iPC_allBtstrp(i_btstrp,:) = activityRL_err(i_sample,:)'*orthonormal_basis(i_sample,i_mode);
    end
    
    subplot(2,4,n_plot); hold on
    func_plot_mean_and_sem(T_cue_aligned_sel(1,:), projErr_iPC_allBtstrp(:,1:size(T_cue_aligned_sel,2)), [.4 .4 1], [.8 .8 1], 'n', 2);
    func_plot_mean_and_sem(T_cue_aligned_sel(1,:), projErr_iPC_allBtstrp(:,size(T_cue_aligned_sel,2)+1:end), [1 .4 .4], [1 .8 .8], 'n', 2);
    
    func_plot_mean_and_sem(T_cue_aligned_sel(1,:), proj_iPC_allBtstrp(:,1:size(T_cue_aligned_sel,2)), 'b', [.6 .6 1], 'n', 2);
    func_plot_mean_and_sem(T_cue_aligned_sel(1,:), proj_iPC_allBtstrp(:,size(T_cue_aligned_sel,2)+1:end), 'r', [1 .6 .6], 'n', 2);
    
    y_scale = mean([proj_iPC_allBtstrp projErr_iPC_allBtstrp]);
    line([-2.6 -2.6],[min(y_scale) max(y_scale)]*1.2,'color','k','linestyle',':');
    line([-1.3 -1.3],[min(y_scale) max(y_scale)]*1.2,'color','k','linestyle',':');
    line([0 0],[min(y_scale) max(y_scale)]*1.2,'color','k','linestyle',':');
    xlim([-3.2 2.2]);
    
    title(['mode ',mode_name{n_plot}]);
end
subplot(2,4,1); ylabel('Activity proj.'); xlabel('Time')





%% visualize activity modes on tSNE
if exist('./Test_Dataset_Yang_et_al_2021_tSNE.mat')
    
    load './Test_Dataset_Yang_et_al_2021_tSNE.mat';
    
    mappedY = mappedX_all(i_cst & i_sel, :);
    orthonormal_basis_sel = orthonormal_basis(i_cst(i_sel), :);
    
    figure
    n_plot = 0;
    for i_mode = mode_ID
        
        n_plot = n_plot+1;
        disp(['ploting mode ',num2str(n_plot),' on tSNE']);
        
        activity_mode = orthonormal_basis_sel(:,i_mode);
        
        my_dot_size = abs(activity_mode)*80+.5;     % selectity amplitude
        my_dot_size(my_dot_size>10)=10;             % cap at 10
        
        subplot(2,4,n_plot); hold on
        for i_tmp = find(activity_mode>0)'
            plot(mappedY(i_tmp,1), mappedY(i_tmp,2),'ow','markerfacecolor','r','markersize',my_dot_size(i_tmp));
        end
        for i_tmp = find(activity_mode<0)'
            plot(mappedY(i_tmp,1), mappedY(i_tmp,2),'ow','markerfacecolor','b','markersize',my_dot_size(i_tmp));
        end
        title(['mode ',mode_name{n_plot}]);
        
    end

end



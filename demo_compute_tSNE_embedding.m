% This demo script computes the activity modes using all the trials
% also see demo_activity_modes_independentTrials.m

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

disp('File loaded. This will run for a long time...')


%% define time points
sample_start_t = -2.6;      % start of sample epoch
delay_start_t = -1.3;       % start of delay epoch
response_start_t = 0;       % start of response epoch

start_t = sample_start_t-.4;                    % trials start (0.4s prior to sample start)
end_t = response_start_t+ 1.8;                  % trial end (1.8s after go cue) 


%% Sort trials into training and testing data, compute PSTH
PSTH_yes_correct_all = [];
PSTH_no_correct_all = [];
PSTH_yes_error_all = [];
PSTH_no_error_all = [];

N_trial_all = [];
for i_cell = 1:size(spk_times_yes_correct_ALL,1)

    if rem(i_cell,1000)==0
        disp(['Computing PSTHs, processing cell ',num2str(i_cell)]);
    end
    
    % yes correct trials
    spk_times_tmp = spk_times_yes_correct_ALL{i_cell,1};

    [psth0 t] = func_getPSTH(spk_times_tmp,-3.5,2);
    
    PSTH_yes_correct_all(i_cell,:) = psth0;
    
    n_trials_yes_correct = length(spk_times_tmp);

    
    % no correct trials
    spk_times_tmp = spk_times_no_correct_ALL{i_cell,1};

    [psth0 t] = func_getPSTH(spk_times_tmp,-3.5,2);
    
    PSTH_no_correct_all(i_cell,:) = psth0;

    n_trials_no_correct = length(spk_times_tmp);

    
    % yes error trials
    spk_times_tmp = spk_times_yes_error_ALL{i_cell,1};
    if ~isempty(spk_times_tmp)
        [psth0 t] = func_getPSTH(spk_times_tmp,-3.5,2);
        
        n_trials_yes_error = length(spk_times_tmp);
    else
        psth0 = nan(1,5101);

        n_trials_yes_error = 0;
    end
    
    PSTH_yes_error_all(i_cell,:) = psth0;
    

    
    % no error trials
    spk_times_tmp = spk_times_no_error_ALL{i_cell,1};
    if ~isempty(spk_times_tmp)
        [psth0 t] = func_getPSTH(spk_times_tmp,-3.5,2);

        n_trials_no_error = length(spk_times_tmp);
        
    else
        psth0 = nan(1,5101);
        psth1 = nan(1,5101);
        
        n_trials_no_error = 0;
    end
    PSTH_no_error_all(i_cell,:) = psth0;

    
    
    N_trial_all = [N_trial_all;...
        n_trials_yes_correct,...
        n_trials_no_correct,...
        n_trials_yes_error,...
        n_trials_no_error,...
        ];
       
end



%% 
pearson_corr = func_get_corr(spk_times_yes_correct_ALL, spk_times_no_correct_ALL, start_t, end_t);
i_cst = pearson_corr > 0.5;

PSTH_yes_correct_sel = PSTH_yes_correct_all(pearson_corr > 0.5, :);
PSTH_no_correct_sel = PSTH_no_correct_all(pearson_corr > 0.5, :);
T_cue_aligned_sel = t;
N_trials_sel = N_trial_all(pearson_corr > 0.5,:);


%% Down Sampling
bin_size = 50;  % 0.05 sec

% chop out the (start_t,end_t) fragment %  Debug: the right end never reaches 2
i_trial = find(T_cue_aligned_sel(1,:)>=start_t & T_cue_aligned_sel(1,:)<=end_t);
T_cue_aligned_sel = T_cue_aligned_sel(:,i_trial);
PSTH_yes_correct_sel = PSTH_yes_correct_sel(:,i_trial);
PSTH_no_correct_sel = PSTH_no_correct_sel(:,i_trial);


% down sample
T_cue_aligned_tmp = T_cue_aligned_sel(:,bin_size/2:bin_size:end);
PSTH_yes_correct_tmp = PSTH_yes_correct_sel(:,bin_size/2:bin_size:end);
PSTH_no_correct_tmp = PSTH_no_correct_sel(:,bin_size/2:bin_size:end);



%% Normalization
i_baseline = find(T_cue_aligned_tmp(1,:)>start_t & T_cue_aligned_tmp(1,:)<(sample_start_t-.1));
FR_baseline = mean([PSTH_yes_correct_tmp(:,i_baseline) PSTH_no_correct_tmp(:,i_baseline)],2);
FR_baseline = repmat(FR_baseline,1,size(PSTH_yes_correct_tmp,2));
R = PSTH_yes_correct_tmp-FR_baseline;
L = PSTH_no_correct_tmp-FR_baseline;
for i=1:size(R,1)
    norm_tmp = norm([R(i,:) L(i,:)]);
    R(i,:)=R(i,:)/norm_tmp;
    L(i,:)=L(i,:)/norm_tmp;
end
activityRL = [R L];




%% run tSNE 10 times, pick the run with the lowest loss
mappedX_allRun = {};
loss_tSNE_allRun = [];
for i_run = 1:10
    [mappedX_iRun loss_tSNE_iRun]=tsne(activityRL,'Algorithm','exact','Distance','cosine','NumDimensions',2,'NumPCAComponents',50,'Perplexity',50,'Verbose',2);
    mappedX_allRun{i_run,1} = mappedX_iRun;
    loss_tSNE_allRun(i_run,1) = loss_tSNE_iRun;
end
[dummy i_best_run] = min(loss_tSNE_allRun);
mappedX = mappedX_allRun{i_best_run};
clear mappedX_allRun  mappedX_allRun  mappedX_iRun  loss_tSNE_allRun  dummy i_best_run


mappedX_all = nan(size(spk_times_yes_correct_ALL,1),2);
mappedX_all(i_cst,1) = mappedX(:,1);
mappedX_all(i_cst,2) = mappedX(:,2);

figure; hold on
scatter(mappedX_all(:,1), mappedX_all(:,2), 3, [0.7 0.7 0.7],'filled');
xlabel('Dim 1');
ylabel('Dim 2');
title('t-SNE');


save Test_Dataset_Yang_et_al_2021_tSNE mappedX_all i_cst


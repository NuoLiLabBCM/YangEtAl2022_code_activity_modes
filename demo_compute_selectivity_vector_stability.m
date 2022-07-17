% This demo script computes the similarity of stimulus vs. choice
% selectivity vector across time

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
time_epochs = [-2.6 -1.3 0];        % start of sample delay response
start_t = time_epochs(1)-.4;        % trials start (0.4s prior to sample start)


%% compute PSTH on split half data
PSTH_yes_correct_set1 = [];
PSTH_no_correct_set1 = [];
PSTH_yes_error_set1 = [];
PSTH_no_error_set1 = [];

PSTH_yes_correct_set2 = [];
PSTH_no_correct_set2 = [];
PSTH_yes_error_set2 = [];
PSTH_no_error_set2 = [];

N_trial_all = [];
for i_cell = 1:size(spk_times_yes_correct_ALL,1)

    if rem(i_cell,1000)==0
        disp(['process cell ',num2str(i_cell)]);
    end
    
    % yes correct trials
    spk_times_tmp = spk_times_yes_correct_ALL{i_cell,1};
    n_trials_yes_correct = length(spk_times_tmp);
    i_sort = 1:length(spk_times_tmp);
    iset1 = randsample(i_sort,round(length(spk_times_tmp)/2));
    iset2 = i_sort(~ismember(i_sort,iset1));

    [psth0 t] = func_getPSTH(spk_times_tmp(iset1),-3.5,2);
    [psth1 t] = func_getPSTH(spk_times_tmp(iset2),-3.5,2);
    
    PSTH_yes_correct_set1(i_cell,:) = psth0;
    PSTH_yes_correct_set2(i_cell,:) = psth1;
    
    
    
    % no correct trials
    spk_times_tmp = spk_times_no_correct_ALL{i_cell,1};
    n_trials_no_correct = length(spk_times_tmp);
    i_sort = 1:length(spk_times_tmp);
    iset1 = randsample(i_sort,round(length(spk_times_tmp)/2));
    iset2 = i_sort(~ismember(i_sort,iset1));

    [psth0 t] = func_getPSTH(spk_times_tmp(iset1),-3.5,2);
    [psth1 t] = func_getPSTH(spk_times_tmp(iset2),-3.5,2);
    
    PSTH_no_correct_set1(i_cell,:) = psth0;
    PSTH_no_correct_set2(i_cell,:) = psth1;

    
    
    % yes error trials
    spk_times_tmp = spk_times_yes_error_ALL{i_cell,1};
    n_trials_yes_error = length(spk_times_tmp);
    if ~isempty(spk_times_tmp)
        i_sort = 1:length(spk_times_tmp);
        iset1 = randsample(i_sort,round(length(spk_times_tmp)/2));
        iset2 = i_sort(~ismember(i_sort,iset1));
        
        [psth0 t] = func_getPSTH(spk_times_tmp(iset1),-3.5,2);
        [psth1 t] = func_getPSTH(spk_times_tmp(iset2),-3.5,2);
        
    else
        psth0 = nan(1,5101);
        psth1 = nan(1,5101);

    end
    
    PSTH_yes_error_set1(i_cell,:) = psth0;
    PSTH_yes_error_set2(i_cell,:) = psth1;
    

    
    % no error trials
    spk_times_tmp = spk_times_no_error_ALL{i_cell,1};
    n_trials_no_error = length(spk_times_tmp);
    if ~isempty(spk_times_tmp)
        i_sort = 1:length(spk_times_tmp);
        iset1 = randsample(i_sort,round(length(spk_times_tmp)/2));
        iset2 = i_sort(~ismember(i_sort,iset1));
        
        [psth0 t] = func_getPSTH(spk_times_tmp(iset1),-3.5,2);
        [psth1 t] = func_getPSTH(spk_times_tmp(iset2),-3.5,2);

    else
        psth0 = nan(1,5101);
        psth1 = nan(1,5101);
        
    end
    PSTH_no_error_set1(i_cell,:) = psth0;
    PSTH_no_error_set2(i_cell,:) = psth1;

    
    % number of trials
    N_trial_all = [N_trial_all;...
        n_trials_yes_correct,...
        n_trials_no_correct,...
        n_trials_yes_error,...
        n_trials_no_error,...
        ];
       
           
end

i_pts = find(t>-3);     % cut off the trace before -3, some dataset has a clipping artifact
i_sel = min(N_trial_all(:,1:4),[],2)>=10;


T_cue_aligned_sel = t(i_pts);



%% only use units with enough error trials
i_pts = find(t>start_t);     % cut off the trace before -3, some dataset has a clipping artifact
i_sel = min(N_trial_all(:,1:4),[],2)>=6;

PSTH_yes_correct_sel1 = PSTH_yes_correct_set1(i_sel,i_pts);
PSTH_no_correct_sel1 = PSTH_no_correct_set1(i_sel,i_pts);
PSTH_yes_error_sel1 = PSTH_yes_error_set1(i_sel,i_pts);
PSTH_no_error_sel1 = PSTH_no_error_set1(i_sel,i_pts);

PSTH_yes_correct_sel2 = PSTH_yes_correct_set2(i_sel,i_pts);
PSTH_no_correct_sel2 = PSTH_no_correct_set2(i_sel,i_pts);
PSTH_yes_error_sel2 = PSTH_yes_error_set2(i_sel,i_pts);
PSTH_no_error_sel2 = PSTH_no_error_set2(i_sel,i_pts);

T_cue_aligned_sel = t(i_pts);


save Test_Dataset_Yang_et_al_2021_selectivity_vector PSTH* T_cue_aligned_sel i_sel time_epochs
load Test_Dataset_Yang_et_al_2021_selectivity_vector PSTH* T_cue_aligned_sel i_sel time_epochs


%% plot selectivity vector correlation
i_baseline = find(T_cue_aligned_sel(1,:)<time_epochs(1));

i_sample = find(T_cue_aligned_sel(1,:)>=time_epochs(1));
i_delay = find(T_cue_aligned_sel(1,:)>=time_epochs(2));
i_response = find(T_cue_aligned_sel(1,:)>=time_epochs(3));
i_sample = i_sample(1);
i_delay = i_delay(1);
i_response = i_response(1);


% stimulus selective dimension
selectivity_stim1 = (PSTH_yes_correct_sel1+PSTH_yes_error_sel1)/2-(PSTH_no_correct_sel1+PSTH_no_error_sel1)/2;
selectivity_stim2 = (PSTH_yes_correct_sel2+PSTH_yes_error_sel2)/2-(PSTH_no_correct_sel2+PSTH_no_error_sel2)/2;


figure
subplot(1,2,1);
x = corr(selectivity_stim1,selectivity_stim2);
x(x<0)=0;
x(x>.45)=0.45;
imagesc(x); hold on
line([0 4800],[i_sample i_sample],'LineStyle',':','color','k')
line([0 4800],[i_delay i_delay],'LineStyle',':','color','k')
line([0 4800],[i_response i_response],'LineStyle',':','color','k')
line([i_sample i_sample],[0 4800],'LineStyle',':','color','k')
line([i_delay i_delay],[0 4800],'LineStyle',':','color','k')
line([i_response i_response],[0 4800],'LineStyle',':','color','k')
title('stim sel')
colormap parula


% choice selective dimension
selectivity_choice1 = (PSTH_yes_correct_sel1+PSTH_no_error_sel1)/2-(PSTH_no_correct_sel1+PSTH_yes_error_sel1)/2;
selectivity_choice2 = (PSTH_yes_correct_sel2+PSTH_no_error_sel2)/2-(PSTH_no_correct_sel2+PSTH_yes_error_sel2)/2;

subplot(1,2,2);
x = corr(selectivity_choice1, selectivity_choice2);
x(x<0)=0;
x(x>.7)=.7;
imagesc(x); hold on
line([0 4800],[i_sample i_sample],'LineStyle',':','color','k')
line([0 4800],[i_delay i_delay],'LineStyle',':','color','k')
line([0 4800],[i_response i_response],'LineStyle',':','color','k')
line([i_sample i_sample],[0 4800],'LineStyle',':','color','k')
line([i_delay i_delay],[0 4800],'LineStyle',':','color','k')
line([i_response i_response],[0 4800],'LineStyle',':','color','k')
title('choice sel')
colormap parula



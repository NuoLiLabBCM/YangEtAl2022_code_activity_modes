
% Quantifies activity variance (across activity mdoes) for each subpopulation 
% To run this script, you must first run "demo_compute_activity_modes.m"

clc
clear all
close all

addpath('./func')

load Test_Dataset_Yang_et_al_2021_activity_modes PSTH* T_cue_aligned_sel orthonormal_basis

%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% orthonormal_basis - Neuron x modes  
%    - this should be the output from function "func_compute_activity_modes_DRT(...)"
% PSTH_yes_correct  - Neuron x timestep
% PSTH_no_correct   - Neuron x timestep
% T_cue_aligned_sel - Neuron x timestep
% i_cluster         - contains the cluster index for each subplot population
% 
% Fill in the i_cluster below. In this example, the population is randomly 
% split into four, so each will carry roughly a quarter of the variance

N = size(PSTH_yes_correct_sel,1);
i_cluster = zeros(N,1);

n_sub = round(N/4);
i_cluster(1:n_sub) = 1;
i_cluster((n_sub+1):(n_sub*2)) = 2;
i_cluster((2*n_sub+1):(n_sub*3)) = 3;
i_cluster((3*n_sub+1):end) = 4;
i_cluster = randsample(i_cluster, length(i_cluster));

%%%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cut off the trace before -3, some dataset has a clipping artifact
i_pts = find(T_cue_aligned_sel(1,:)>-3);     
PSTH_yes_correct_sel = PSTH_yes_correct_sel(:,i_pts);
PSTH_no_correct_sel = PSTH_no_correct_sel(:,i_pts);
T_cue_aligned_sel = T_cue_aligned_sel(:,i_pts);


%% variance calculation
activityRL = [PSTH_yes_correct_sel PSTH_no_correct_sel];
activityRL = activityRL-repmat(mean(activityRL,2),1,size(activityRL,2));        % remove mean
% activityRL = activityRL./repmat(var(activityRL,0,2),1,size(activityRL,2));      % uni-variance

proj_allDim = activityRL'*orthonormal_basis;
var_allDim = sum(proj_allDim.^2);
total_var = sum(var_allDim);

mode_ID = [1 2 6 3 7 8 9];
mode_name = {'stimulus', 'choice', 'action', 'outcome', 'ramping', 'go', 'response'};

figure
bar(var_allDim(mode_ID)/total_var);
xlabel('Activity modes')
ylabel('Frac var.');
title(['Total Var Explained: ',num2str(sum(var_allDim(mode_ID))/total_var)]);


% identify clusters
clustID_sorted = unique(i_cluster);


% calculate activity variance carried by each cluster
wClust_var_allDim = [];         % fraction of total variance carried by each cluster
wClust_var_eachDim = [];        % fraction of each activity mode (quantified as var) carried by each cluster
for i_clust = 1:length(clustID_sorted)
    
    clust_ID = clustID_sorted(i_clust);
    
    i_sel = i_cluster==clust_ID;
    
    proj_allDim = activityRL(i_sel,:)'*orthonormal_basis(i_sel,:);
    
    var_allDim = sum(proj_allDim.^2);
    wClust_var_allDim(i_clust,1) = sum(var_allDim)/total_var;      % express as fraction of var. of the full population from above
    

    mode_ID = [1 2 6 3 7 8 9];        
    wClust_var_eachDim(i_clust,:) = var_allDim(mode_ID);        % express as var. within this poulation and within each of the top 7 activity mode, here using L2 norm(),  activity mode x the sum of several components (i.e. wClusters) x = x1+x2+x3... 

end
wClust_var_eachDim = wClust_var_eachDim./repmat(sum(wClust_var_eachDim),size(wClust_var_eachDim,1),1);


% plot by cluster
mycolor = [{'g'},{'m'},{'b'},{[1 .4 .1]},{[1 .6 .9]},{'c'}];
figure;
for i_clust = 1:length(clustID_sorted)
    subplot(6,1,i_clust);
    h = bar(wClust_var_eachDim(i_clust,:));
    set(h,'facecolor',mycolor{i_clust});
    xlabel('Activity modes')
    ylabel('Frac var.');
    title(['Total Var of full population Explained: ',num2str(sum(wClust_var_allDim(i_clust)))]);
    ylim([0 1])
    
end

% plot by mode
mode_name = {'stim','choice','action','outcome','ramp','go','response'};
figure;
for i_mode = 1:size(wClust_var_eachDim,2)
    subplot(2,4,i_mode); hold on
    bar(wClust_var_eachDim(:, i_mode),'k');
    xlabel('cluster');
    ylabel('frac var');
    title(mode_name{i_mode})
    ylim([0 1])
end





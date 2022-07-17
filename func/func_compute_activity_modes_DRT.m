function [orthonormal_basis var_allDim] = func_compute_activity_modes_DRT(PSTH_yes_correct, PSTH_no_correct, PSTH_yes_error, PSTH_no_error, T_cue_aligned_sel, time_epochs)

%
% This function takes in PSTHs of each trial type and computes activity modes
%   PSTH_yes_correct -- correct lick right trial
%   PSTH_no_correct -- correct lick left trial
%   PSTH_yes_error -- error lick right trial
%   PSTH_no_error -- error lick left trial
%   T_cue_aligned_sel -- time stamps
%   time_epochs -- [Onset time of Sample, Delay, and Response epoches]
%
%   orthonormal_basis -- neuron weights for each activity modes. 
%   var_allDim -- activity variance explained by each dimension
%

t_sample = time_epochs(1);
t_delay = time_epochs(2);
t_response = time_epochs(3);

activityRL = [PSTH_yes_correct PSTH_no_correct];
activityRL = activityRL-repmat(mean(activityRL,2),1,size(activityRL,2));        % remove mean
[u s v] = svd(activityRL');
proj_allDim = activityRL'*v;

var_s = diag(s(1:size(proj_allDim,2),:)).^2;        % s is the std, s^2 is equal to "var_allDim = sum(proj_allDim.^2)" 
var_s = var_s/sum(var_s);                           % express as fraction of var.


%% compute behavior-related activity modes
CD_stim_mode = [];
CD_choice_mode = [];
CD_outcome_mode = [];
CD_sample_mode = [];    % wt during the first 400 ms of the sample epoch
CD_delay_mode = [];     % wt during the last 400 ms of the delay epoch
CD_go_mode = [];        % wt during the first 400 ms of the response epoch
Ramping_mode = [];      % In Figure S1H, ramping direction (RD) was defined as a vector maximally distinguishing the mean activity before the trial onset (0.6 s window) and the mean activity before the Go cue (0.1 s window). 
GoDirction_mode = [];   % To calculate the go direction (GD), we subtracted (rlick-right, t + rlick-left, t)/2 after the Go cue (Tgo < t < Tgo + 0.1 s) from that before the Go cue (Tgo - 0.1 s < t < Tgo), followed by normalization by its own norm. 


wt = (PSTH_yes_correct+PSTH_yes_error)/2 - (PSTH_no_correct+PSTH_no_error)/2;
i_t = find(T_cue_aligned_sel(1,:)>t_sample & T_cue_aligned_sel(1,:)<t_delay);
CD_stim_mode = mean(wt(:,i_t),2);

wt = (PSTH_yes_correct+PSTH_no_error)/2 - (PSTH_no_correct+PSTH_yes_error)/2;
i_t = find(T_cue_aligned_sel(1,:)>t_delay & T_cue_aligned_sel(1,:)<t_response);
CD_choice_mode = mean(wt(:,i_t),2);

wt = (PSTH_yes_correct+PSTH_no_correct)/2 - (PSTH_yes_error + PSTH_no_error)/2;
i_t = find(T_cue_aligned_sel(1,:)>t_response & T_cue_aligned_sel(1,:)<(t_response+1.3));
CD_outcome_mode = mean(wt(:,i_t),2);



wt = PSTH_yes_correct - PSTH_no_correct;

i_t = find(T_cue_aligned_sel(1,:)>(t_sample+.2) & T_cue_aligned_sel(1,:)<(t_sample+.4));
CD_sample_mode = mean(wt(:,i_t),2);

i_t = find(T_cue_aligned_sel(1,:)>(t_response-.3) & T_cue_aligned_sel(1,:)<(t_response-.1));
CD_delay_mode = mean(wt(:,i_t),2);

i_t = find(T_cue_aligned_sel(1,:)>(t_response+.1) & T_cue_aligned_sel(1,:)<(t_response+.3));
CD_go_mode = mean(wt(:,i_t),2);


wt = (PSTH_yes_correct + PSTH_no_correct)/2;

i_t1 = find(T_cue_aligned_sel(1,:)>(t_sample-.3) & T_cue_aligned_sel(1,:)<(t_sample-.1));
i_t2 = find(T_cue_aligned_sel(1,:)>(t_response-.3) & T_cue_aligned_sel(1,:)<(t_response-.1));
Ramping_mode = mean(wt(:,i_t2),2)-mean(wt(:,i_t1),2);

i_t1 = find(T_cue_aligned_sel(1,:)>(t_response-.1) & T_cue_aligned_sel(1,:)<t_response);
i_t2 = find(T_cue_aligned_sel(1,:)>t_response & T_cue_aligned_sel(1,:)<(t_response+0.1));
GoDirction_mode = mean(wt(:,i_t2),2)-mean(wt(:,i_t1),2);


CD_stim_mode    = CD_stim_mode/norm(CD_stim_mode);
CD_choice_mode  = CD_choice_mode/norm(CD_choice_mode);
CD_outcome_mode = CD_outcome_mode/norm(CD_outcome_mode);
CD_sample_mode  = CD_sample_mode/norm(CD_sample_mode);
CD_delay_mode   = CD_delay_mode/norm(CD_delay_mode);
CD_go_mode      = CD_go_mode/norm(CD_go_mode);
Ramping_mode    = Ramping_mode/norm(Ramping_mode);
GoDirction_mode = GoDirction_mode/norm(GoDirction_mode);


orthonormal_basis = Gram_Schmidt_process([CD_stim_mode CD_choice_mode CD_outcome_mode CD_sample_mode CD_delay_mode CD_go_mode Ramping_mode GoDirction_mode v]);

proj_allDim = activityRL'*orthonormal_basis;
var_allDim = sum(proj_allDim.^2);

var_allDim = var_allDim/sum(var_allDim);                           % express as fraction of var.





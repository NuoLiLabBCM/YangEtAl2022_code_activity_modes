 function [pearson_corr] = func_get_corr(spk_times_yes_correct_ALL, spk_times_no_correct_ALL, start_t, end_t)
%func_get_rate_corr gets the firing rate and firing of each unit
%   Input: cell array of units, each unit contains a cell array of trials
%   and each trial contains a list of spiking times within that trial
%   Output: 
%       pearson_corr: randomly choose half of the trials of each unit to
%       compute the PSTH and derive the pearson_corr for these two PSTHs,
%       which indicates the consistency of the unit

% calculate the Pearson's correlation

pearson_corr = NaN(size(spk_times_yes_correct_ALL,1),1);
for i_unit = 1:size(spk_times_yes_correct_ALL,1)
    i_yes_rand_1 = randsample(length(spk_times_yes_correct_ALL{i_unit}), round(.5 * length(spk_times_yes_correct_ALL{i_unit})));
    i_yes_rand_2 = setdiff(1:length(spk_times_yes_correct_ALL{i_unit}), i_yes_rand_1);
    i_no_rand_1 = randsample(length(spk_times_no_correct_ALL{i_unit}), round(.5 * length(spk_times_no_correct_ALL{i_unit})));
    i_no_rand_2 = setdiff(1:length(spk_times_no_correct_ALL{i_unit}), i_no_rand_1);
    PSTH_yes_correct_1 = func_getPSTH(spk_times_yes_correct_ALL{i_unit}(i_yes_rand_1),start_t, end_t);
    PSTH_yes_correct_2 = func_getPSTH(spk_times_yes_correct_ALL{i_unit}(i_yes_rand_2),start_t, end_t);
    PSTH_no_correct_1 = func_getPSTH(spk_times_no_correct_ALL{i_unit}(i_no_rand_1),start_t, end_t);
    PSTH_no_correct_2 = func_getPSTH(spk_times_no_correct_ALL{i_unit}(i_no_rand_2),start_t, end_t);
    activityRL1 = [PSTH_yes_correct_1 PSTH_no_correct_1];
    activityRL2 = [PSTH_yes_correct_2 PSTH_no_correct_2];
    pearson_corr(i_unit) = corr(activityRL1', activityRL2');
end


    
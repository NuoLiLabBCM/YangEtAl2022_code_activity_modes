function [PSTH time] = func_getPSTH_individualTrials(SpikeTimes, PSTH_StartTime, PSTH_EndTime)

% 
% SpikeTimes -- {n_rep,1}
% 

if nargin == 1
    PSTH_StartTime = -.52;
    PSTH_EndTime = 5.020;
end

time = PSTH_StartTime:.001:PSTH_EndTime;


n_rep = size(SpikeTimes,1);
total_counts = 0;
for i_rep = 1:n_rep
    
    [counts] = hist(SpikeTimes{i_rep,1},PSTH_StartTime:0.001:PSTH_EndTime);
    
%     window = ones(1,25)/0.025;
%     window = ones(1,100)/0.025;
    window = ones(1,200)/0.2;
    
    PSTH(i_rep,:) = conv(counts,window,'same');
    
end

% 
time = time(201:end-200);
PSTH = PSTH(:,201:end-200);
% time = time(26:end-25);
% PSTH = PSTH(26:end-25);
% time = time(101:end-100);
% PSTH = PSTH(101:end-100);

return
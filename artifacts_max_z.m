%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This script is created for automatically 
%detecting bad trials based on a max-Z threshold %%%%%%
%%%%%%    funded by ARC-DP [DP170103148]    %%%%%%
%%%%%%  prepared  Paul Sowman October 2018 borrowing heavily from  ft_rejectvisual 
%Copyright (C) 2005-2006, Markus Bauer, Robert Oostenveld
% Copyright (C) 2006-2016, Robert Oostenveld,     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output is z  (z-values - max over samples) , bad_trials (list of bad trial
%indices , cdata (the data with bad trials replaced by NaN
%input data; an ft struct of epochs and a threshold for max-z value. ~ 10 is
%good usually


function [z,bad_trials,cdata]=artifacts_max_z(data,threshold)

nchan = size(data.trial{1},1);
ntrl  = size(data.trial,2);

level = zeros(nchan, ntrl);

runsum = zeros(nchan, 1);
runss  = zeros(nchan, 1);
runnum = 0;

for i=1:ntrl
    data2  = data.trial{i};
    runsum = runsum + sum(data2, 2);
    runss  = runss  + sum(data2.^2, 2);
    runnum = runnum + size(data2, 2);
end
mval = runsum/runnum;
sd   = sqrt(runss/runnum - (runsum./runnum).^2);


ft_progress('init', 'dial', 'Please wait...') 
for i=1:ntrl
    ft_progress(i/ntrl, 'computing metric %d of %d\n', i, ntrl);
    level(:, i) = max( ( data.trial{i}-repmat(mval, 1, size(data.trial{i}, 2)) )./repmat(sd, 1, size(data.trial{i}, 2)) , [], 2);
    
end

z          = level;
bad_trials = find(max(z)>threshold);

for i =1:length(bad_trials)
data.trial{bad_trials(i)} = NaN(size(data.trial{1}));
end

cdata = data;

end
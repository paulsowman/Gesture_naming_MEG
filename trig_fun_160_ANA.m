function [trl, event] = trig_fun_160_ANA(cfg)

trigger_nums = 161:225;
triggers     = ft_read_data(cfg.dataset,'dataformat','yokogawa_con','chanindx',trigger_nums);

triggers(6,:)  = [0 abs(diff(triggers(6,:)))]; %want to return onsets and offsets of photodetector **FIXME** make this an optional case
triggers(7,:)  = 0; %this channel is the audio and will create loads of triggers
triggers(33,:) = 0; %this channel is always high and not used for triggers

hdr = ft_read_header(cfg.dataset,'dataformat','yokogawa_con');

for j=1:size(triggers,1)
    trig_height(j) = max(triggers(j,:)); %Y = prctile(x,42)
end

trig_thresh = 0.25*max(trig_height);
triggers    = triggers>trig_thresh; %Binarise trigger channels

event = [];
list  = trigger_nums;

for k=1:size(triggers,1)
    channel   = list(k);
    trig      = triggers(k,:);
    pad       = trig(1);
    trigshift = 0;
    begsample = 1;
    
    
    for l=find(diff([pad trig(:)'])>0)
        event(end+1).type   = num2str(channel);
        event(end  ).sample = l + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(l+trigshift);      % assign the trigger value just _after_ going up
    end
    
end

%% Fill up the event struct
event2 = nestedSortStruct(event,'sample');
event3 = {event2.type};
event3 = cell2mat(cellfun(@str2num,event3,'un',0))-161;
event3 = arrayfun(@num2str, event3, 'unif', 0);
event4 = repmat(struct('type','trigger','value',1,'sample',1,'time',1,'duration',''), 1, length(event3));

times2 = [event2.sample];

for m=1:length(event3)
    event4(m).time   = times2(m)/hdr.Fs;
    event4(m).sample = times2(m);
    event4(m).value  = cell2mat(event3(m));
end

event = event4;


%%%% HERE ADD IN TRIGGER ADJUSTMENTS FOR CODE_DELAY, PHOTODETECTOR, AND
%%%% REASSIGN CODES BASED ON TASK - this is all custom for specific
%%%% experiment. Maybe better to put this in a separate function and
%%%% keep the standard stuff separate

%code_delay = 0.03; % 30ms if you use a response key

events1 = event;
events1 = events1(logical([1 abs(diff(cell2mat(cellfun(@str2num,{events1.value},'un',0))))]));%gets rid of repeated triggers - deal with flickering PD
events1 = events1(ismember({events1.value},[{'5'} {'18'} {'20'} {'21'} {'26'} {'27'} {'28'} {'29'}])); % and remove triggers we don't need
events2 = {events1.value};

%Recode bad trials from behaviour here?
if ~isempty(cfg.behav_data)
    behav_data = importdata(cfg.behav_data);
    idx_20     = find(ismember(events2,{'20'}));
    idx_21     = find(ismember(events2,{'21'}));
    bad_trial  = find(~behav_data.data(behav_data.data(:,1)==cfg.block,8));
    
    events2(idx_20(bad_trial)) = {'999'};
    events2(idx_21(bad_trial)) = {'999'};
else
end
%pad with zeros
events2                    = horzcat(repmat ({'0'},1,15),events2);

%look back in time 16 events to find previous match status
for s=16:length(events2)
    if ismember(events2(s),{'20'}) & any(ismember(events2(s-15:s-1),[{'27'} {'28'} {'29'}]))
        events2(s) = {strcat(events2{s},cell2mat(events2(s-min(find(fliplr(ismember(events2(s-15:s-1),[{'27'} {'28'} {'29'}])))))))};
    end
end

for s=16:length(events2)
    if ismember(events2(s),{'21'}) & any(ismember(events2(s-15:s-1),[{'27'} {'28'} {'29'}]))
        events2(s) = {strcat(events2{s},cell2mat(events2(s-min(find(fliplr(ismember(events2(s-15:s-1),[{'27'} {'28'} {'29'}])))))))};
    end
end

events2     = events2(16:end); % remove padding
events2     = cell2mat(cellfun(@str2num,events2,'un',0));
event_times = [events1.time]';

trigger1_idx = min(find(ismember(events2,27:29))); % Start trigger sequence when we have triggers that belong to 1st primes i.e. cut off any triggers that are initial voice key photodetector etc
events2      = events2(trigger1_idx:end)';
event_times  = event_times(trigger1_idx:end);

if isempty(find(ismember(events2,5))) %no photodetector
    photodetectorcodes_n = 0;
    fprintf('No photodetector found\n');
else
    photodetectorcodes_n = sum(ismember(events2,5));
end

if photodetectorcodes_n>20 %isalmost(photodetectorcodes_n,gocodes_n,10)
    photodetector_times = event_times(ismember(events2,5));
    go_times            = event_times(find(ismember(events2,5))-1);
    screen_delay        = mode(photodetector_times-go_times);
else
    screen_delay = 0.056;
end

fprintf('screen_delay is %d\n',screen_delay);

events2 = arrayfun(@num2str, events2, 'unif', 0);
events3 = repmat(struct('type','trigger','value',1,'time',1,'sample',1,'duration',''), 1, length(events2));

for m=1:length(events2)
    events3(m).time   = event_times(m);
    events3(m).sample = event_times(m)*hdr.Fs;
    events3(m).value  = cell2mat(events2(m));
    events3(m).type   = events3(m).type;
end

event = events3;

if ~isempty(cfg.pdchan)
    %Now that its clean find all photodetector events and replace the times of
    %the preceding events **FIXME** need to deal with possibility of missing
    %PDs maybe we should replave PD with trigger code so we have corrected
    %and uncorrected
    event_list        = {event.value};
    event_list        = cell2mat(cellfun(@str2num,event_list,'un',0));
    photodetector_idx = find(ismember(event_list,cfg.pdchan));
    
    %this doesn't deal with missing PDs!
    for n = 1:length(photodetector_idx)
        event(photodetector_idx(n)-1).time   = event(photodetector_idx(n)).time; %replace the time and sample of the trigger with that of the following PD trigger
        event(photodetector_idx(n)-1).sample = event(photodetector_idx(n)).sample;
    end
else
end

% search for "trigger" events
value  = {event(find(strcmp('trigger', {event.type}))).value}';
sample = round([event(find(strcmp('trigger', {event.type}))).sample]');

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig = round(cfg.trialdef.poststim * hdr.Fs);

trl = [];
for n = 1:length(value)
    trg1     = value(n);
    trlbegin = sample(n) + pretrig;
    trlend   = sample(n) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end

end


function ANA_pipe(rawfile,behav_data,blocks)

close all

for block=1:blocks
    
    rawfileb = [rawfile,'_B',num2str(block),'.con'];
    
    hdr = ft_read_header(rawfileb,'dataformat','yokogawa_con'); %hdr to get Fs etc.
    
    if ~isempty(find(ismember(hdr.label,'AG160')))
        fprintf('ADULT MEG')
    else
        fprintf('CHILD MEG')
        error('This function is only for adult MEG')
    end
    % trigger events
    eventcodes = {{'targetmatch'},{'2027'};{'cuematch'},{'2127'};{'targetmismatch'},{'2028'};{'cuemismatch'},{'2128'};{'targetneutral'},...
        {'2029'};{'cueneutral'},{'2129'};{'match'},{'27'};{'mismatch'},{'28'};{'neutral'},{'29'};{'response'},{'26'}};
    
    % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
    cfg                      = [];
    cfg.headerfile           = rawfileb;
    cfg.datafile             = rawfileb;
    cfg.trialdef.triallength = Inf;
    cfg.trialdef.ntrials     = 1; % read in all data as a single segment
    %    cfg.trialfun             = ft_trialfun_general;
    
    cfg = ft_definetrial(cfg);
    
    % ft_preprocessing: reads in MEG data
    cfg.continuous = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [1.0 70]; % bandpass filter
    
    alldata = ft_preprocessing(cfg);
    
    % Select channels 1-160 (i.e. MEG data)
    cfg         = [];
    cfg.channel = alldata.label(1:160);
    
    alldata = ft_selectdata(cfg, alldata);
    
    % deal with 50Hz line noise
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = [49.5 50.5];
    
    alldata = ft_preprocessing(cfg, alldata);
    
    %     Create layout file for later + save
    cfg      = [];
    cfg.grad = alldata.grad; % struct containing gradiometer definition
    lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations
    
    %     Define trials using custom trialfun
    cfg                   = [];
    cfg.dataset           = rawfileb;
    cfg.continuous        = 'yes';
    cfg.trialfun          = 'trig_fun_160_ANA';
    cfg.trialdef.prestim  = 1;         % pre-stimulus interval
    cfg.trialdef.poststim = 2;        % post-stimulus interval
    cfg.pdchan            = [];
    cfg.trigchannels      = [165:190];
    cfg.behav_data        = behav_data;
    cfg.block             = block;
    
    trialinfo_b = ft_definetrial(cfg);
    alldata = ft_redefinetrial(trialinfo_b, alldata);
    
    cfg         = [];
    cfg.demean  = 'yes';
    cfg.detrend = 'yes';
    
    alldata = ft_preprocessing(cfg, alldata);
    
    eventnames = eventcodes(:,1); % extract a list of all event names
    eventnames = [eventnames{:}]; % convert into strings
    trialsgone = 0;
    
    % each cycle is one "block" (i.e. one '.con' file)
    for i = 1:length(trialinfo_b)
        for j = 1:length(eventcodes)
            events.(eventnames{j}) = find(strcmp({trialinfo_b(i).event.value}, eventcodes{j,2})); % 9 fields representing the 9 types of events
            % each field contains a list of all events belonging to this type
            % (by matching event code)
            % NB. this is a tmp var, it gets overwritten in each cycle
        end
        if i == 1 % first block only
            for j = 1:length(eventcodes)
                events_b.(eventnames{j}) = events.(eventnames{j})'; % save the lists to a perm var, also transpose each list
            end
        else % all other blocks
            trialsinblock = length(trialinfo_b(i-1).event); % how many "trials" (i.e. events) were identified in previous block
            trialsgone    = trialsgone + trialsinblock; % add this number to the total number of "past" trials
            
            for j = 1:length(eventcodes)
                events_b.(eventnames{j}) = [events_b.(eventnames{j}); events.(eventnames{j})' + trialsgone]; % continue to append to the perm lists (stored in "events_b")
                % in the end, each perm list will contain all events of that type from all blocks
            end
        end
    end
    
    %%%%***THIS not necessary as no NaNs at this point?***
    for g=1:length(alldata.trial)
        good(g) = ~isempty(find(isnan(alldata.trial{g})));
    end
    
    % ft_redefine the "response" trials & calc erf
    cfg           = [];
    cfg.trials    = setdiff(events_b.response,find(good)); % list of "response" events
    cfg.toilim    = [-0.5 0.5];
    cfg.minlength = 'maxperlen';
    response      = ft_redefinetrial(cfg, alldata);
    
    cfg          = [];
    response_erf = ft_timelockanalysis(cfg, response);
    
    %Run ICA on the "response" trials
    disp('About to run ICA using the SVD method')
    cfg           = [];
    cfg.method    = 'svd';
    response_comp = ft_componentanalysis(cfg, response_erf);
    
    %Change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colour = colormap(flipud(brewermap(64, 'RdBu'))); % change the colormap
    
    %Display Components - change layout as needed
    cfg          = [];
    cfg.viewmode = 'component';
    cfg.layout   = lay;
    cfg.channel  = {'svd001','svd002','svd003'};
    cfg.colormap = colour;
    ft_databrowser(cfg, response_comp);
    
    %**unmix here...is that necessary?***
    cfg           = [];
    cfg.component = 1; % reject top x components
    alldata_clean = ft_rejectcomponent(cfg, response_comp, alldata); % reject these comps from all trials
    
    %reject trials Z max > 10
    [z,bad_trials,data_clean] = artifacts_max_z(alldata_clean,10); 
    good_trials_idx = setdiff(1:length(trialinfo_b.event),bad_trials);   
    
    % ft_redefine all other event types & calc erf
    for j = 1:length(eventcodes)
        cfg                 = [];
        cfg.trials          = setdiff(events_b.(eventnames{j}),bad_trials); %setdiff(events_b.response,find(good))
        cfg.minlength       = 'maxperlen';
        tmp.(eventnames{j}) = ft_redefinetrial(cfg, alldata);
        
    end
    
    % do the same for the cleaned data
    for j = 1:length(eventcodes)
        cfg                       = [];
        cfg.trials                = setdiff(events_b.(eventnames{j}),bad_trials);
        cfg.minlength             = 'maxperlen';
        tmp_clean.(eventnames{j}) = ft_redefinetrial(cfg, alldata_clean);
    end
    
    if block==1     
        for j = 1:length(eventcodes)
            trials.(eventnames{j})       = tmp.(eventnames{j});
            trials_clean.(eventnames{j}) = tmp_clean.(eventnames{j});
        end    
    else
        for j = 1:length(eventcodes)
            trials.(eventnames{j})       = ft_appenddata([],trials.(eventnames{j}),tmp.(eventnames{j}));
            trials_clean.(eventnames{j}) = ft_appenddata([],trials_clean.(eventnames{j}),tmp_clean.(eventnames{j}));
        end
    end
    %%%
    %%HERE we concat across blocks 
end
%%end block loop

% calc erf
for j = 1:length(eventcodes)
    cfg                 = [];
    cfg.nanmean         = 'yes';
    erf.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
    erf.(eventnames{j}) = ft_timelockanalysis(cfg, erf.(eventnames{j})); % Do this to create average field and timelock struct
end

% fill in some erf struct
for j = 1:length(eventcodes)
    tmp                        = trials.(eventnames{j}).trial;
    tmp_mean                   = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    erf.(eventnames{j}).avg    = tmp_mean;
    erf.(eventnames{j}).time   = erf.(eventnames{j}).time;
    erf.(eventnames{j}).dimord = 'chan_time';
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    tmp                              = trials_clean.(eventnames{j}).trial;
    tmp_mean                         = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    fprintf([eventnames{j},' has%6.0f good trials\n'], length(tmp));
    good_trials{j}                   = [eventnames{j},',',num2str(length(tmp))];
    erf_clean.(eventnames{j})        = trials_clean.(eventnames{j});
    erf_clean.(eventnames{j}).avg    = tmp_mean;
    erf_clean.(eventnames{j}).time   = erf_clean.(eventnames{j}).time{1};
    erf_clean.(eventnames{j})        = rmfield(erf_clean.(eventnames{j}), 'trial');
    erf_clean.(eventnames{j}).dimord = 'chan_time';
end

% SAVE all relevant variables in the workspace
save([cd '/data.mat'], 'eventcodes', 'lay', 'response_comp', 'erf', 'erf_clean', 'good_trials');

%% FOR PLOTTING (can use this to regen all plots from saved variables)
% required configurations b4 any calls to ft_multiplotER()
cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = lay;
cfg.baseline   = [-0.2 0];
%cfg.baselinetype = 'absolute';

figure('Name','ft_multiplotER: target-match, target-mismatch, target-neutral'); %{'target'},{'20'};{'cue'},{'21'};{'match'},{'27'};{'mismatch'},{'28'};{'neutral'},{'29'};{'response'},{'26'}};
ft_multiplotER(cfg, erf_clean.(eventnames{1}), ...
    erf_clean.(eventnames{3}), erf_clean.(eventnames{5})); %%%FIXME - USE list for variable names
legend([eventcodes{1:2:5, 1}])

% Calc global averages across all sensors (GFP = global field potentials)
cfg        = [];
cfg.method = 'power';
for j = 1:length(eventnames)
    erf_clean_GFP.(eventnames{j}) = ft_globalmeanfield(cfg, erf_clean.(eventnames{j}));
end

% plot global averages for target-locked
cfg = [];
figure('Name','target_GFP'); hold on
for j = 1:2:5
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{1:2:5, 1}])

% plot global averages for cue-locked
cfg = [];
figure('Name','cue_GFP'); hold on
for j = 2:2:6
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{2:2:6, 1}])

trials_cue_active  = ft_appenddata([],trials_clean.(eventnames{2}),trials_clean.(eventnames{4}));
trials_cue_neutral = tmp_clean.(eventnames{6});

cfg                 = [];
cfg.output          = 'pow';
cfg.channel         = 'MEG';
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.foi             = 5:0.5:70;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin       = 5./cfg.foi;
cfg.tapsmofrq       = 0.4 *cfg.foi;
cfg.toi             = -0.5:0.05:2.0;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.pad             = 'nextpow2';
TFRhann_cue_active  = ft_freqanalysis(cfg, trials_cue_active);
TFRhann_cue_neutral = ft_freqanalysis(cfg, trials_cue_neutral);
%
cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relchange';
cfg.zlim         = [-0.3 0];
cfg.xlim         = [-0.5 1];
cfg.ylim         = [5 30];
cfg.channel      = 'all'; % top figure

figure
subplot(2,1,1)
ft_singleplotTFR(cfg, TFRhann_cue_active);
title('Active')
subplot(2,1,2)
ft_singleplotTFR(cfg, TFRhann_cue_neutral);
title('Neutral')

cfg              = [];
cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'relchange';
cfg.zlim         = [-0.3 0];
cfg.xlim         = [-0.5 1];
cfg.ylim         = [30 70];
cfg.channel      = 'all'; % top figure

figure
subplot(2,1,1)
ft_singleplotTFR(cfg, TFRhann_cue_active);
title('Active_gamma')
subplot(2,1,2)
ft_singleplotTFR(cfg, TFRhann_cue_neutral);
title('Neutral_gamma')

save trials_cue_neutral trials_cue_neutral
save trials_cue_active trials_cue_active
save TFRhann_cue_active TFRhann_cue_active
save TFRhann_cue_neutral TFRhann_cue_neutral

end
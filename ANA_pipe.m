% Delta wave ? (0.1 ? 3 Hz)
% Theta wave ? (4 ? 7 Hz)
% Alpha wave ? (8 ? 15 Hz)
% Mu wave ? (7.5 ? 12.5 Hz)
% SMR wave ? (12.5 ? 15.5 Hz)
% Beta wave ? (16 ? 31 Hz)
% Gamma wave ? (32 ? 100 Hz)

%channel for photodetector
%channel for audio recording
%
% % fixation_event1.set_port_code(17);  # intertrial interval
% % video_event.set_port_code(26+congruency[which_stim1[K]]); onset of video primes# channels 88, 89, 89 for congruency w/ the picture that appears next in the trial
% % fixation_event2.set_port_code(18); # trial 100ms SOA
% % image_event.set_port_code(20); # onset arget picture no cue
% % image_event2.set_port_code(21); # onset target picture with cue for naming # micro on channel 86
%
%
%
% cfg=[];
% cfg.dataset='3055_LR_PI174_2018_06_21_B1.con';
% cfg.trialdef.prestim=0.5;
% cfg.trialdef.poststim=1.5;
% [trl, event] = trig_fun_160_basic(cfg);
% eval(cell2sym({event.value}))';

function ANA_pipe(rawfile,blocks)

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
    %rawfile    = '2525_IF_ME126_2017_04_12_B1.con';
    
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
    cfg.bpfreq     = [1.0 40]; % bandpass filter
    
    alldata = ft_preprocessing(cfg);
    
    % Select channels 1-125 (i.e. MEG data)
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
    cfg.trialfun          = 'trig_fun_160_basic';
    cfg.trialdef.prestim  = 1;         % pre-stimulus interval
    cfg.trialdef.poststim = 2;        % post-stimulus interval
    cfg.pdchan            = [];
    cfg.trigchannels      = [165:190];
    
    trialinfo_b = ft_definetrial(cfg);
    
%     eventx        = eval(cell2sym({trialinfo_b.event.value}))
%     trialinfo_b.event = trialinfo_b.event(logical(abs(diff(eventx)))); %clean up the trial struct to get rid of incorrect photodetectors
%     trialinfo_b.trl   = trialinfo_b.trl(logical(abs(diff(eventx))),:);
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
    
    for i=1:length(alldata.trial)
        good(i) = ~isempty(find(isnan(alldata.trial{i})));
    end
    
    % ft_redefine the "response" trials & calc erf
    cfg           = [];
    cfg.trials    = setdiff(events_b.response,find(good)); % list of "response" events
    %cfg.trials    = cfg.trials(1:length(cfg.trials-10)); % deal with possible short trials at end
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
    ft_databrowser(cfg, response_comp)
    
    %**unmix here...is that necessary?***
    cfg           = [];
    cfg.component = 1; % reject top x components
    alldata_clean = ft_rejectcomponent(cfg, response_comp, alldata); % reject these comps from all trials
    
    % Reject Outlier Trials
    % Display visual trial summary to reject deviant trials.
    cfg             = [];
    cfg.method      = 'summary';
    cfg.keepchannel = 'no';
    cfg.keeptrial   = 'nan';%necessary?
    alldata_clean   = ft_rejectvisual(cfg, alldata_clean);
    
    % ft_redefine all other event types & calc erf
    for j = 1:length(eventcodes)
        cfg                 = [];
        cfg.trials          = setdiff(events_b.(eventnames{j}),find(good)); %setdiff(events_b.response,find(good))
        cfg.minlength       = 'maxperlen';
        tmp.(eventnames{j}) = ft_redefinetrial(cfg, alldata);
        
    end
    
    % do the same for the cleaned data
    for j = 1:length(eventcodes)
        cfg                       = [];
        cfg.trials                = setdiff(events_b.(eventnames{j}),find(good));
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


% calc erf?
for j = 1:length(eventcodes)
    cfg                 = [];
    cfg.nanmean         = 'yes';
    erf.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
    erf.(eventnames{j}) = ft_timelockanalysis(cfg, erf.(eventnames{j})); % Do this to create average field and timelock struct
end

% fill in some erf struct?
for j = 1:length(eventcodes)
    tmp                        = trials.(eventnames{j}).trial;
    tmp_mean                   = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    erf.(eventnames{j}).avg    = tmp_mean;
    erf.(eventnames{j}).time   = erf.(eventnames{j}).time;
    %erf.(eventnames{j})        = rmfield(erf.(eventnames{j}), 'trial');
    erf.(eventnames{j}).dimord = 'chan_time';
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    tmp                              = trials_clean.(eventnames{j}).trial;
    tmp_mean                         = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    erf_clean.(eventnames{j})        = trials_clean.(eventnames{j});
    erf_clean.(eventnames{j}).avg    = tmp_mean;
    erf_clean.(eventnames{j}).time   = erf_clean.(eventnames{j}).time{1};
    erf_clean.(eventnames{j})        = rmfield(erf_clean.(eventnames{j}), 'trial');
    erf_clean.(eventnames{j}).dimord = 'chan_time';
end


% SAVE all relevant variables in the workspace
%save([SubjectFolder 'data.mat'], 'eventcodes', 'lay', 'response_comp', 'erf', 'erf_clean');

%% FOR PLOTTING (can use this to regen all plots from saved variables)
% required configurations b4 any calls to ft_multiplotER()

cfg            = [];
cfg.showlabels = 'yes';
cfg.fontsize   = 6;
cfg.layout     = lay;
cfg.baseline   = [-0.2 0];
%cfg.baselinetype = 'absolute';

% all pairwise comparisons btwn raw data (_erf) & after artefact removal (_erf_clean)
% just to see how good the artefact removal is (atm we reject components 1:5, we can adjust this for more/less removal)
% for j = 1:length(eventcodes)
%     figure;
%     ft_multiplotER(cfg, erf.(eventnames{j}), erf_clean.(eventnames{j}));
% end

% % to compare the 3 conds
% % (requires the list of 'cfg' assignments above, if running in console)
% figure('Name','ft_multiplotER: target, cue, match, mismatch, neutral'); %{'target'},{'20'};{'cue'},{'21'};{'match'},{'27'};{'mismatch'},{'28'};{'neutral'},{'29'};{'response'},{'26'}};
% ft_multiplotER(cfg, erf_clean.(eventnames{1}), erf_clean.(eventnames{2}), ...
%     erf_clean.(eventnames{3}), erf_clean.(eventnames{4}), erf_clean.(eventnames{5})); %%%FIXME - USE list for variable names

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

% % plot global averages for cue-locked
% cfg = [];
% figure('Name','GFP'); hold on
% for j = 1:length(eventnames)
%     plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
% end
% legend([eventcodes{1:length(eventnames), 1}])

% plot global averages for target-locked
cfg = [];
figure('Name','target_GFP'); hold on
for j = 1:2:5
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{1:2:5, 1}])

%     cfg           = [];
%     cfg.output    = 'pow';
%     cfg.channel   = 'MEG';
%     cfg.method    = 'mtmconvol';
%     cfg.taper     = 'hanning';
%     cfg.foi       = 12:0.5:30;                         % analysis 2 to 30 Hz in steps of 2 Hz
%     %cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.2;   % length of time window = 0.5 sec
%     cfg.t_ftimwin = 5./cfg.foi;
%     cfg.tapsmofrq = 0.4 *cfg.foi;
%     cfg.toi       = -0.5:0.05:2.0;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
%     cfg.pad       = 'nextpow2';
%     TFRhann       = ft_freqanalysis(cfg, trials_clean.cue);
%


% cfg              = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'relchange';
% cfg.zlim         = 'minzero';
% cfg.channel      = 'all'; % top figure

% figure
% ft_singleplotTFR(cfg, TFRhann);
%
% cfg              = [];
% cfg.feedback     = 'no';
% cfg.method       = 'triangulation';
% cfg.planarmethod = 'sincos';
% cfg.channel      = {'MEG'};
% cfg.trials       = 'all';
% cfg.neighbours   = ft_prepare_neighbours(cfg, trials_clean.cue);
% data_planar      = ft_megplanar(cfg,trials_clean.cue);

%%tf here

% cfg         = [];
% cfg.channel = 'MEG';
% cfg.method  = 'wavelet';
% cfg.width   = 7;
% cfg.output  = 'pow';
% cfg.foi     = 1:0.5:30;
% cfg.toi     = -0.5:0.05:2.0;
% cfg.pad     = 'nextpow2';
% TFRwave     = ft_freqanalysis(cfg, trials_clean.one);

%     cfg               = [];
%     cfg.combinemethod = 'sum';
%     freq_combined     = ft_combineplanar(cfg,TFRwave);
%
%     freq_combined_avg = ft_freqdescriptives([],freq_combined);

% cfg              = [];
% cfg.layout       = lay;
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'relchange';
% cfg.xlim         = [0.1 0.6];
% %cfg.zlim         = 'minzero';
% cfg.ylim         = [15 30];
% cfg.marker       = 'on';

%     figure
%     ft_topoplotTFR(cfg, TFRwave);
%     title('sensors TFR_beta');colorbar;
%
%     cfg.ylim         = [8 15];
%     figure
%     ft_topoplotTFR(cfg, TFRwave);
%     title('sensors TFR_alpha');colorbar;
%
%     cfg.ylim         = [4 7];
%     figure
%     ft_topoplotTFR(cfg, TFRwave);
%     title('sensors TFR_theta');colorbar;

% cfg.ylim = [1 30];
% figure
% ft_multiplotTFR(cfg, TFRwave);
% title('sensors TFR_broadband');colorbar;
%
% cfg              = [];
% cfg.baseline     = [-0.5 -0.1];
% cfg.baselinetype = 'relchange';
% %cfg.zlim         = 'minzero';
% cfg.channel      = 'all'; % top figure
%
% figure
% ft_singleplotTFR(cfg, TFRwave);
% title('average_sensors');colorbar;
%%%%
% save TFR TFRwave
% save erf erf
% save erf_clean erf_clean
% save trials trials
% save trials_clean trials_clean
%%%%
end
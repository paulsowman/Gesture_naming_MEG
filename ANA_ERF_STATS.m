%%%%%%%%%%%%%%%%%%
%SET TIME WINDOW
%%%%%%%%%%%%%%%%%%%

time_1 = -0.25; %in seconds
time_2 = 0.75;

orig               = cd; %need to start in directory that contains subject folders - could hardcode if you prefer
folders            = dir;
ERF_targetmatch    = []; %initialize just in case you run multiple times
ERF_targetmismatch = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop over subject folders and extract the data you want - in this case
%targettmatch vs target mismatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e                  = 1;
for i=1:length(folders)
    if ~isempty(str2num(folders(i).name))
        cd([orig,'/',folders(i).name])
        if exist('data.mat','file')
            load ('data.mat')
            
            N                        = erf_clean.targetmatch.time;
            V                        = [N' repmat(time_1,length(N),1) repmat(time_2,length(N),1)];
            [minValue(1),time_1_idx] = min(abs(V(:,1)-V(:,2)));
            [minValue(2),time_2_idx] = min(abs(V(:,1)-V(:,3)));
            
            eval(['ERF_targetmatch(e,:,:)=erf_clean.targetmatch.avg(:,time_1_idx:time_2_idx);'])%%%%%%%Put into single array
            eval(['ERF_targetmismatch(e,:,:)=erf_clean.targetmismatch.avg(:,time_1_idx:time_2_idx);'])
            e = e+1;
        else
        end
    else
    end
    cd(orig)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set-up for ept-toolbox https://github.com/Mensen/ept_TFCE-matlab/tree/master/TFCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newlocs     = fieldtripchan2eeglab_paul(erf.targetmatch.grad); %grad structure in EEGLAB style
ept_tfce_nb = ept_ChN2(newlocs,1); %Calculate the neighbourhood

Results = ept_TFCE(ERF_targetmatch,ERF_targetmismatch, newlocs, 'flag_ft',0,'flag_save',1,'plots',1,'type','d'); %Do the stats type "d" = dependent - within subjects


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot with EEGvis toolbox https://github.com/behinger/eegvis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fgPlot = [];
cfgPlot.pvalues = Results.P_Values;
cfgPlot.colormap = {{'div','RdYlBu'},{'div','RdBu'},'seq'};
cfgPlot.topoalpha = 0.05; % where to put the significance dots?
cfgPlot.individualcolorscale = 'row'; % different rows have very different interpretation
cfgPlot.time = [0 0.75]; % we will zoom in
plot_topobutter(cat(3,squeeze(mean(Data{1}(:,1:160,:),1))-squeeze(mean(Data{2}(:,1:160,:),1)),Results.TFCE_Obs,Results.P_Values),-.25:0.001:0.75,Info.Electrodes.e_loc(1:160),cfgPlot) %plot rows 1Difference 2t-values 3p-values

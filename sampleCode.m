%% Odean2023sampleCode
%
% This file explains how to access the data from Odean, Sanayei, and
% Shadlen, "Transient oscillations of neural firing rate associated
% with routing of evidence in a perceptual decision", Journal of
% Neuroscience, 2023.
%
% The spike data from the cued attention task (two random dot
% motion patches were displayed following a cue indicating which location
% would have informative motion) are stored as cuedAttentionData.mat
%
% The spike data from the variable location task (one random dot motion
% patch appeared at one of two locations. Which location was used varied
% from trial to trial and was not cued) are stored as
% variableLocationData.mat
%
% Each file contains a table, where every row is a trial. The final column,
% 'SpCell', has all of the spike times for that trial, relative to
% beginning of the recording session. Other columns include task
% parameters and event times (also relative to the beginning of the
% recording session). A full description of each column is included at the
% end of this script, after the sample code.
%
% Contact naomi.neva.odean@gmail.com for assitance.
%

%% Using the cued attention spike data: figure 4
% This section requires the function BigSBigTloop.m
load('cuedAttentionData')

% we start by setting the bin width
binwidth=10;
% and the x-axis of the figures we are going to plot
edges = -200:binwidth:400;

% First we will make panels A and B (example cell Dm49)
% We begin by selecting data for the raster plots
% Select trials with the appropriate dot location, exclude some trials with
% a different study design, and select trials from the cell of interest
rasterConditionDown=cuedAttentionData.dotLoc==0&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==49;
rasterConditionUp=cuedAttentionData.dotLoc==1&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==49;
% Get the spike times for each condition
allSpikesDown=cuedAttentionData.spCell(rasterConditionDown);
allSpikesUp=cuedAttentionData.spCell(rasterConditionUp);
% Get the time to align to (cue onset)
cueDown=cuedAttentionData.cueOn(rasterConditionDown);
cueUp=cuedAttentionData.cueOn(rasterConditionUp);
% subtract the cue onset from the cue offset to get the cue duration
cueDurDown=cuedAttentionData.cueOff(rasterConditionDown)-cueDown;
cueDurUp=cuedAttentionData.cueOff(rasterConditionUp)-cueUp;


% Before we plot the rasters, we will prepare the PSTHs for this cell
% set the reference time, which everything will be aligned to
refTime = cuedAttentionData.cueOn;
% set a minimum time to include in the plots (not very important in this
% data)
tmin= cuedAttentionData.targetsOn-cuedAttentionData.cueOn;
% set a maximum time to include in the plots. This is important because dot
% duration varies between trials and we don't want to include activity from
% later epochs
tmax= cuedAttentionData.dotsOn-refTime;
% Select trials from the cell of interest with the motion patch above the 
% fixation point
relevantTrials=find(cuedAttentionData.dotLoc==1&cuedAttentionData.cellNumber==49);
[bigS,bigT]=BigSBigTloop(relevantTrials,cuedAttentionData,refTime,tmin,tmax);
% Caluclate the number of spikes in each bin
[nsp] = hist(bigS,edges);
% Calculate the number of ms that contributed to each bin
[nms,be] = hist(bigT,edges);
% divide the number of spikes by the number of ms and multiply by 1000 to
% get the spike rate per second
cueUpAll=(1000*nsp./nms)';

% Do the same thing for trials with the cue below fixation
relevantTrials=find(cuedAttentionData.dotLoc==0&cuedAttentionData.cellNumber==49);
[bigS,bigT]=BigSBigTloop(relevantTrials,cuedAttentionData,refTime,tmin,tmax);
[nsp] = hist(bigS,edges); % number of spikes
[nms,be] = hist(bigT,edges); % number of ms.
cueDownAll=(1000*nsp./nms)'; % spike rate

% Now we plot everything, starting with panel A, trials where the cue is
% above the fixation point
figure
subplot(3,2,1)
hold on
% plot the PSTH
plot(edges/1000,cueUpAll,'color',[.6 0 .4],'linewidth',2)
box off
set(gca,'TickDir','out');
% plot the rasters
for iUp=1:length(allSpikesUp)
    rasterUp=allSpikesUp{iUp}-cueUp(iUp);
    ys=ones(length(rasterUp),1)*iUp+100;
    plot(rasterUp,ys,'k.','markersize',1)
    hold on
end
% make the figure look nice
xlim([-.1,.4])
ax=gca;
ax.FontSize=13;
ylabel('Firing rate (sp/s)             Trial','fontsize',16)
xlabel('Time from cue onset (s)','fontsize',16)
box off
line([0 0],[0 250],'linewidth',2,'color','k')
set(gca,'TickDir','out');
yticks([0 50 100 150 200 250])
yticklabels([0 50 0 50 100 150])
ylim([0 275])

% do the same thing for trials with the cue below fixation
subplot(3,2,3)
hold on
plot(edges/1000,cueDownAll,'color',[.6 0 .4],'linewidth',2)
box off
set(gca,'TickDir','out');
for iDown=1:length(allSpikesDown)
    rasterDown=allSpikesDown{iDown}-cueDown(iDown);
    ys=ones(length(rasterDown),1)*iDown+100;
    plot(rasterDown,ys,'k.','markersize',1)
    hold on
end
xlim([-.1,.4])
ax=gca;
ax.FontSize=13;
ylabel('Firing rate (sp/s)             Trial','fontsize',16)
xlabel('Time from cue onset (s)','fontsize',16)
set(gca,'TickDir','out');
line([0 0],[0 262],'linewidth',2,'color','k')
yticks([0 50 100 150 200 250])
yticklabels([0 50 0 50 100 150])
ylim([0 275])

% Now we do the same thing for example cell Dm35

% Get the rasters
rasterConditionDown=cuedAttentionData.dotLoc==0&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==35;
rasterConditionUp=cuedAttentionData.dotLoc==1&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==35;
allSpikesDown=cuedAttentionData.spCell(rasterConditionDown);
allSpikesUp=cuedAttentionData.spCell(rasterConditionUp);
cueDown=cuedAttentionData.cueOn(rasterConditionDown);
cueUp=cuedAttentionData.cueOn(rasterConditionUp);
cueDurDown=cuedAttentionData.cueOff(rasterConditionDown)-cueDown;
cueDurUp=cuedAttentionData.cueOff(rasterConditionUp)-cueUp;

% Make the PSTHs
refTime = cuedAttentionData.cueOn;
tmin= cuedAttentionData.targetsOn-cuedAttentionData.cueOn;
tmax= cuedAttentionData.dotsOn-refTime;
relevantTrials=find(cuedAttentionData.dotLoc==1&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==35);
[bigS,bigT]=BigSBigTloop(relevantTrials,cuedAttentionData,refTime,tmin,tmax);
[nsp] = hist(bigS,edges);
[nms,be] = hist(bigT,edges);
cueUpAll=(1000*nsp./nms)';
relevantTrials=find(cuedAttentionData.dotLoc==0&cuedAttentionData.cueOn<cuedAttentionData.dotsOn&cuedAttentionData.cellNumber==35);
[bigS,bigT]=BigSBigTloop(relevantTrials,cuedAttentionData,refTime,tmin,tmax);
[nsp] = hist(bigS,edges);
[nms,be] = hist(bigT,edges);
cueDownAll=(1000*nsp./nms)';

% Plot it
subplot(3,2,2)
hold on
plot(edges/1000,cueUpAll,'color',[.6 0 .4],'linewidth',2)
box off
set(gca,'TickDir','out');
for iUp=1:length(allSpikesUp)
    rasterUp=allSpikesUp{iUp}-cueUp(iUp);
    ys=ones(length(rasterUp),1)*iUp+100;
    plot(rasterUp,ys,'k.','markersize',1)
    hold on
end
xlim([-.1,.4])
ax=gca;
ax.FontSize=13;
ylabel('Firing rate (sp/s)             Trial','fontsize',16)
xlabel('Time from cue onset (s)','fontsize',16)
box off
line([0 0],[0 length(allSpikesUp)+100],'linewidth',2,'color','k')
set(gca,'TickDir','out');
yticks([0 50 100 150 200])
yticklabels([0 50 0 50 100])
ylim([0 275])

subplot(3,2,4)
hold on
plot(edges/1000,cueDownAll,'color',[.6 0 .4],'linewidth',2)
box off
set(gca,'TickDir','out');
for iDown=1:length(allSpikesDown)
    rasterDown=allSpikesDown{iDown}-cueDown(iDown);
    ys=ones(length(rasterDown),1)*iDown+100;
    plot(rasterDown,ys,'k.','markersize',1)
    hold on
end
xlim([-.1,.4])
ax=gca;
ax.FontSize=13;
ylabel('Firing rate (sp/s)             Trial','fontsize',16)
xlabel('Time from cue onset (s)','fontsize',16)
set(gca,'TickDir','out');
line([0 0],[0 length(allSpikesDown)+100],'linewidth',2,'color','k')
yticks([0 50 100 150 200])
yticklabels([0 50 0 50 100])
ylim([0 275])

% Now we will generate the cue aligned PSTHs of averaged activity across
% neurons (panels E and F)
% We will use a loop that makes panel E and then panel F
monkeys = [1 4]; % monkey 1 is Dm, monkey 4 is Np
for m=1:2
    monkey = monkeys(m);
    % As for the single neuron PSTH, we set the bin width
    binwidth=5;
    % And the x-axis of the figure
    edges = -100:binwidth:500;
    % set the time we are aligning the data to
    refTime = cuedAttentionData.cueOn;
    % set a minimum time to include in the plots
    tmin = cuedAttentionData.targetsOn-refTime + .2;
    % set a maximum time to include in the plots, because dot duration 
    % varies between trials and we don't want to include activity from
    % later epochs
    tmax = cuedAttentionData.dotsOn-cuedAttentionData.cueOn;
    % select trials from a single monkey
    relevantTrials=find(cuedAttentionData.monkey==monkey);
    [bigS,bigT]=BigSBigTloop(relevantTrials,cuedAttentionData,refTime,tmin,tmax);
    % Caluclate the number of spikes in each bin
    [nsp] = hist(bigS,edges);
    % Calculate the number of ms that contributed to each bin
    [nms,be] = hist(bigT,edges);
    % Calculate the spike rate from the number of spikes and the amount of
    % time that contributes to each bin
    cueAlignedAll = (1000*nsp./nms)';
    % make the figure
    subplot(3,2,4+m)
    plot(edges,cueAlignedAll,'linewidth',2,'color',[.6 0 .4])
    set(gca,'fontsize',14)
    xlabel('Time from cue onset (s)','fontsize',18)
    ylabel('Firing rate (sp/s)','fontsize',18)
    xlim([-100 400])
    ylim([0 20])
    box off
    line([0 0],[0 20],'color','k','linewidth',2,'linestyle',':')
    set(gca,'TickDir','out');
    xticks([-100 0 100 200 300 400])
    xticklabels([-.1 0 .1 .2 .3 .4])

end


%% Using the variable location spike data: figure 9A--C

% single cell plots (panels A-D)

load('variableLocationData')

% we start by setting the bin width
binwidth=10;
% and the x-axis of the figures we are going to plot
edges = -200:binwidth:500;
% select only data from example cell Ap19
cellData=variableLocationData(variableLocationData.cellNumber==19,:); % Ex cell from monkey A in the experimenter controlled duration task

% Generate the PSTHs for panels A and B
% set the reference time, which everything will be aligned to
refTime = cellData.dotsOn;
% set a minimum time to include in the plots
tmin = cellData.targetsOn-refTime;
% set a maximum amount of time to include in the plots. This is important 
% because RDM duration varies between trials and we don't want to include 
% activity from later epochs
tmax = cellData.duration +.1;
% Select trials with the motion above the fixation point
relevantTrials=find(cellData.dotLoc==1);
[bigS,bigT]=BigSBigTloop(relevantTrials,cellData,refTime,tmin,tmax);
% Caluclate the number of spikes in each bin
[nsp] = hist(bigS,edges);
% Calculate the number of ms that contributed to each bin
[nms,be] = hist(bigT,edges);
% divide the number of spikes by the number of ms and multiply by 1000 to
% get the spike rate per second
dotsUpAll = (1000*nsp./nms)';

%  repeat for trials with motion stimulus below fixation
relevantTrials=find(cellData.dotLoc==0);
[bigS,bigT]=BigSBigTloop(relevantTrials,cellData,refTime,tmin,tmax);
[nsp] = hist(bigS,edges); % number of spikes
[nms,be] = hist(bigT,edges); % number of ms.
dotsDownAll = (1000*nsp./nms)';


% Generate rasters for this cell in both conditions
% Select trials with the appropriate dot location
rasterConditionDown=cellData.dotLoc==0;
rasterConditionUp=cellData.dotLoc==1;
% Get spikes for each condition
allSpikesDown=cellData.spCell(rasterConditionDown);
allSpikesUp=cellData.spCell(rasterConditionUp);
% And the time to align to
dotsDown=cellData.dotsOn(rasterConditionDown);
dotsUp=cellData.dotsOn(rasterConditionUp);

% make the figure
figure
subplot(1,4,1)
% plot the psth
plot(edges/1000,dotsUpAll,'linewidth',2)
xlabel('Time from motion onset (s)','fontsize',12)
ylabel('Firing rate (sp/s)','fontsize',12)
hold on
% plot the rasters
for trial=1:150
    rasterUp=allSpikesUp{trial}-dotsUp(trial);
    ys=ones(length(rasterUp),1)*trial+100;
    plot(rasterUp,ys,'k.','markersize',1)
    hold on
end
xlim([-.1,.4])
ax=gca;
ax.FontSize=12;
ylabel('Firing rate (sp/s)       Trial','fontsize',14)
xlabel('Time from motion onset (s)','fontsize',14)
box off
line([0 0],[0 250],'color','k')
set(gca,'TickDir','out');
yticks([0 50 100 150 200 250])
yticklabels([0 50 0 50 100 150])

% plot the other RDM location
subplot(1,4,2)
plot(edges/1000,dotsDownAll,'linewidth',2)
hold on
for trial=1:150
    rasterDown=allSpikesDown{trial}-dotsDown(trial);
    ys=ones(length(rasterDown),1)*trial+100;
    plot(rasterDown,ys,'k.','markersize',1)
    hold on
end
xlim([-.1,.4])
ax=gca;
ax.FontSize=12;
ylabel('Firing rate (sp/s)       Trial','fontsize',14)
xlabel('Time from motion onset (s)','fontsize',14)
box off
line([0 0],[0 250],'color','k')
ylim([0 250])
set(gca,'TickDir','out');
yticks([0 50 100 150 200 250])
yticklabels([0 50 0 50 100 150])


% PSTH of mean activity for all neurons for monkeys Ap and Dz (panel C)

load('variableLocationData')
% exclude data from monkey Dm, who did not show signs of motion integration
% in LIP in this task
variableLocationData=variableLocationData(variableLocationData.monkey~=1,:);
% As for the single neuron PSTH, we set the bin width
binwidth=5;
% And the x-axis of the figure
edges = -100:binwidth:500;
% set the time we are aligning the data to
refTime = variableLocationData.dotsOn;
% set the earliest time to include for each trial
tmin = variableLocationData.targetsOn-refTime + .2;
% and the lastest time to include for each trial
tmax = variableLocationData.duration;
% exclude ODR tasks and incomplete trials
relevantTrials=find(variableLocationData.trialType==20&variableLocationData.complete==1);
[bigS,bigT]=BigSBigTloop(relevantTrials,variableLocationData,refTime,tmin,tmax);
% get the number of spikes in each bin
[nsp] = hist(bigS,edges);
% get the number of ms contributing to each bin
[nms,be] = hist(bigT,edges);
% divide the number of spikes by the number of ms and multiply by 1000 to
% get the spike rate as sp/s
motionAlignedAll = (1000*nsp./nms)';
% make the figure
subplot(1,4,3)
plot(edges./1000,motionAlignedAll,'linewidth',1.5)
xlim([-0.1 0.4])
box off
line([0 0],[0 50],'color','k')
set(gca, 'FontSize', 12);
set(gca,'TickDir','out');
xlabel('Time from motion onset (s)','fontsize',14)
ylabel('Firing rate (sp/s)','fontsize',14)

%% Table fields for cued attention task (variable location data is similar)
%
% 1 Cell number -- each neuron recorded within a task was asigned a number,
% in the order recorded. Trials with the same "cell number" are from the
% same neuron. Cells were numbered beginning at 1 for each of the two
% tasks, and no cells were recorded in both tasks.
%
% 2 Monkey -- each monkey was assigned a number. 1 -- monkey Dm; 2 --
% monkey Ap; 3 monkey -- Dz; 4 monkey -- Np. Monkey Dm participated in both
% tasks and is monkey 1 in both. However, different (adjacent) grid holes
% were used for the two tasks.
%
% 3 Complete -- aborted trials are not included in this table
%
% 4 Direction -- direction of motion of RDM patch 1
%
% 5 Direction2 -- direction of motion of RDM patch 2
%
% 6 coherence -- unsigned fraction coherent motion x1000 for RDM patch 1.
%
% 7 coherence2 -- unsigned fraction coherent motion x1000 for RDM patch 2.
% In this data set coherence == coherence2
%
% 8 correct -- 1 if monkey was rewarded, 0 if monkey was not rewarded
%
% 9 incorrect -- 0 if monkey was rewarded, 2 if monkey was not rewarded
%
% 10 duration -- duration of the RDM stimulus
%
% 11 saccadeTime -- time relative to session start of response saccade
%
% 12 rewardTime -- time relative to session start that reward delivery
% began
%
% 13 and 14 noChoice and FixBreak refer to error types not included in this
% table
%
% 15 and 16 seedBase and seedVar allow for reconstruction of the dot motion
% stimulus on each trial. Contact naomi.neva.odean@gmail.com if you need
% to do this
%
% 17 dotPos gives the position of the cued dot motion stimulus
%
% 18 dotChange is not relevant here
%
% 19 dotsOn is the time relative to session start of the onset of RDM
%
% 20 fixOn is the time relative to session start of the appearance of the
% fixation point
%
% 21 targetsOn is the time relatvie to the session start of the appearance
% of the targets
%
% 22 fixOff is the time relative to the session start of the disappearance
% of the fixation point
%
% 23 saccadeComplete is the time relative to the session start that the
% monkey's response was detected
%
% 24 and 25 targ1Pos and targ2Pos give the position of the two targets in
% visual degrees*10. The first column is the x coordinate and the second
% column is the y coordinate
%
% 26 trialType is not relevant here
%
% 27 signedCoherence motion strength and direction. Positive values are
% towards the RF, negative values are away from the RF. Divide by 1000 to
% get the fraction coherent motion.
%
% 28 dotLoc is 1 if the cued motion patch is above the fixation point and 0
% if it is below
%
% 29 dotsOff is the time relative to the session start that the RDM
% disappears
%
% 30 cueOn is the time relative to the session start that the cue appears
%
% 31 cueOff is the time relative to the session start that the cue is
% extinguished
%
% 32-34 trialsSinceSwitch, meanFR, and stdFR are not relevant here
%
% 35 right is 1 if the direction of motion was to the right and 0 if motion
% was to the left
%
% 36 spCell is a cell array containing the time relative to the session
% start of all spikes in the trial
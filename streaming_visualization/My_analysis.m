% For Data_Sorting.m testing

% clear all;
close all;
%clc;

threshold = 0.8; % spike detection threshold

Channels = 1:16;  % number of channels, here we have 16
Events = 1 : 6;   % number of events, here we have 7 events in total
Positions = 1 : 9; % number of robot positions

BinSiz = 20;   % bin widths for PSTH (s)
MinTim = 200;  % minimum time for PSTH before event
MaxTim = 600; % maximum time for PSTH after event
edges = 0 : BinSiz : MinTim + MaxTim;

SelectChn = 10;
SelectEvt = 2;
TrialType = 'All';
Directory = ('C:\Users\admin\Desktop\Data'); % Directory of all the data files
Name = 'Test_24-Sep-2017_event_0023.mat';
% Data = Data_Sorting(threshold, Name, Directory, File_Number); % read the neural data

% while(1)
Data = Data_Sorting(threshold, Name, Directory); % read the neural data
Positions = unique(Data.event_times(1,(Data.event_times(1,:)~=0 & Data.event_times(1,:)<10)));
PosSize = size(Positions);
count = 1;
% pause(10)

% Data.spike_waves
Data.spike_times
Data.event_times
% end

% EventData = [Data.event_times(1,:); Data.event_times(2:end, :)/1000];
if strcmp(TrialType, 'All')
    Ind =  ~isnan( Data.event_times(end,:) );
    EventData = Data.event_times(:,  Ind);
    SpikeData = Data.spike_times(SelectChn, Ind);
    SpikeWaveData = Data.spike_waves(Channels, Ind );
elseif strcmp(TrialType, 'Correct')
    Ind = ~isnan( Data.event_times(end,:) ) & Data.event_times(end,:)~=0;
    EventData = Data.event_times(:, Ind );
    SpikeData = Data.spike_times(SelectChn, Ind);
    SpikeWaveData = Data.spike_waves(Channels, Ind );
elseif strcmp(TrialType, 'Incorrect')
    Ind = ~isnan( Data.event_times(end,:) ) & Data.event_times(end,:)==0;
    EventData = Data.event_times(:,  Ind);
    SpikeData = Data.spike_times(SelectChn, Ind);
    SpikeWaveData = Data.spike_waves(Channels, Ind );
end

% EventData = Data.event_times;
% SpikeData = Data.spike_times(SelectChn, :);

% Sorting the data based on event time
if size(EventData, 2)>0
    PosCnt = ones(size(Positions)); % position counter
    for trl = 1 : length(SpikeData)
        
        SpkSortData{EventData(1, trl), PosCnt(EventData(1, trl))} = SpikeData{trl}(SpikeData{trl}>=EventData(SelectEvt,trl)-MinTim & SpikeData{trl}<=EventData(SelectEvt,trl)+MaxTim)-EventData(SelectEvt,trl) + MinTim; % spike time
        %         SpkSorWavetData{EventData(1, trl), PosCnt(EventData(1, trl))} = SpikeWaveData{trl}(SpikeData{trl}>=EventData(SelectEvt,trl)-MinTim & SpikeData{trl}<=EventData(SelectEvt,trl)+MaxTim, :); % spike time
        PosCnt(EventData(1, trl)) = PosCnt(EventData(1, trl)) + 1 ; % Position Ind counter
        
    end
else
    error('You do not have any recorded trials')
end
tic
%%
Ylim = 2;
for pos = 1 : 9
    
    if(PosSize>=count)
    if(Positions(count)== pos)
        subplot(3, 3, pos)
        Yscl = 0.1;
        trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
        for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))

            h1 = plot(SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 'sk', 'markersize', 2, 'markerfacecolor', 'k');
            Yscl = Yscl + 0.1;
            hold on

        end
        plot(MinTim*[1 1], [0 Ylim], ':k')
        axis([0 MaxTim+MinTim 0 Ylim])
        axis off
        set(gcf, 'color', 'w')
    count = count + 1;
    end
    end
end
toc

figure
Ylim = 3500;
for pos = 1 : 9
    
    subplot(3, 3, pos)
    psth1 = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
    ph = bar(edges(1:end-1), 1000*psth1(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
    set(ph,'edgecolor', 'k', 'facecolor', 'k');
    hold on
    plot(MinTim*[1 1], [0 Ylim], ':k')
    axis([0 MaxTim+MinTim 0 Ylim])
    %     axis off
    set(gcf, 'color', 'w')
    
end

%%
FontSiz = 6;

NumSpktoPlot = 100;
NumTrl = size(EventData, 2);
if NumTrl >= 5
    for chInd = 1 : length(Channels)
        figure(500)
        SpkWavFig = subplot(4,4,chInd);
        SWave = cell2mat(SpikeWaveData(chInd, NumTrl-4: NumTrl)');
        RndInd = randperm(size(SWave,1));
        if size(SWave,1) >= NumSpktoPlot
            plot(SpkWavFig, SWave(RndInd(1:NumSpktoPlot), 2:end)', 'k')
        else
            plot(SpkWavFig, SWave(:, 2:end)', 'k')
        end
        set(SpkWavFig, 'Xlim', [0 size(SWave,2)])
        set(SpkWavFig, 'box', 'off')
        set(SpkWavFig, 'tickdir', 'out')
        set(SpkWavFig, 'fontsize', FontSiz)
    end
else
    for chInd = 1 : length(Channels)
        figure(500)
        SpkWavFig = subplot(4,4,chInd);
        SWave = cell2mat(SpikeWaveData(chInd, :)');
        RndInd = randperm(size(SWave,1));
        if size(SWave,1) >= NumSpktoPlot
            plot(SpkWavFig, SWave(RndInd(1:NumSpktoPlot),  2:end)', 'k')
        else
            plot(SpkWavFig, SWave(:,  2:end)', 'k')
        end
        set(SpkWavFig, 'Xlim', [0 size(SWave,2)])
        set(SpkWavFig, 'box', 'off')
        set(SpkWavFig, 'tickdir', 'out')
        set(SpkWavFig, 'fontsize', FontSiz)
    end
end
set(gcf,'color','w')
set(gcf, 'Position', [0 0 480 480]) % Plos Comp supp
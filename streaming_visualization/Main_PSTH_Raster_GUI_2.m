function varargout = Main_PSTH_Raster_GUI_2(varargin)

% MAIN_PSTH_RASTER_GUI_2 MATLAB code for Main_PSTH_Raster_GUI_2.fig
%      MAIN_PSTH_RASTER_GUI_2, by itself, creates a new MAIN_PSTH_RASTER_GUI_2 or raises the existing
%      singleton*.
%
%      H = MAIN_PSTH_RASTER_GUI_2 returns the handle to a new MAIN_PSTH_RASTER_GUI_2 or the handle to
%      the existing singleton*.
%
%      MAIN_PSTH_RASTER_GUI_2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_PSTH_RASTER_GUI_2.M with the given input arguments.
%
%      MAIN_PSTH_RASTER_GUI_2('Property','Value',...) creates a new MAIN_PSTH_RASTER_GUI_2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Main_PSTH_Raster_GUI_2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Main_PSTH_Raster_GUI_2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Main_PSTH_Raster_GUI_2

% Last Modified by GUIDE v2.5 29-Jul-2017 19:31:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Main_PSTH_Raster_GUI_2_OpeningFcn, ...
    'gui_OutputFcn',  @Main_PSTH_Raster_GUI_2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Main_PSTH_Raster_GUI_2 is made visible.
function Main_PSTH_Raster_GUI_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Main_PSTH_Raster_GUI_2 (see VARARGIN)

% Choose default command line output for Main_PSTH_Raster_GUI_2
handles.output = hObject;

%Timer for refreshing the screen
%-----------------------------------

allTimers = timerfindall;

if(~isempty(allTimers))
    delete(allTimers)
end

% START USER CODE
% Create a timer object to fire at 1/10 sec intervals
% Specify function handles for its start and run callbacks
handles.timer = timer(...
    'ExecutionMode', 'fixedRate', ...       % Run timer repeatedly
    'Period', 20, ...                        % Initial period is 1 sec.
    'TimerFcn', {@update_display,hObject}); % Specify callback function

% defaul values
Tresh_Def = str2double(get(handles.Tresh,'String'));
if isnan(Tresh_Def)
    set(handles.Tresh,'String','0.8');
end
BinSiz_Def = str2double(get(handles.Bin_Width, 'String')); % Bin size for PSTH
if isnan(BinSiz_Def)
    set(handles.Bin_Width,'String','20');
end
MinTim_Def = str2double(get(handles.Pre_Event_Wind, 'String'));  % minimum time for PSTH before event
if isnan(MinTim_Def)
    set(handles.Pre_Event_Wind,'String','100');
end
MaxTim_Def = str2double(get(handles.Post_Event_Wind, 'String'));  % amximum time for PSTH before event
if isnan(MaxTim_Def)
    set(handles.Post_Event_Wind,'String','300');
end
SelectChn_Def = str2double(get(handles.Channel, 'String'));  %   Selected Channel
if isnan(SelectChn_Def)
    set(handles.Channel,'String','1');
end

SelectEvt_Def = get(handles.Event_Mark, 'Value'); % selected event, all spikes time will be aligned to this
if SelectEvt_Def>=2
    set(handles.Event_Mark,'Value',2);
end

guidata(hObject, handles);

% UIWAIT makes Main_PSTH_Raster_GUI_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Main_PSTH_Raster_GUI_2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Start_Timer.
function Start_Timer_Callback(hObject, eventdata, handles)
% hObject    handle to Start_Timer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% START USER CODE
% Only start timer if it is not running
if strcmp(get(handles.timer, 'Running'), 'off')
    start(handles.timer);
end
% END USER CODE

% --- Executes on button press in StopTimer.
function StopTimer_Callback(hObject, eventdata, handles)
% hObject    handle to StopTimer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% START USER CODE
% Only stop timer if it is running
clc
fprintf('*** You just stopped the GUI! ***\n\n')

if strcmp(get(handles.timer, 'Running'), 'on')
    stop(handles.timer);
end
% END USER CODE

% START USER CODE
function update_display(hObject,eventdata,hfigure)
% Timer timer1 callback, called each time timer iterates.
% Gets surface Z data, adds noise, and writes it back to surface object.
% clc
handles = guidata(hfigure);
% get(handles.timer, 'TasksExecuted')
if get(handles.timer, 'TasksExecuted')==1
    clc
    fprintf('*** The GUI has started working baby! ***\n\n')
else
    clc
    fprintf('*** The GUI is analysing the data honey! ***\n\n')
end
Channels = 1:4;  % number of channels, here we have 16
Events = 1 : 6;   % number of events, here we have 7 events in total
Positions = 1 : 9; % number of robot positions


Path_Def = get(handles.Data_Path,'String');
if isempty(Path_Def)
    [File_Name, File_Path] = uigetfile ({'*.*', 'All files'});
    set(handles.Data_Path,'String',[File_Path File_Name]);
else
    My_Path = get(handles.Data_Path,'String');
    File_Name = My_Path(find(My_Path=='\', 1, 'last'):end);
    File_Path = My_Path(1:find(My_Path=='\', 1, 'last'));
end

% get the updated values
BinSiz = str2double(get(handles.Bin_Width, 'String')); % Bin size for PSTH
MinTim = str2double(get(handles.Pre_Event_Wind, 'String'));  % minimum time for PSTH before event
MaxTim = str2double(get(handles.Post_Event_Wind, 'String'));  % amximum time for PSTH before event
edges = 0 : BinSiz : MinTim + MaxTim;
SelectChn = str2double(get(handles.Channel, 'String'));  %   Selected Channel
SelectEvt = get(handles.Event_Mark, 'Value'); % selected event, all spikes time will be aligned to this
threshold = str2double(get(handles.Tresh, 'String'));  %   Treshold for spike detection
SelectedTrialType = get(handles.TrialForAnalysing, 'Value');  %   Selected Trial Type
TrialTypes = get(handles.TrialForAnalysing, 'string');
TrialType = TrialTypes{SelectedTrialType};

% load the data

Data = Data_Sorting(threshold, File_Name, File_Path); % read the neural data
Positions = unique(Data.event_times(1,(Data.event_times(1,:)~=0 & Data.event_times(1,:)<10)))
[m ,PosSize] = size(Positions);
if(PosSize<9)
PosSize = PosSize - 1;
end
count = 1;

if strcmp(TrialType, 'All')
    EventData = Data.event_times;
    SpikeData = Data.spike_times(SelectChn, :);
    SpikeWaveData = Data.spike_waves(Channels, : );
elseif strcmp(TrialType, 'Correct')
    Ind = ~isnan( Data.event_times(end,:) );
    EventData = Data.event_times(:, Ind );
    SpikeData = Data.spike_times(SelectChn, Ind);
    SpikeWaveData = Data.spike_waves(Channels, Ind );
elseif strcmp(TrialType, 'Incorrect')
    Ind = isnan( Data.event_times(end,:) );
    EventData = Data.event_times(:,  Ind);
    SpikeData = Data.spike_times(SelectChn, Ind);
    SpikeWaveData = Data.spike_waves(Channels, Ind );
else
    error('Trial type is not correct') % stupid error
end

% Sorting the data based on event time
if size(EventData, 2)>0
    PosCnt = ones(size(Positions)); % position counter
    for trl = 1 : length(SpikeData)
        
        SpkSortData{EventData(1, trl), PosCnt(EventData(1, trl))} = SpikeData{trl}(SpikeData{trl}>=EventData(SelectEvt,trl)-MinTim & SpikeData{trl}<=EventData(SelectEvt,trl)+MaxTim)-EventData(SelectEvt,trl) + MinTim; % spike time
        PosCnt(EventData(1, trl)) = PosCnt(EventData(1, trl)) + 1 ; % Position Ind counter
        
    end
else
    error('You do not have any recorded trials')
end

% some axes properties
YlimRaster = [0 2];
YlimPSTH = [0 2000];
Col = [0 0 0];
Xlim = [0 MaxTim+MinTim];
FontSiz = 8;
YlabPSTH = 'Fring Rate';
XlabPSTH = 'Time (ms)';
YlabRaster = 'Trials';
XAxSpac = 100;

axes(handles.Raster1)
pos = 1 ; % Position 1
Yscl = 0.1;

% Raster
if(PosSize>=count)
if(Positions(count)==1)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster1, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster1, 'on')
end

plot(handles.Raster1, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster1, 'Ylim', YlimRaster)
set(handles.Raster1, 'Xlim', Xlim)
set(handles.Raster1, 'box', 'off')
set(handles.Raster1, 'xticklabel', [])
set(handles.Raster1, 'yticklabel', [])
set(handles.Raster1, 'tickdir', 'out')
set(handles.Raster1, 'fontsize', FontSiz)
hold(handles.Raster1, 'off')
% PSTH
axes(handles.PSTH1)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);

Ph = bar(handles.PSTH1, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH1, 'on')
plot(handles.PSTH1, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH1, 'Ylim', YlimPSTH)
set(handles.PSTH1, 'Xlim', Xlim)
set(handles.PSTH1, 'box', 'off')
set(handles.PSTH1, 'xticklabel', [])
set(handles.PSTH1, 'tickdir', 'out')
set(handles.PSTH1, 'fontsize', FontSiz)
hold(handles.PSTH1, 'off')
count = count + 1;
disp('..1')
end
end

%
axes(handles.Raster2)
pos = 2 ; % Position 2
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==2)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));

for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster2, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster2, 'on')
end

plot(handles.Raster2, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster2, 'Ylim', YlimRaster)
set(handles.Raster2, 'Xlim', Xlim)
set(handles.Raster2, 'box', 'off')
set(handles.Raster2, 'xticklabel', [])
set(handles.Raster2, 'yticklabel', [])
set(handles.Raster2, 'tickdir', 'out')
set(handles.Raster2, 'fontsize', FontSiz)
hold(handles.Raster2, 'off')
% PSTH
axes(handles.PSTH2)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH2, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH2, 'on')
plot(handles.PSTH2, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH2, 'Ylim', YlimPSTH)
set(handles.PSTH2, 'Xlim', Xlim)
set(handles.PSTH2, 'box', 'off')
set(handles.PSTH2, 'xticklabel', [])
set(handles.PSTH2, 'tickdir', 'out')
set(handles.PSTH2, 'fontsize', FontSiz)
hold(handles.PSTH2, 'off')
count = count + 1;
disp('..2')
end 
end 
%
axes(handles.Raster3)
pos = 3 ; % Position 3
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==3)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));

for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster3, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster3, 'on')
end
plot(handles.Raster3, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster3, 'Ylim', YlimRaster)
set(handles.Raster3, 'Xlim', Xlim)
set(handles.Raster3, 'box', 'off')
set(handles.Raster3, 'xticklabel', [])
set(handles.Raster3, 'yticklabel', [])
set(handles.Raster3, 'tickdir', 'out')
set(handles.Raster3, 'fontsize', FontSiz)
hold(handles.Raster3, 'off')
% PSTH
axes(handles.PSTH3)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH3, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH3, 'on')
plot(handles.PSTH3, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH3, 'Ylim', YlimPSTH)
set(handles.PSTH3, 'Xlim', Xlim)
set(handles.PSTH3, 'box', 'off')
set(handles.PSTH3, 'xticklabel', [])
set(handles.PSTH3, 'tickdir', 'out')
set(handles.PSTH3, 'fontsize', FontSiz)
hold(handles.PSTH3, 'off')
count = count + 1;
disp('..3')
end
end 
%
axes(handles.Raster4)
pos = 4 ; % Position 4
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==4)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster4, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster4, 'on')
end

plot(handles.Raster4, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster4, 'Ylim', YlimRaster)
set(handles.Raster4, 'Xlim', Xlim)
set(handles.Raster4, 'box', 'off')
set(handles.Raster4, 'xticklabel', [])
set(handles.Raster4, 'yticklabel', [])
set(handles.Raster4, 'tickdir', 'out')
set(handles.Raster4, 'fontsize', FontSiz)
hold(handles.Raster4, 'off')
% PSTH
axes(handles.PSTH4)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH4, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH4, 'on')
plot(handles.PSTH4, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH4, 'Ylim', YlimPSTH)
set(handles.PSTH4, 'Xlim', Xlim)
set(handles.PSTH4, 'box', 'off')
set(handles.PSTH4, 'xticklabel', [])
set(handles.PSTH4, 'tickdir', 'out')
set(handles.PSTH4, 'fontsize', FontSiz)
hold(handles.PSTH4, 'off')
count = count + 1;
disp('..4')
end 
end 
%
axes(handles.Raster5); 
pos = 5 ; % Position 5
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==5)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));

for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster5, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster5, 'on')
end
plot(handles.Raster5, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster5, 'Ylim', YlimRaster)
set(handles.Raster5, 'Xlim', Xlim)
set(handles.Raster5, 'box', 'off')
set(handles.Raster5, 'xticklabel', [])
set(handles.Raster5, 'yticklabel', [])
set(handles.Raster5, 'tickdir', 'out')
set(handles.Raster5, 'fontsize', FontSiz)
hold(handles.Raster5, 'off')
% PSTH
axes(handles.PSTH5)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH5, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold on
plot(handles.PSTH5, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH5, 'Ylim', YlimPSTH)
set(handles.PSTH5, 'Xlim', Xlim)
set(handles.PSTH5, 'box', 'off')
set(handles.PSTH5, 'xticklabel', [])
set(handles.PSTH5, 'tickdir', 'out')
set(handles.PSTH5, 'fontsize', FontSiz)
hold(handles.PSTH5, 'off')
count = count + 1;
disp('..5')
end
end

%
axes(handles.Raster6)
pos = 6 ; % Position 6
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==6)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster6, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster6, 'on')
end
plot(handles.Raster6, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster6, 'Ylim', YlimRaster)
set(handles.Raster6, 'Xlim', Xlim)
set(handles.Raster6, 'box', 'off')
set(handles.Raster6, 'xticklabel', [])
set(handles.Raster6, 'yticklabel', [])
set(handles.Raster6, 'tickdir', 'out')
set(handles.Raster6, 'fontsize', FontSiz)
hold(handles.Raster6, 'off')
% PSTH
axes(handles.PSTH6)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH6, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH6, 'on')
plot(handles.PSTH6, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH6, 'Ylim', YlimPSTH)
set(handles.PSTH6, 'Xlim', Xlim)
set(handles.PSTH6, 'box', 'off')
set(handles.PSTH6, 'xticklabel', [])
set(handles.PSTH6, 'tickdir', 'out')
set(handles.PSTH6, 'fontsize', FontSiz)
hold(handles.PSTH6, 'off')
count = count + 1;
disp('..6')
end
end

%
axes(handles.Raster7)
pos = 7 ; % Position 7
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==7)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster7, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster7, 'on')
end
plot(handles.Raster7, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster7, 'Ylim', YlimRaster)
set(handles.Raster7, 'Xlim', Xlim)
set(handles.Raster7, 'box', 'off')
set(handles.Raster7, 'xticklabel', [])
set(handles.Raster7, 'yticklabel', [])
set(handles.Raster7, 'tickdir', 'out')
set(handles.Raster7, 'fontsize', FontSiz)
ylabel(handles.Raster7, YlabRaster)
hold(handles.Raster7, 'off')
% PSTH
axes(handles.PSTH7)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH7, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold on
plot(handles.PSTH7, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH7, 'Ylim', YlimPSTH)
set(handles.PSTH7, 'Xlim', Xlim)
set(handles.PSTH7, 'box', 'off')
set(handles.PSTH7, 'tickdir', 'out')
set(handles.PSTH7, 'xtick', [Xlim(1) MinTim:XAxSpac:Xlim(2)])
set(handles.PSTH7, 'xticklabel', strsplit(num2str([Xlim(1) MinTim:XAxSpac:Xlim(2)]-MinTim)))
set(handles.PSTH7, 'fontsize', FontSiz)
ylabel(handles.PSTH7, YlabPSTH)
xlabel(handles.PSTH7, XlabPSTH)
hold(handles.PSTH7, 'off')
count = count + 1;
disp('..7')
end
end

%
axes(handles.Raster8)
pos = 8 ; % Position 8
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==8)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster8, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster8, 'on')
end
plot(handles.Raster8, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster8, 'Ylim', YlimRaster)
set(handles.Raster8, 'Xlim', Xlim)
set(handles.Raster8, 'box', 'off')
set(handles.Raster8, 'xticklabel', [])
set(handles.Raster8, 'yticklabel', [])
set(handles.Raster8, 'tickdir', 'out')
set(handles.Raster8, 'fontsize', FontSiz)
hold(handles.Raster8, 'off')
% PSTH
axes(handles.PSTH8)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH8, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH8, 'on')
plot(handles.PSTH8, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH8, 'Ylim', YlimPSTH)
set(handles.PSTH8, 'Xlim', Xlim)
set(handles.PSTH8, 'box', 'off')
set(handles.PSTH8, 'xticklabel', [])
set(handles.PSTH8, 'tickdir', 'out')
set(handles.PSTH8, 'fontsize', FontSiz)
hold(handles.PSTH8, 'off')
count = count + 1;
disp('..8')
end
end
%
axes(handles.Raster9)
pos = 9 ; % Position 9
Yscl = 0.1;
% Raster
if(PosSize>=count)
if(Positions(count)==9)
trlInd = find(~cellfun(@isempty,SpkSortData(pos,:)));
for trlpos = 1 : sum( ~cellfun(@isempty,SpkSortData(pos,:)))
    plot(handles.Raster9, SpkSortData{pos, trlInd(trlpos)}, Yscl*ones(size(SpkSortData{pos, trlInd(trlpos)})), 's',...
        'markersize', 2, 'markerfacecolor', Col, 'markeredgecolor', Col);
    Yscl = Yscl + 0.1;
    hold(handles.Raster9, 'on')
end
plot(handles.Raster9, MinTim*[1 1], YlimRaster, ':k')
set(handles.Raster9, 'Ylim', YlimRaster)
set(handles.Raster9, 'Xlim', Xlim)
set(handles.Raster9, 'box', 'off')
set(handles.Raster9, 'xticklabel', [])
set(handles.Raster9, 'yticklabel', [])
set(handles.Raster9, 'tickdir', 'out')
set(handles.Raster9, 'fontsize', FontSiz)
hold(handles.Raster9, 'off')
% PSTH
axes(handles.PSTH9)
psth_data = sum(cell2mat(cellfun(@(x) histc(x, edges), SpkSortData(pos, :), 'UniformOutput',false)), 2);
Ph = bar(handles.PSTH9, edges(1:end-1), 1000*psth_data(1:end-1)/(sum( ~cellfun(@isempty,SpkSortData(pos,:)))*BinSiz),'histc');
set(Ph,'edgecolor', Col, 'facecolor', Col);
hold(handles.PSTH9, 'on')
plot(handles.PSTH9, MinTim*[1 1], YlimPSTH, ':k')
set(handles.PSTH9, 'Ylim', YlimPSTH)
set(handles.PSTH9, 'Xlim', Xlim)
set(handles.PSTH9, 'box', 'off')
set(handles.PSTH9, 'xticklabel', [])
set(handles.PSTH9, 'tickdir', 'out')
set(handles.PSTH9, 'fontsize', FontSiz)
hold(handles.PSTH9, 'off')
count =1;
disp('..9')
end
end

% Plot the spike wave form
% stop(handles.timer)
handles.SpkWave = figure(500);
NumSpktoPlot = 100;
NumTrl = size(EventData, 2);
if NumTrl >= 5
    for chInd = 1 : length(Channels)
        
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
        set(SpkWavFig, 'yticklabel', [])
    end
else
    for chInd = 1 : length(Channels)
        SpkWavFig = subplot(4,4,chInd);
        SWave = cell2mat(SpikeWaveData(chInd, :)');
        RndInd = randperm(size(SWave,1));
        if size(SWave,1) >= NumSpktoPlot
            plot(SpkWavFig, SWave(RndInd(1:NumSpktoPlot), 2:end)', 'k')
        else
            plot(SpkWavFig, SWave(:, end)', 'k')
        end
        set(SpkWavFig, 'Xlim', [0 size(SWave,2)])
        set(SpkWavFig, 'box', 'off')
        set(SpkWavFig, 'tickdir', 'out')
        set(SpkWavFig, 'fontsize', FontSiz)
        set(SpkWavFig, 'yticklabel', [])
    end
end
set(handles.SpkWave,'color','w')
set(handles.SpkWave, 'Position', [0 100 700 700]) % Plos Comp supp
pause(2)
% END USER CODE


function Bin_Width_Callback(hObject, eventdata, handles)
% hObject    handle to Bin_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bin_Width as text
%        str2double(get(hObject,'String')) returns contents of Bin_Width as a double


% --- Executes during object creation, after setting all properties.
function Bin_Width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bin_Width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Pre_Event_Wind_Callback(hObject, eventdata, handles)
% hObject    handle to Pre_Event_Wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pre_Event_Wind as text
%        str2double(get(hObject,'String')) returns contents of Pre_Event_Wind as a double


% --- Executes during object creation, after setting all properties.
function Pre_Event_Wind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pre_Event_Wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Post_Event_Wind_Callback(hObject, eventdata, handles)
% hObject    handle to Post_Event_Wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Post_Event_Wind as text
%        str2double(get(hObject,'String')) returns contents of Post_Event_Wind as a double


% --- Executes during object creation, after setting all properties.
function Post_Event_Wind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Post_Event_Wind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Channel_Callback(hObject, eventdata, handles)
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Channel as text
%        str2double(get(hObject,'String')) returns contents of Channel as a double


% --- Executes during object creation, after setting all properties.
function Channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Event_Mark.
function Event_Mark_Callback(hObject, eventdata, handles)
% hObject    handle to Event_Mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Event_Mark contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Event_Mark


% --- Executes during object creation, after setting all properties.
function Event_Mark_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Event_Mark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Tresh_Callback(hObject, eventdata, handles)
% hObject    handle to Tresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tresh as text
%        str2double(get(hObject,'String')) returns contents of Tresh as a double


% --- Executes during object creation, after setting all properties.
function Tresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Data_Path_Callback(hObject, eventdata, handles)
% hObject    handle to Data_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Data_Path as text
%        str2double(get(hObject,'String')) returns contents of Data_Path as a double


% --- Executes during object creation, after setting all properties.
function Data_Path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Data_Path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in TrialForAnalysing.
function TrialForAnalysing_Callback(hObject, eventdata, handles)
% hObject    handle to TrialForAnalysing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns TrialForAnalysing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from TrialForAnalysing


% --- Executes during object creation, after setting all properties.
function TrialForAnalysing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialForAnalysing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% START USER CODE
% Necessary to provide this function to prevent timer callback
% from causing an error after GUI code stops executing.
% Before exiting, if the timer is running, stop it.
if strcmp(get(handles.timer, 'Running'), 'on')
    stop(handles.timer);
end
% Destroy timer
delete(handles.timer)
% END USER CODE

% Hint: delete(hObject) closes the figure
delete(hObject);

function [data_out] = Data_Sorting(threshold, Input_File_Name, Directory)
%% Authors' notes
% Data sorting code
% Created by Aaron Choong and Jesslyn Sutanti
% Last modified 15/8/2017
% Supervised by Yan, Kostas, and Masoud
%
% Version 1.0
%----------------------------
% Online data sorting 
%
% Controllable variables
%-----------------------------
% sp_wave_before = the number of points retrieved before the spike
% sp_wave_after = the number of points retrieved after the spike
% min_trials = the number of new trials required to be recorded before sorting the next batch
%
% Small notes
%-----------------------------
% 1) Running this code for the first time during an online streaming will be
% the slowest speed as all the new data has not been sorted yet.
% 2) The data sorting will not run until there are 3 new trials that have
% been recorded. To change the number of new trials, change the "min_trials" variable.


%% Checking for files
% Extracting date from input name
%---------------------------
% Last 4 strings are the digits of the experiment on the day
file_number = str2num(Input_File_Name(find(Input_File_Name=='_', 1, 'last')+1:end-4));
if(file_number)
    nth_test = file_number;
else 
    nth_test = 1; % Assumes file test is the first one
end 

% Finding the second last _ in the file name
if(length(find(Input_File_Name=='_'))>2)
    File_Name = Input_File_Name(1:find(Input_File_Name=='_', 2, 'last')-1);    
else 
    File_Name = Input_File_Name;
end 

% To use the input file name as the error message
if(find(Input_File_Name == '.')>0)
    Input_File_Name = [File_Name '_data_' sprintf('%.4d',nth_test) '.filetype'];
 else 
    Input_File_Name = [Input_File_Name '.filetype']; 
end 

% Allocating all the file names
%---------------------------
% For all the event times
event_name = [File_Name '_event_' sprintf('%.4d',nth_test) '.mat'];
file_event = fullfile(Directory, event_name);

% For the raw neural data
data_name = [File_Name  '_data_' sprintf('%.4d',nth_test) '.dat'];
file_data = fullfile(Directory,data_name);

% For the timestamps for the neural data
time_name = [File_Name '_timestamps_' sprintf('%.4d',nth_test) '.bin'];
file_time = fullfile(Directory, time_name);

% To write spike information into
spike_name = [File_Name '_sp_' sprintf('%.4d',nth_test) '.mat'];
file_spike = fullfile(Directory, spike_name);
fid_spike = matfile(file_spike,'Writable',true);

% Opening time first because the time.bin file is the first one to be written into 
% In Data_Stream.m which means this file will be indicative of how many
% trials has been performed (by the reset of the time-stamps)
fid_time = fopen(file_time,'r');

% Error written back to the user
if(fid_time<0)
    error('FileName:IncorrectNames',['Could not find existing .bin timestamps file with the corresponding name: ' Input_File_Name '\n' ...
    'Please check:\n1. if you have selected the correct files to analyse.\n'...
    '2. if all the files are in the same directory.\n'...
    '3. if all files name are in this configuration : ''Name''_''Date''_''data/timestamps/event''_''4 digits number''.''dat/bin/mat''\n'...
    'Example: MyName_28-Jul-2017_timestamps_0005.bin\n']);
end 

fid_event = matfile(file_event,'Writable',false);
% Error written back to the user
if(fid_event<0)
    error('FileName:IncorrectNames',['Could not find existing .mat event file with the corresponding name: ' Input_File_Name '\n' ...
    'Please check:\n1. if you have selected the correct files to analyse.\n'...
    '2. if all the files are in the same directory.\n'...
    '3. if all files name are in this configuration : ''Name''_''Date''_''data/timestamps/event''_''4 digits number''.''dat/bin/mat''\n'...
    'Example: MyName_28-Jul-2017_timestamps_0005.bin\n']);
end 

clear('fid_event'); % Using matfile to read the file results in constantly reading the updated file
fid_event = open(file_event); % Static values read for event.mat


fid_all_temp = fopen(file_data);
% Error written back to the user
if(fid_all_temp<0)
    error('FileName:IncorrectNames',['Could not find existing .dat data file with the corresponding name: ' Input_File_Name '\n' ...
    'Please check:\n1. if you have selected the correct files to analyse.\n'...
    '2. if all the files are in the same directory.\n'...
    '3. if all files name are in this configuration : ''Name''_''Date''_''data/timestamps/event''_''4 digits number''.''dat/bin/mat''\n'...
    'Example: MyName_28-Jul-2017_timestamps_0005.bin\n']);
end

% Checks to see if the sp.mat file already exists
if(exist(file_spike,'file')==2)
    Outcome = load(file_spike,'Outcome'); % Based on a previous version
    temp_total_trial = load(file_spike,'time_ptr'); % Based on a previous version
    if(isfield(temp_total_trial,'time_ptr'))
        total_prev_trial = fid_spike.total_trial;   % Provides number of trials that have already been sorted
        fprintf('Total previous trial read: %d\n',total_prev_trial);
        time_ptr = fid_spike.time_ptr;      % Reads up to the byte which the sorting was last up to
        data_ptr = fid_spike.data_ptr;      % Reads up to the byte which the sorting was last up to
    else
        total_prev_trial = 2;   % Based on a previous version
        time_ptr = 0;       % If the .sp file is new, start the sorting from scratch
        data_ptr = 0;
    end 
 else
    Outcome.Outcome = 'Not_done'; % Based on a previous version
    time_ptr = 0;   % If the .sp file does not exist, start the sorting from scratch
    data_ptr = 0;
end 

%% The data sorting component
if(fid_all_temp>=0) % If data files exist

    % To see if the file has been sorted already
    if(strcmp(Outcome.Outcome,'Not_done')) % Based on a previous version 

    % Predefining any variables that need to be defined
    %---------------------------
    yScale = 5/(2^15); % Scaling the integer by resolution (16 signed bits)
    min_trials = 3; % Minimum number of trials before sorting the data
    sp_wave_before = 20; % Number of points before the spike
    sp_wave_after = 30; % Number of points after the spike

    % Reading all the variables and files
    %---------------------------
    fid_data_info_temp = fid_event.data_buffer_info; % To retrieve the data buffer info

    % For the waveform channels
    no_of_chans = fid_data_info_temp.Number_of_Channels;   
    Time_Resolution_in_ms = fid_data_info_temp.Time_Resolution_in_ms;  % Sampling speed

    % Reading the waveform data file
    fid_data = fopen(file_data,'r');
    fseek(fid_data,data_ptr,'bof'); % Where to begin reading from
    data = fread(fid_data,inf,'int16')*yScale;  

    % Reading the timestamps data file 
    fseek(fid_time,time_ptr,'bof'); % Where to begin reading from
    time = fread(fid_time,inf,'single'); % In milliseconds
    
    % Finding the starting index for each trial
    %---------------------------
    time_scaling = Time_Resolution_in_ms; % To find the first time value
    time_index = find(time==time_scaling);  % To find where each new trial starts 
    total_trial = length(time_index)-1; % Total number of trials in the chunk of data read (last trial is still being recorded)

    % Checking if new data has come in
    %----------------------------
    if((total_trial)<min_trials)
        if(time_ptr==0)    
          fprintf('No relevant data in test #%d file\n',nth_test)
        else
           fprintf('No new data yet in #%d file\n',nth_test)
        end         
     fid_spike.Outcome = 'No_data';
    else        
        
        % For writing the new time and data pointers for the next data
        % sorting
        if(exist(file_spike,'file')==2)
            fid_spike.time_ptr = fid_spike.time_ptr + (time_index(end))*4-4; % Because time is in single format therefore 4 bytes
            fid_spike.data_ptr = fid_spike.data_ptr + (time_index(end))*2*no_of_chans-2*no_of_chans; % Because 2 bytes per point
        else
            fid_spike.time_ptr = (time_index(end))*4-4; % Because time is in single format therefore 4 bytes
            fid_spike.data_ptr = (time_index(end))*2*no_of_chans-2*no_of_chans; % Because 2 bytes per point
        end 
        
        % Setting up the sp.mat file for writing data
        if(~exist('total_prev_trial','var'))
            total_prev_trial=2;
            fid_spike.spike_time(no_of_chans,1000) = 0;
            fid_spike.spike_wave(no_of_chans,1000) = 0;
        end 

    spike_data{no_of_chans,1} = []; % Holds spike times data 
    spike_wave{no_of_chans,1} = []; % Holds the spike waveform data
    
    
    % Finding spike waveforms for each channel
    %---------------------------
    for index = 1:no_of_chans
        
 
        %  Sorting out the interleaved waveform data & finding index for spike times
        %---------------------------
        sorted_data = data(index:no_of_chans:end); % Sorting out the interleaved data (basic sorting)
        [spike_index] = spike_detection(sorted_data,threshold); % Spike detection code

        trial = 2; % Dummy variable to index the trials to be read 
                   % Starts by reading all the values between trial 1 and 2

        temp_total_trial = total_trial+1; % Indicates which time_index to read up to          
        
        % Separating the spike times into usable information for
        % the spike waveforms
        while(trial <= temp_total_trial)            

            relevant_spikes = (spike_index((time_index(trial-1)<spike_index) & (spike_index<time_index(trial)))); % Finding relevant spikes for each trial
            spike_data{index,trial-1} = time(relevant_spikes); % Finding time stamps
            spike_wave{index,trial-1} = zeros(length(relevant_spikes),sp_wave_before+sp_wave_after+2); % 2 represents time value + max spike value

            for wave_index = 1:length(relevant_spikes)           

                % Checks to see if the spike waveform exceeds the
                % dimensions of the spike data. 
                if(relevant_spikes(wave_index)>sp_wave_before && relevant_spikes(wave_index)<time_index(trial)-sp_wave_after)
                    % Retrieves the raw neural data for the corresponding
                    % spike waveform
                    spike_wave{index,trial-1}(wave_index,:) = [time(relevant_spikes(wave_index))...
                        sorted_data(relevant_spikes(wave_index)-sp_wave_before:relevant_spikes(wave_index)+sp_wave_after)'];      
                end 
                
            end 

            trial = trial + 1;
            clear('relevant_spikes');
        end 

    end

    % Writing the sorted out data to a .mat file
    %---------------------------
    if(total_prev_trial~=2)
        fid_spike.spike_time = [fid_spike.spike_time spike_data];  % Appending new data
        fid_spike.spike_wave = [fid_spike.spike_wave spike_wave];
        fid_spike.total_trial = fid_spike.total_trial + total_trial - 1; % Last trial is still being recorded
     elseif(total_prev_trial == 2)
        fid_spike.spike_time = spike_data(:,2:end); % Writing the first set of data
        fid_spike.spike_wave = spike_wave(:,2:end);
        fid_spike.total_trial = total_trial - 1; % Last trial is still being recorded
        
        % Transferring information from event.mat to sp.mat
        fid_spike.wave_info = {'Rows are different spikes.', '1st column is the time of the spike', 'Rest are the spike waveform'};
        fid_spike.event_info = fid_event.event_buffer_info;
        fid_spike.data_info = fid_event.data_buffer_info;
        fid_spike.event_markers = {'1 = Robot Position', '2 = FootSwitch On', '3 = Green LED on', '4 = Red LED on', '5 = TouchSensor On', '6 = TouchSensor Off', '7 = First Reward'};
    end 
    
    fprintf('Number of trials sorted: %d\n',fid_spike.total_trial);
    fid_spike.raw_event_time = fid_event.event_time; % Raw event time (unsorted)
    fid_spike.Outcome = 'Data_sorted'; % Based on previous version
    fprintf('Data sorting successful for test #%d file\n',nth_test)

    clear('spike_data');
    clear('spike_wave');

    end % End of data sorting code code

    end % Nested IF to check if the file has already been sorted (based on previous version)
end % End of IF statement to see if the file exists


%% Writing spike data as an output

% Based on previous version
if(file_number)
    nth_test = file_number;
else 
    nth_test = 1;
end 

fid_spike.Outcome = 'Not_done'; % Based on previous version

event_name = [File_Name '_event_' sprintf('%.4d',nth_test) '.mat']; % Based on previous version
file_event = fullfile(Directory,event_name); % Based on previous version
event_info = open(file_event); % Based on previous version

% should we modify the event sorter then? since it will sort different
% errors maybe


event_time = event_sorter(fid_spike.raw_event_time); % Sorting the event data (correcting where the events are wrongly placed)
temp_spike_time = load(file_spike,'spike_time'); % Based on previous version

if(isfield(temp_spike_time,'spike_time'))
    data_out.event_times = event_time(:,1:end-1); % Last event is still being recorded 
    data_out.spike_waves = fid_spike.spike_wave;
    data_out.spike_times = fid_spike.spike_time;
    data_out.event_markers = event_info.event_markers;
    data_out.data_info = event_info.data_buffer_info;
    data_out.event_info = event_info.event_buffer_info;
    
else 
    data_out.spike_waves = NaN;
    data_out.spike_times = NaN;
    data_out.event_times = NaN;
    data_out.event_markers = NaN;
    data_out.data_info = NaN;
    data_out.event_info = NaN;
end

end 
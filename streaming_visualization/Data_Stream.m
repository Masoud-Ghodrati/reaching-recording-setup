% This code will stream data off the 1401 in real time 
%======================================================
% BY Aaron Choong and Jesslyn Sutanti
% Data Modified: 15/8/2017
% Supervised by: Yan, Masoud, Kostas
%======================================================
%
% Version 1.1 
%--------------------------------------------------------
% Data_Stream is for streaming experimental data from the Power1401 in real time.
% No input or output arguments required.
% To start the Data_Stream, type Data_Stream into the MATLAB command
% window.
% More instructions can be found within the code itself.

%=====================================================
% version 1.1 changes:
% using if statement instead of while statement to retrieve the channel buffer
% change buffer reset and clock at footswitch only
% fixed major bug in the channel streaming (fixed the retrieved buffer address in the channel streaming if statement)


%%
function Data_Stream
%% Initialises connection with 1401
%-----------------------------------
initialise_1401();  

%% Initilisation of variables to setup transfer buffers in the 1401 for transfer to host 
%-----------------------------------
% Sets up data buffer variables (for channel inputs)
% To change number of channel inputs, go into the data_buffer_variables...
% and change channel_array to specify which channels is required
[channel_array, single_data_length, no_of_chans, ...
    all_data_length, buffer_size, rpt, pre, cnt, scaling_factor] = data_buffer_variables();

% Sets up event buffer variables (for event inputs)
% To change the number of digital inputs, change "no_of_dig_inputs" ...
% It is currently set to 8
% 0->2 are robot positions
% 3->7 are events (reward, touchsensor, red LED, green LED and footswitch
% in that respective order)
[Units, Time_Length, event_buffer_size, no_of_dig_inputs] = event_buffer_variables();

%% Setting up event and data buffers
%-----------------------------------

% Communicates with the 1401 via matced64c - sets up event buffer
%-----------------------------------
events_buffer_setup(event_buffer_size, buffer_size);


%% Setting up the files to be saved
%-----------------------------------
% Sets up the file names and file IDs. 
% Can do more than one set of experimental trials in one day.

Name = 'Test';  % File name that is written to the start of the file
Directory = ('C:\Users\admin\Desktop\Data\'); % Directory to save all the data
[fid_data, fid_time, fid_event] = file_setup(Name, Directory); % Sets up files


%% Setting up data streaming variables
%-----------------------------------
robot_buffer_size = event_buffer_size;
no_of_data_buffers = 3; % Number of buffers space for waveform data reading and retrieving

% Empty arrays (preallocation)
%-----------------------------------
event_data = nan(event_buffer_size*(no_of_dig_inputs+3)/2,1); % Event data array (transferred from the 1401)
empty_event = nan(no_of_dig_inputs * event_buffer_size + (3*robot_buffer_size),1);% Empty event data to send to the buffer

% Setting up the time arrays to write the timestamps
%-----------------------------------
time_scaling = 4/(pre*cnt*no_of_chans); % Time resolution
time_start=(1:single_data_length)*time_scaling;  
time=time_start; % Time values
timeout = 30; % In seconds

%% Writing event markers 
%-----------------------------------
event_markers(channel_array, no_of_chans, no_of_data_buffers, buffer_size, single_data_length...
    ,time_scaling, pre , cnt ,Units, Time_Length, event_buffer_size, timeout, fid_event);

% Clears the buffer in the 1401 before commencing main stream
%-----------------------------------
matced64c('cedTo1401',buffer_size,0,zeros(buffer_size,1),1);
matced64c('cedTo1401',no_of_dig_inputs*event_buffer_size,buffer_size+2,empty_event);

%% Final setup for the data streaming
%-----------------------------------

% Variables used in the data streamming
%-----------------------------------
toggle_foot = 0; % Toggles for the rising and falling edges for the footswitch
toggle_touchsensor = 0;
toggle_red = 0;
toggle_green = 0;
time_toggle = 0;
toggle_robot_position = 0;

green_led_tracker = 0; % Tracks the current number of green LEDs on and off during one trial
max_green_led = event_buffer_size / (8*2); % Max number of green LEDs on and off for one trial
current_trial = 0; % For writing 
fid_event.event_time(no_of_dig_inputs+1,1) = 0; % For event time writing
dummy_reward = 0;
buffer_point = 1;

prev_max_event(no_of_dig_inputs + 3) = 0; % all digital inputs
event_time(no_of_dig_inputs) = 0;
event_thres = 10;
pre_ptr = floor(200/scaling_factor);
prev_robot_pos = 1;


%% Main data streaming
%-----------------------------------
% Waits for user to start the streaming
fprintf('Press any key to start to stream the data\n');
pause();

% Starts the data streaming
%-----------------------------------
fprintf('Main data streaming for the channels has started\n')
matced64c('cedSendString',['ADCMEM,I,2,0,',num2str(buffer_size),',',num2str(channel_array,'%1d '),',',num2str(rpt),',H,',num2str(pre),',',num2str(cnt),';']);


% Starts the event clock
%-----------------------------------
fprintf('\nStarted the clock for the event buffer\n')
matced64c('cedSendString',['AUDATM,G,0,',num2str(Time_Length),',',Units,',-;']);

tic % For restarting the clock if no footswitch has been pressed for x seconds

%% Actual streaming of the code
%-----------------------------------
fprintf('Ready to start experiment\n');
prev_max_event(1) = 1;
while(1)


    matced64c('cedSendString','ADCMEM,P;');
    res = eval(matced64c('cedGetString'));
    
    % Streaming first buffer
    %-----------------------------------
    if(mod(buffer_point,3)==1)
        % Checks if the bottom third of the buffer has filled. 
        % The -ve component accounts for the delay while finding the pointer,
        % hence it can start transferring at minus "x" amount 
        if res>(buffer_size/3-pre_ptr)  % -200 as there is delay when the pointer has been retrieved
        buffer_point = buffer_point + 1;
        fwrite(fid_data,matced64c('cedToHost',all_data_length,0),'int16');
        time_toggle = time_toggle + 1;
        end
    end 

    % Streaming second buffer
    %-----------------------------------
    if(mod(buffer_point,3)==2)
        if ((res>(2*buffer_size/3-pre_ptr)))
        buffer_point = buffer_point + 1;
        fwrite(fid_data,matced64c('cedToHost',all_data_length,buffer_size/3),'int16');
        time_toggle = time_toggle + 1;
        elseif(res<(buffer_size/3-pre_ptr))
        buffer_point = 1;
        fwrite(fid_data,matced64c('cedToHost',all_data_length*2,buffer_size/3),'int16');
        time_toggle = time_toggle + 1;    
        fwrite(fid_time,time,'single');
        time = time + (time_scaling * single_data_length);
        
        end
    end 

    % Streaming third buffer
    %-----------------------------------
    if(mod(buffer_point,3)==0)
        if res<(2*buffer_size/3-pre_ptr)
        buffer_point = buffer_point + 1;
        fwrite(fid_data,matced64c('cedToHost',all_data_length,buffer_size*2/3),'int16');
        time_toggle = time_toggle + 1;
        end
    end 


      % Retrieves the event data
      event_data(1:end,1)=matced64c('cedToHost',event_buffer_size*(no_of_dig_inputs+3)/2 ,buffer_size+2);
            
      % Below is the processing of the event data for each event
      % Processing of the event data for footswitch (indicated by the first
      % chunk of event_data
      
      for event_idx = 1:no_of_dig_inputs+3
          temp = unique(event_data(1+event_buffer_size/2*(event_idx-1):event_buffer_size/2*(event_idx)));   
          event_time(event_idx) = temp(end);
      end 
    
          
    % Reading the robot position BUFFER
    %-----------------------------------                           

    % Checks to see if a value was recorded for the event buffer
    if(event_time(6)>prev_max_event(6)+event_thres)
        if(toggle_robot_position)
              prev_max_event(6) = event_time(6);
              robot_counter = robot_counter + 1;
         end % End of robot position code 
         toggle_robot_position~=toggle_robot_position;
    end % End of max statement
    

    % For all other digital inputs. Mostly the same as the above (footswitch)
    %-----------------------------------   
    for index = 2:no_of_dig_inputs


                if(event_time(2)>prev_max_event(2)+event_thres)
                    prev_max_event(2) = event_time(2);
                    prev_max_event(1) = 0; 
                    matced64c('cedTo1401',event_buffer_size/2,buffer_size+2,empty_event(1:event_buffer_size/2));
                    % If GREEN LED is on, it means the footswitch is held, hence
                    % toggle is set to 1. green_LED_tracker also increments
                    %-----------------------------------
                    green_led_tracker = green_led_tracker + 1; % Tracks max number of green LEDs
                    toggle_green = ~toggle_green;
                    toggle_red=0;
                    toggle_touchsensor=0;
                    if(toggle_green)
                        fid_event.event_time(3,current_trial) = event_time(2); % Rising edge
                        fprintf('Green LED is on\n');

                        dummy_reward = 0;
                        if(robot_save)
                            
                            
                            if(robot_counter~=0)
                            fprintf('Saved robot condition value :%d\n',robot_counter);
                            fid_event.event_time(1,current_trial) = robot_counter; % Writing the robot counter
                            prev_robot_pos = robot_counter;
                            else
                            fprintf('Saved robot condition value :%d\n',prev_robot_pos);
                            fid_event.event_time(1,current_trial) = prev_robot_pos;
                            disp('in here')
                            end 
                            robot_counter = 0;
                            robot_save=0;
                            toggle_robot_position = 1;
                        end 
                    end  
                end

            % Writes to .mat file for the event times (red LED on)
            %-----------------------------------            
          if(event_time(3)>prev_max_event(3)+event_thres)
                prev_max_event(3) = event_time(3);
                toggle_green = 0;
                toggle_red = ~toggle_red;
                    if(toggle_red)
                        fid_event.event_time(4,current_trial) = event_time(3); % Rising edge
                        fprintf('Red LED is on\n');                     
                    end
          end 


                % Writes to .mat file for the event times (touchsensor on and off)
                %-----------------------------------
                if(event_time(4)>prev_max_event(4)+event_thres)
                    prev_max_event(4) = event_time(4);
                    prev_max_event(5) = 0;
                   % Writes to .mat file for the event times (touchsensor on)
                   %-----------------------------------
                    if(toggle_touchsensor==0)
                        fprintf('Touchsensor is on\n');
                        fid_event.event_time(5,current_trial) = event_time(4); % Rising edge %Touch sensor

                    % Writes to .mat file for the event times (touchsensor off)
                    %-----------------------------------
                    else
                        fid_event.event_time(6,current_trial) = event_time(4); % Rising edge %Touch sensor
                        fprintf('Touchsensor is off\n');
                    end 
 
                    toggle_touchsensor = ~ toggle_touchsensor;
                end 


                % Writes to .mat file for the event times (1st reward)
                %-----------------------------------
                if(event_time(5)>prev_max_event(5)+event_thres)
                    prev_max_event(5) = event_time(5);
                    
                    % To ensure that the reward from previous trial is not
                    % recorded instead based on a threshold value
                    if(event_time(5)>1000 &&  event_time(5)>dummy_reward+500) 
                        fid_event.event_time(7,current_trial) = event_time(5); % Rising edge
                        dummy_reward = event_time(5);
                        fprintf('Reward is given\n');
                    end 
                end             

    end % End of checking the event inputs

    
    % Checks if there's any new value - both edges (for footswitch)
    %-----------------------------------    
    if(event_time(1)>prev_max_event(1)+event_thres)
        
        if(prev_max_event(1) ~= 0)
          toggle_foot = ~toggle_foot;
        end 
        
        prev_max_event(1) = event_time(1);     
    end 
  
    for i = 1:time_toggle
    % Adds more time to the current time stamp
    % Write to time file
    %-----------------------------------
    fwrite(fid_time,time,'single');
    time = time + (time_scaling * single_data_length);
    time_toggle = 0;
    end 
    
    
   % Checks if the event clock should be reset based on the rising edge of
   % the foot switch 
   % This resets all the event buffers for each event
   %-----------------------------------  
   if((toggle_foot) || toc>timeout || green_led_tracker>max_green_led)        
        fprintf('-------------------\n');
        fprintf('Trial has reset\n'); 
        fprintf('-------------------\n');
        current_trial = current_trial + 1;  % Next pointer value to be written
        fid_event.event_time(2,current_trial) = 1; % Rising edge reset aligned to the footswitch

        % Resets each event channel buffer
        %-----------------------------------   
         for dig_index = 7:-1:0

            % Sets up each digital event input's buffer
            %-----------------------------------   
            matced64c('cedSendString',['AUDATM,A' num2str(dig_index) ','...
            num2str(buffer_size+2+event_buffer_size*(7-dig_index)) ',' num2str(event_buffer_size) ';']);

        end   % End of event channel resets

        % Resets the event clock for all events
        %-----------------------------------   
        matced64c('cedSendString',['AUDATM,G,0,',num2str(Time_Length),',',Units,';']);

        % Empties out the event buffer 
        matced64c('cedTo1401',(no_of_dig_inputs+3)*event_buffer_size,buffer_size+2,empty_event);

        matced64c('cedSendString',['ADCMEM,I,2,0,',num2str(buffer_size),',',num2str(channel_array,'%1d '),',',num2str(rpt),',H,',num2str(pre),',',num2str(cnt),';']);
        buffer_point = 1;

        tic;    % Resets the timeout 
        green_led_tracker=0;    % Resets the green LED timeout
        time = time_start;  % Resets the timestamp
        
        toggle_green = 0;
        toggle_red = 0;
        toggle_touchsensor = 0;
        toggle_robot_position = 1;
        robot_save = 1;
        robot_counter = 0;
        prev_max_event(2:end) = 0;
        toggle_foot = 0;
   end


       

end % End of data streaming code

end   % End of main function

%%----------------------------------
% Rest of the function declarations
%-----------------------------------

% Sets up connection with the Power1401
function initialise_1401()
%% Initialising variables
U14TYPEPOWER=3;  % Power1401  


%% Checking connection
res=matced64c('cedOpenX',0);  % v3+ can specify 1401
if (res < 0)
   disp(['1401 not opened, error number ' int2str(res)]);
else %opened OK  
    
%% Retrieving typeOf1401 details
typeOf1401=matced64c('cedTypeOf1401');  

%% Resetting the settings/buffer of 1401
res=matced64c('cedResetX');
res=0;

if (res > 0)
    disp('error with call to cedWorkingSet - try commenting it out');
    return
end


matced64c('cedLdX','C:\1401\','ADCMEM','AUDATM','DIG');   % To specify which commands we want to use : ADCMEM & AUDATM & DIG
end 
end 


% Sets up the event buffers
function events_buffer_setup(event_buffer_size, buffer_size)

% To set the number of digital inputs
no_of_dig_inputs = 7;   % 4 + 1 (7 -> 3 has 5 dig inputs)


    % This starts at digital input 15 (which is 7) due to the footswitch being
    % connected to digital input 7
    for dig_index = 7:-1:(7-no_of_dig_inputs)

      % Sets up each digital event input's buffer
      matced64c('cedSendString',['AUDATM,A' num2str(dig_index) ','...
      num2str(buffer_size+2+event_buffer_size*(7-dig_index)) ',' num2str(event_buffer_size) ';']);

    end   % End of for loop

end   % End of function

% Sets up the files 
function [fid_data, fid_time, fid_event] = file_setup(Name, Directory)
%% Setting up the files 
%-----------------------------------
% Retrieves the current data
date_time = date; 

% Sets up file names
nth_test = 0;
fid_temp = 1;

% To find the current experimental trial today
while(fid_temp>=0)
    nth_test = nth_test + 1;
    data_name = [Name '_' date_time  '_data_' sprintf('%.4d',nth_test) '.dat'];
    file_data = fullfile(Directory,data_name);
    fid_temp = fopen(file_data);
end 

event_name = [Name '_' date_time '_event_' sprintf('%.4d',nth_test)  '.mat'];
timestamps_name = [Name '_' date_time '_timestamps_' sprintf('%.4d',nth_test)  '.bin'];

% Retrieving the folder IDs
file_data = fullfile(Directory,data_name);
fid_data = fopen(file_data,'a');
file_time = fullfile(Directory,timestamps_name);
fid_time = fopen(file_time,'a');

file_event = fullfile(Directory,event_name);
fid_event = matfile(file_event,'Writable',true);

end

% Saving all the information of the event streaming setup
function event_markers(channel_array, no_of_chans, no_of_data_buffers, buffer_size, single_data_length ...
    ,time_scaling, pre, cnt, Units, Time_Length, event_buffer_size, timeout, fid_event)

    % Information being stored in structs
    data_info = struct('Channel_Array',channel_array,'Number_of_Channels',no_of_chans,'Number_of_Data_Buffers',no_of_data_buffers...
        ,'Data_buffer_size_in_bytes',buffer_size,'Single_Data_Length', single_data_length...
        ,'Sampling_Speed_per_Channel_in_Hz',4*10^6/(no_of_chans*pre*cnt)...
        ,'Time_Resolution_in_ms',time_scaling,'Time_out_for_trial_in_seconds',timeout);
    
    event_info = struct('Sampling_Unit',{Units},'Sampling_Rate_per_Unit',Time_Length, 'Event_Buffer_Size_in_bytes', event_buffer_size);
    
    event_markers = {'1 = Robot Position', '2 = FootSwitch On', '3 = Green LED on', '4 = Red LED on' ...
        ,'5 = TouchSensor On', '6 = TouchSensor Off', '7 = First Reward'};

    % Recording into the .mat file
    fid_event.data_buffer_info = data_info;
    fid_event.event_buffer_info = event_info;
    fid_event.event_markers = event_markers;

end 

% The place to setup all the data buffer variables
function [channel_array, single_data_length, no_of_chans, all_data_length, buffer_size, rpt, pre, cnt, scaling_factor] = data_buffer_variables()
%% Important comments which indicate how to control the data transfer (1401 -> MATLAB)
% Details on clock
%==========================================================================
% collect 100 pts at 5000 Hz from both channels 0 and 1 = 400 bytes,
% and implies a overall sampling rate of 10000 Hz for ADCMEM = 16x25
% but just set desired sample rate for ADCBST (5000Hz) = 32x25 (H clock)

%% Description for ADCMEM
% cedSendString is used to set up the buffer size and sample rates between
% HOST and 1401. ADCMEM is an inbuilt 1401 command to sample for data channels. 
% Refer to the Programming family manual (PROGMANW)

% Function input: ADCMEM,kind,byte,st,sz,chan,rpt,clock,pre,cnt Clock set up
% matced64c('cedSendString','ADCMEM,I,2,0,800,0 2 4 9,1,H,16,1000;');
% kind = I (multitasking) OR F (sequential
% byte = 2 
% "H" starts the clock immediately to start (4 MHz clock). C (1MHZ) or T (10MHZ)
% st = start point of the user space in the 1401 for the array (start
% address of the byte - must be a multiple of the argument [ 1 or 2]) 
% sz = number of bytes in the array that is stored
% rpt = the number of times to cycle around the array (buffer)
% pre & cnt works together to create a clock division (pre * cnt = clock
% division). I.E: if we want 1000 KHz for sampling rate off a 1Mhz cycle, 
% we can get pre = 2 and cnt = 500 (therefore 1 MHz / (pre*cnt) = 1 kHz)

%% Description for AUDATM
% MATLAB has to use cedSendString to communicate to the 1401 (and it is
% used to set up the commands within 1401). AUDATM is an inbuilt 1401
% command to sample the event input channels. Refer to PROGMANW manual

% Function input: AUDATM,An,st,sz
% matced64c('cedSendString','AUDATM,A7,0,100;');
% "An" = the digital input event channel that will be accessed (n = channel)
% Please note, this does not correspond to "E0 and E1" at the front of the
% 1401. It refers to the digital inputs at the back of the 1401! 
% st = the start address of the buffer
% sz = the size of the buffer

% Please note for AUDATM, it does not record any event inputs if no events
% have occurred and the buffer does not wipe any previous events that were
% already written into the buffer. The streaming code must get rid of any
% repeated values. 

%% Instructions:
%==========================================================================
% These are values for setting up the buffer transfer from the 1401 to host
% The only thing that should change is "pre" and "cnt" to get sampling rate

% Can also change the no_of_dig_inputs for the event inputs streaming 


%% Define channels and data lengths
%-----------------------------------

% Data buffer variables
%-----------------------------------
channel_array = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; % Defines each individual channels
% channel_array = [0 1 2 3]; % Defines each individual channels

no_of_chans = length(channel_array); % Defines the total number of channels that are being analysed
scaling_factor = floor(16 / no_of_chans);   % Scaling factor is to make sure the channels are sampled at ~31.25khz and the buffer size 

single_data_length=1000; % For each channel, how many sample points (cannot be too small otherwise buffer will fill up too quick)
all_data_length=single_data_length*no_of_chans; % Converts from # of sample points for a single channel to 16 channels



bytes = 2;
buffer_size = all_data_length*bytes*3;  % The *2 is for 2 buffers that can store one complete loop of data...
% ... Each buffer point only stores 1 byte so it must assign 2 bytes for each sample

rpt = 0;    % Maximum number of repeats
pre = 2;    % Divisor for the clock
cnt = 4*scaling_factor;   % Divisor for the clock

end 

% Setting up the event buffer variables
function [Units, Time_Length, event_buffer_size, no_of_dig_inputs] = event_buffer_variables()
%% Defines event buffer variables
%-----------------------------------

% Events testing 
Units = 'U'; % Sampling unit: U = microseconds, M = milliseconds
Time_Length = 1000; % Sampling rate
event_buffer_size = 400; 
no_of_dig_inputs = 5; % The main events such as footswitch->reward


end 

% Closing connection with the 1401
function  res = close_1401()
%% To close the connection with the 1401
%-----------------------------------
res=matced64c('cedCloseX'); % Disconnects host from 1401
fprintf('\n\nConnection with the 1401 has closed\n');
end 

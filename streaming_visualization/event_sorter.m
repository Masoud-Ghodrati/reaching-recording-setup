function [event_time] = event_sorter(event_time)
%%
% Event correction code 
% Correcting the robot position in sequential order
% and the touch sensor on/off if touch sensor on was recorded as touch sensor off
%
% Created by Aaron Choong and Jesslyn Sutanti
% Last modified 10/8/2017
% Supervised by Yan, Kostas, and Masoud
%
% Version 1.0
%------------------------------------

% Checks first position when it is equal to 1
index = 1;
while(event_time(1,index)~=1 && length(event_time)<index)
   index = index + 1;
end 

% Backtracks from the first value of the robot position 1 to fix any incorrect
% robot positions
if(index>2)
    for temp_index = index:-1:2
        if(event_time(7,temp_index-1)~=0 && event_time(1,temp_index-1)~=event_time(1,temp_index))
            if(event_time(1,temp_index)==1)
              event_time(1,temp_index-1) = 9;
            else 
              event_time(1,temp_index-1) = event_time(1,temp_index)-1;
            end
        end
    end 
end 

% For each event recorded (last one is still being recorded)
for index = 1:length(event_time(1,:))-1

    % If the current one is successful, next one should be plus one
    if(event_time(7,index)~=0 && event_time(1,index+1)~=event_time(1,index)+1)
        if(event_time(1,index)<9)
            event_time(1,index+1) = event_time(1,index)+1;
        else 
           event_time(1,index+1) = 1;
        end
    end 

     % If the current one is unsuccessful, the next one should equal
     % the previous one
    if(event_time(7,index)==0 && event_time(1,index+1)~=event_time(1,index))
        if(event_time(1,index)>0)
           event_time(1,index+1) = event_time(1,index);
        end
    end 

    % If touchsensor off became touchsensor on (rewritten, therefore
    % write back)
    if(event_time(7,index)~=0 && event_time(6,index)==0)
        if(event_time(7,index) - event_time(5,index)<300) % Threshold to change touchsensor on to off
            event_time(6,index) = event_time(5,index);
            event_time(5,index) = 0;              
        end 
    end 

end 

    event_time(1,event_time(1,:)<1) = 1; % Sets any unsorted robot positions as 1 (contingency)
    event_time(event_time==0) = NaN; % Returns any zeroes as NaNs
    
end % End of event sorting
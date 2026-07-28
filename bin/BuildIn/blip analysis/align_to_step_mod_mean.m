function [fwd_mod_step_data, bwd_mod_step_data] = align_to_step_mod_mean(data, min_step_length)
% JS 24/06/19 to take all values around the mean and see if you see dips
% data - steptrace.data
% time_mean_modded (array of histogram time borders) i.e. 0:0.0005:0.05;
% min_step_length (total time length of step needed to be considered) i.e. 0.05;

if nargin < 2
    min_step_length = 0.05;
end

% create meaned out data trace
mean_mod = data.trace(:,1) - data.trace(:,3);

% now go through look for all steps
idx = find(data.trace(:,5)>0);
idx = [1;idx;size(data.trace(:,5),1)]; % pad with ends

fwd_mod_step_data = cell(0,2);
bwd_mod_step_data = cell(0,2);

for i = 2:length(idx)-1
    
    %forward step longer than min_step_length
    if data.trace(idx(i),3) - data.trace(idx(i)-1,3) > 0 && data.time(idx(i)) - data.time(idx(i-1)) > min_step_length
        % get times within range
        timelog = and(data.time < data.time(idx(i)), data.time >= (data.time(idx(i)) - min_step_length));
        tidx = find(timelog>0);
        fwd_mod_step_data{size(fwd_mod_step_data,1)+1,1} = data.time(tidx) - data.time(idx(i));
        fwd_mod_step_data{size(fwd_mod_step_data,1),2} = mean_mod(tidx);
    %backward step longer than min_step_length
    elseif data.trace(idx(i),3) - data.trace(idx(i)-1,3) < 0 && data.time(idx(i+1)) - data.time(idx(i)) > min_step_length
        % get times within range
        timelog = and(data.time >= data.time(idx(i)), data.time < (data.time(idx(i)) + min_step_length));
        tidx = find(timelog>0);
        bwd_mod_step_data{size(bwd_mod_step_data,1)+1,1} = data.time(tidx) - data.time(idx(i));
        bwd_mod_step_data{size(bwd_mod_step_data,1),2} = mean_mod(tidx);
    end
    

end



end
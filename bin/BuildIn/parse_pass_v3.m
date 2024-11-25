function [startpoint,endpoint] = parse_pass_v3(xytrace)

% Based on parse_and_review_v2 from QD_analysis

%then, if prompted, it plots each trace, and prompts the user to review it,
%accept it, or reject it.
%Specifically, it asks first if each trace shall be reviewed, and then
%asks to specify a new beginning and end to the trace when prompted

% make sure we write this on a new figure
h =  findobj('type','figure');
n = length(h);
figure(n+1)

title("Select [endpoint,startpoint] OR just ENTER for entire")
scatter(1:length(xytrace), xytrace(:,1)', 'filled');
hold on
scatter(1:length(xytrace), xytrace(:,2)', 'filled');
set(gca,'XLim',[-0.05*length(xytrace),1.05*length(xytrace)]) %give us some room to crop
%prompt = ['Trace ' num2str(i) ': keep (k), modify(m), separate(s), or delete(d) [k]'];
%review = input(prompt,'s');
endpoint = length(xytrace(:,1));
startpoint = 1;
x = ginput(2);

if ~isempty(x)
    endpoint = min(round(x(1)),endpoint);
end
if length(x) == 2
    startpoint = max(round(x(2)),startpoint);
end

% JS Edit 2023/01/24 round to nearest non NaN value to make it happy for replotting
idx = find(~isnan(xytrace(:,1)));
if ~ismember(endpoint,idx)
    [~,a] = min(abs(idx - endpoint));
    endpoint = idx(a);
end
if ~ismember(startpoint,idx)
    [~,a] = min(abs(idx - startpoint));
    startpoint = idx(a);
end

% if ( strcmp(review,'m') )       %modify by clicking on end and beginning of the trace
%     disp(['select endpoint, startpoint [' num2str(endpoint) ',' num2str(startpoint) ']']);
%     %x = ginput(2);          %we assume more likely to keep existing startpoint, so default to that if one point is picked
%     if (length(x)>=1)
%         endpoint = min(round(x(1)),endpoint);
%         if (length(x)>=2)
%             startpoint = round(x(2));
%         end
%     end
%     disp(['setting end, start to ' num2str(endpoint) ',' num2str(startpoint) ]);
% else
%     disp('keeping trace');
% end        
close(figure(n+1))

end


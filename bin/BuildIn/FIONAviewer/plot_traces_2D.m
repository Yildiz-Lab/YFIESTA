function plot_traces_2D(varargin)

%This function plots traces objects in an Ahmet-friendly manner
%the first plot is green and the second is blue
%it makes a 16 nm grid, and  uses dots for the trace (like Excel)


color(1) = 'g';            %set an array of different colors
color(2) = 'b';
color(3) = 'r';
color(4) = 'k';

usage = varargin{nargin};

if (nargin < 2)
    disp('please input as traces, usage, usage is last');
end

if (nargin >= 5) 
   color (5:nargin) = 'c';      %all plots are cyan past 5 (rarely used)
end



trace = varargin{1};
min_x = min(trace(:,1));      %initialize minimum and maximum trace sizes
max_x = max(trace(:,1));
min_y = min(trace(:,2));
max_y = max(trace(:,2));

%first plot the x axis traces

figure; hold;              %make a new figure

for k = 1:nargin-1                %loop over input traces
    clear x_cut;
    trace = varargin{k};        
    plot (trace(:,1),strcat(color(k),'-'));    %plot trace with dots and dashed lines 
    x_cut(1:length(trace(:,3))) = NaN;
    for i=1:length(trace(:,3))
       if (trace(i,6) >= usage) 
           x_cut(i) = trace(i,3);
       end
    end
    plot (x_cut,color(k));                %plot fit with lines
    if (min(trace(:,1)) <= min_x )
       min_x = min(trace(:,1));           %find min and max of all traces
    end
    if (max(trace(:,1)) >= max_x )
       max_x = max(trace(:,1));
    end
end

set (gca,'YTickMode', 'manual');            
set (gca, 'YTick', (min_x:16:max_x));   %put ticks in every 16 nm, using min and max variables
set (gca, 'YGrid', 'on');

%repeat with the y components:

figure; hold;              %make a new figure

for k = 1:nargin-1                %loop over input traces
    clear y_cut;
    trace = varargin{k};        
    plot (trace(:,2),strcat(color(k),'-'));    %plot trace with dots and dashed lines
    y_cut(1:length(trace(:,4))) = NaN;
    for i=1:length(trace(:,4))
       if (trace(i,6) >= usage) 
           y_cut(i) = trace(i,4);
       end
    end
    plot (y_cut,color(k));                %plot fit with lines
    if (min(trace(:,2)) <= min_y )
       min_y = min(trace(:,2));           %find min and max of all traces
    end
    if (max(trace(:,2)) >= max_y )
       max_y = max(trace(:,2));
    end
end

set (gca,'YTickMode', 'manual');            
set (gca, 'YTick', (min_y:16:max_y));   %put ticks in every 16 nm, using min and max variables
set (gca, 'YGrid', 'on');

end 
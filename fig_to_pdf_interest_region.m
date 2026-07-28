function fig_to_png_pdf_interest_region(write_pdf)
%JS 220204
% MATLAB function to write all fig files selected to png with the option to
% also do pdf. Is useful for trying to present rather than having to open
% and save every individual figure

%LL 20240706
%removed grid line

%LL 20240709 Zoom in
if nargin < 1
    write_pdf = 0;
else
    write_pdf = 1;
end

[fnames, dir] = uigetfile({'*.fig'},'MultiSelect','on');
%Multiselect allows user to grab more than one file and parse accordingly

% if only one file is selected, need to package into a cell array so that
% the for loop will work
if class(fnames) == 'char'
    fnames = {fnames};
end

% loop through each selected file, assuming they are all in the same
% directory
    
for i = 1:length(fnames)
    try
        fname = fnames{i};
    catch
        fprintf("Finished combining xy and yx traces. Aborting... \n")
        return
    end
    % new option to grab the data from yx and put it on the original
    % figure, but only if it is included in the set selected already
    fyxname = strcat(fnames{i}(1:end-9), 'yx', fnames{i}(end-9:end));
    idx = find(contains(fnames,fyxname));
    if ~isempty(idx)
        figgyx = openfig(strcat(dir,fyxname));
        axchild = get(gca,'Children');
        yx = findobj(axchild, 'DisplayName', 'data1');
        X = yx.XData; Y = yx.YData;

        close(figgyx);
        fnames(idx) = [];
        
    end
    
    

    figgs = openfig(strcat(dir,fname));
    % only if yx is included should you copy it over
    if ~isempty(idx)
        axchild = get(gca,'Children');
        delete(findobj(axchild, 'DisplayName', 'data2'));
        hold on
        plot(X,Y,'Color',[0.3,0.3,0.9],'LineStyle','-','LineWidth',0.5,'MarkerSize',6,'Marker','none','MarkerFaceColor','None','DisplayName','data2')
    end
    hold on
    XL = max(xlim);
    YL = max(ylim);
    %handles.display.gridLines = 0; % ignore this
    plot([0.05*XL,0.05*XL],[0.6*YL,0.6*YL+200],'k','LineWidth',1,'tag','200 nm') %200ms, 200nm - Change based on need
    plot([0.05*XL,0.05*XL+0.2],[0.6*YL,0.6*YL],'k','LineWidth',1,'tag','200 ms')
    text(0.05*XL,0.6*YL+240,'200 nm')
    text(0.05*XL+0.21,0.6*YL,'200 ms') % 200ms
    yticks('auto') %haha one line fix it
    ylabel('')
    xlabel('')
    box off 
    grid off
        
    
  
    all_lines = findobj(gca, 'Type', 'Line');
    yx = [];
    max_len = 0;

    for k = 1:length(all_lines)
        this_line = all_lines(k);
        
        % Check 1: Ignore the "Off-Axis" line (which we named 'data2')
        is_off_axis = strcmp(this_line.DisplayName, 'data2');
        
        % Check 2: Ignore Scale Bars (which have Tags)
        is_scale_bar = ~isempty(this_line.Tag);
        
        % If it's NOT off-axis and NOT a scale bar, it's our On-Axis data
        if ~is_off_axis && ~is_scale_bar
            % Safety: Pick the longest line just in case there's noise
            if length(this_line.XData) > max_len
                max_len = length(this_line.XData);
                yx = this_line;
            end
        end
    end

    % --- INSET PLOTTING ---
    if ~isempty(yx)
        yx.LineWidth = 0.5; % Make the main trace thin
        X = yx.XData; 
        Y = yx.YData;
        
        % Create inset axes
        axes('position',[.65 .3 .2 .2])
        box on 
        
        % --- TIME-BASED SLICING (Auto-detects duration) ---
        % This works even if lengths vary (e.g., 1s vs 2s traces)
        t_start = min(X);
        t_total = max(X) - t_start;
        
        % User Settings: 30% to 40% of the trace time
        t_cut_start = t_start + (0.30 * t_total); 
        t_cut_end   = t_start + (0.40 * t_total);
        
        % Find indices within this time window
        indices = find(X >= t_cut_start & X <= t_cut_end);
        
        if ~isempty(indices)
            % Plot with the same color as the original line
            plot(X(indices), Y(indices), 'LineWidth', 0.5, 'Color', yx.Color) 
            axis tight
            % Optional: Remove ticks for cleaner look
            % set(gca,'XTick',[],'YTick',[])
        end
    else
        fprintf("Warning: On-axis data trace not found in %s.\n", fname);
    end
    
 
    
    saveas(figgs, strcat(dir,fname(1:end-4),'.pdf'))
    fprintf("Wrote " + strcat(fname(1:end-4)) + '.pdf' + '\n')   
    
    close(figgs)

end

end


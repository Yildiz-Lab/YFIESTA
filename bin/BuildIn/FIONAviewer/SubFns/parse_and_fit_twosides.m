function fitted_trace = parse_and_fit_twosides(x,y,usage,thresh,window)

%test function for running the step fitter

%usage: x, y, usage, and thresh are suppied as HORIZTONAL vectors
%for now, it returns a vertical six-column vector, in the traditional
%format

%to adapt to FIONAviewer down the line, will return horizontal vectors for
%replacing the input vectors

%breaks up a single trace (actually columns 1, 2, and 6, passed as HORIZONTAL vectors) into
%sub-traces based on usage. Each piece is de-NaN-ified and fit using the
%SIC fitter. Returned chunks from the SIC fitter are stiched back into a
%step vector. For now return a 6-column trace including a changepoint. In
%future could just return off and on-axis steps

%for off-axis steps: First, I made an off-axis fit using the on-axis cp and
%then found extra steps in the residual, but that appeared to overfit. For
%this version, "two sides", I fit both traces separately. Then if a cp is
%seen within "window" of the off-axis, take it out, and use that on-axis
%step as a guide for both axes

%for making dwell locations good, add changepoints right after usage
%windows and fit as one long mean

%for cosmetic purposes might just want to only plot using a user-supplied
%usage threshold...but that will happen in plot_traces

usage_mask = (usage >= thresh);

counter = 0;
if usage_mask(1) == 1
    start = 1;
else
    start = 0;
end
sub_traces_x = [];
sub_traces_y = [];
sub_traces_ind = [];
for (i=2:length(usage_mask))   %gonna do some ugly matlab jujitsu here...sorry mom!
    if usage_mask(i) == 1 && usage_mask(i-1) == 0  %if we transition from non-use to use, set start of sub_trace
        start = i;
    end
    if usage_mask(i) == 0 && usage_mask(i-1) == 1  %if we transition from use to non-use, make the sub_traces
        counter = counter + 1;  
        sub_traces_x{counter} = x(start:i-1);
        sub_traces_y{counter} = y(start:i-1);
        sub_traces_ind{counter} = (start:i-1);
    end
end

%now strip out the NaNs from all the sub_traces
for (i=1:length(sub_traces_x))
    for (j = 4:length(sub_traces_x{i}) )
        if (isnan(sub_traces_x{i}(j)))
            sub_traces_x{i}(j) = mean(sub_traces_x{i}(j-3:j-1));        %use average of previous three points as the fake value
            sub_traces_y{i}(j) = mean(sub_traces_y{i}(j-3:j-1));
        end
    end
    disp (sum(isnan(sub_traces_x{i})));
end

% fitted_trace.x = sub_traces_x;
% fitted_trace.y = sub_traces_y;
% fitted_trace.ind = sub_traces_ind;


%so...now we have a set of sub traces with NaNs removed
%run through step fitter, get changepoints
%then apply cps to off-axis trace, subtract fit and then fit residuals

for (i = 1:length(sub_traces_x))  %fit each sub trace with the SIC step fitter
    fit = SICstepFinder(sub_traces_x{i});       %run the step finder on both the on and off axes
    changepoints_x = [ 0  ( ( fit.StepFit(1:length(fit.StepFit)-1) - fit.StepFit(2:length(fit.StepFit)) ) ~= 0) ]; %make a changepoints array (1=step, 0=no step)
    yfit = SICstepFinder(sub_traces_y{i});
    changepoints_y = [ 0  ( ( yfit.StepFit(1:length(yfit.StepFit)-1) - yfit.StepFit(2:length(yfit.StepFit)) ) ~= 0) ]; %make a cp array of the off-axis too
    %now harmonize the two chanegpoints arrays, using for loops (sorry!)
    for (j = 1:length(changepoints_x))
       if changepoints_x(j) == 1
           changepoints_y(max(1,j-window):min(length(changepoints_x),j+window)) = 0;  %use only information from the on-axis within "window" of each on-axis step
       end
    end

    changepoints_all{i} = (changepoints_x + changepoints_y) > 0;
    last = 1;   %Use the main changeopoints array to construct the fit of the sub_trace
    fit_traces_x{i} = NaN(1,length(changepoints_all{i}));
    fit_traces_y{i} = NaN(1,length(changepoints_all{i}));
    for (j=2:length(changepoints_all{i}))
       if changepoints_all{i}(j) == 1
           fit_traces_x{i}(last:j-1) = nanmean(sub_traces_x{i}(last:j-1));
           fit_traces_y{i}(last:j-1) = nanmean(sub_traces_y{i}(last:j-1));
           last = j;
       end
    end
    fit_traces_x{i}(last:j) = nanmean(sub_traces_x{i}(last:j));
    fit_traces_y{i}(last:j) = nanmean(sub_traces_y{i}(last:j));
    disp (sum(isnan(fit_traces_x{i})));
    disp (sum(isnan(fit_traces_y{i})));
end

%pre-build four of the columns of the trace
fitted_trace = NaN(length(x),6);
fitted_trace(1:length(x),[1 2 5 6]) = [x' y' zeros(length(x),1) usage' ];

%then use the sub_traces and their indices to write out a 6-column trace
%object
%future versions of this function will re-structure output based on
%FIONAviewer's handles object, right?
for (i=1:length(sub_traces_x))
    fitted_trace(sub_traces_ind{i},3:6) = [ fit_traces_x{i}' fit_traces_y{i}' changepoints_all{i}' usage(sub_traces_ind{i})'];
    pre_cp_ind = sub_traces_ind{i}(1) - 1; %after each usage zone, add in a changepoint to make sure fitted values can't change once its keyed in (might ditch this some point)
    post_cp_ind = sub_traces_ind{i}(length(sub_traces_ind{i})) + 1;
    if (pre_cp_ind > 1)
        fitted_trace(pre_cp_ind,5) = 1;
    end
    if (post_cp_ind < length(x) )
        fitted_trace(post_cp_ind,5) = 1;
    end
end








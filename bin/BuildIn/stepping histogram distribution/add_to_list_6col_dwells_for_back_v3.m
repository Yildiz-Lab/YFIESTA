function [forward,backward] = add_to_list_6col_dwells_for_back_v3(trace,tchoice,framerate)
xsteps = trace(:,4); xsteps = xsteps(~isnan(xsteps));
ysteps = trace(:,3); ysteps = ysteps(~isnan(ysteps));
stepsis = [ysteps xsteps];
L = length(stepsis);
forsteps = zeros(L,1);
for i = 2:L
    if ysteps(i)> ysteps(i-1)
        forsteps(i) = 1;
    elseif ysteps(i)< ysteps(i-1)
        forsteps(i) = -1;
    else
        forsteps(i)= 0;
    end
end

% Indices are shifted to start at 1 rather than the index that tchoice
% recognizes if one were to not start index counting at 1 from the trace.
% To counteract that, we find the first index that is not NaN
shiftidx = find(~isnan(trace(:,3)));
shiftidx = shiftidx(1)-1;

dwell_for=[];
dwell_back=[];
cnt=0; %changed count to 1, but should we be careful with our previous function?
for i=1:length(forsteps)
    %check that the beginning point is within that region
    % if we want to do end, then check i rather than i-1
    if ismember(i-1+shiftidx,tchoice)
        if(forsteps(i)==0)
            cnt=cnt+1;
        elseif(forsteps(i)==1)
            dwell_for=[ dwell_for cnt];
            cnt=1;
        else
            dwell_back=[ dwell_back cnt];
            cnt=1;
        end
    end

end
forward = dwell_for .* framerate;
forward = forward';
backward = dwell_back .* framerate;
backward = backward';


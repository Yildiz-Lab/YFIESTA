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
dwell_for=[];
dwell_back=[];
cnt=1; %changed count to 1, but should we be careful with our previous function?
for i=1:length(forsteps)
    %check that the beginning point is within that region
    % if we want to do end, then check i rather than i-1
    if ismember(i-1,tchoice) 
        if(forsteps(i)==0)
            cnt=cnt+1;
        elseif(forsteps(i)==1)
            dwell_for=[ dwell_for cnt]; cnt=1;
        else
            dwell_back=[ dwell_back cnt]; cnt=1;
        end
    end

end
forward = dwell_for .* framerate;
forward = forward';
backward = dwell_back .* framerate;
backward = backward';



function [forward,backward] = add_to_list_6col_dwells_for_back(trace,framerate)
xsteps = trace(:,4);
xsteps = xsteps(~isnan(xsteps));
ysteps = trace(:,3);
ysteps = ysteps(~isnan(ysteps));
stepsis = [ysteps xsteps];
L = length(stepsis);
forsteps = zeros(L,1);
    for i = 2:L
        if ysteps(i)> ysteps(i-1)
            forsteps(i) = 1;
        elseif ysteps(i)< ysteps(i-1)
            forsteps(i) = -1;
        else forsteps(i)= 0;
        end
    end
dwell_for=[];
dwell_back=[];
cnt=0;
    for i=1:length(forsteps)
        if(forsteps(i)==0)
            cnt=cnt+1;
        elseif(forsteps(i)==1)
                dwell_for=[ dwell_for cnt]; cnt=0;
                    
        else
                dwell_back=[ dwell_back cnt]; cnt=0;
        end
          
    end
    forward = dwell_for .* framerate;
    forward = forward';
    backward = dwell_back .* framerate;
    backward = backward';
end

function [forward,backward] = add_to_list_6col_dwells_for_back(trace,framerate)
xsteps = trace(:,4);
xsteps = xsteps(~isnan(xsteps));
ysteps = trace(:,3);
ysteps = ysteps(~isnan(ysteps));
stepsis = [ysteps xsteps];
L = length(stepsis);

if length(framerate) < 2 %if framerate is an integer, then make it an array of times
    time = framerate*0:length(trace(:,3))-1;
else
    time = framerate;
end

forsteps = zeros(L,1);
    for i = 2:L
        if ysteps(i)> ysteps(i-1)
            forsteps(i) = 1;
        elseif ysteps(i)< ysteps(i-1)
            forsteps(i) = -1;
        else forsteps(i)= 0;
        end
    end
dwell_for=[]; dwell_for_t=[];
dwell_back=[]; dwell_back_t=[];
cnt=0; t=0;
    for i=1:length(forsteps)
        if(forsteps(i)==0)
            cnt=cnt+1;
            %JS Edit 2024/03/07 % rather than count add dt
            t = t + time(i+1)-time(i);
        elseif(forsteps(i)==1)
                dwell_for=[ dwell_for cnt]; cnt=0;
                dwell_for_t = [dwell_for_t; t]; t=0;
                    
        else
                dwell_back=[ dwell_back cnt]; cnt=0;
                dwell_back_t = [dwell_back_t; t]; t=0;
        end
          
    end
    % forward = dwell_for .* framerate;
    % forward = forward';
    % backward = dwell_back .* framerate;
    % backward = backward';
    
    % JS Edit 2024/03/07 for dt MINFLUX
    forward = dwell_for_t;
    backward = dwell_back_t;

end
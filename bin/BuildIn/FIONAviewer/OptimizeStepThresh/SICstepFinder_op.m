function StepStatistics=SICstepFinder(x)

% NOTE THAT IN FIESTA > FIONA this is unused. Instead look into parent
% folder for same filename

%% 
%----Info------------------------------------------------------------------
% Implementation of SIC Step-Finding Routine as published by B. Kalafut and 
% K. Visscher, 2008, Computer Physics Communications
% Originally written May 2011 by Thomas Bilyard and Sheng-Min Shih, 
% University of California, Berkeley
%
% Modified to minimize the Schwarz Information Criterion (SIC) October 2012
% by Vladislav Belyy
%
%
%-------------Step-finding input parameters--------------------------------
% x - input data
%
%--------------Outputs in 'StepStatistics' structure-----------------------
% StepsChi2 - Chi2 decrease for each step added
% NumberOfStepsFound - number of steps found
% StepSizeStats - size of each step
% StepLengths - length of each step
% StepFit - stepwise fit of data
% FinalSIC - final value of the Schwarz Information Criterion
%--------------------------------------------------------------------------

%% Initialization routine

% if the correction factor is included, one-point steps can not be fitted
IfCorrect = 1;
MinPointsInStep = 2;

%first and last transitions defined at 0 and end of data array
NP=length(x);
StepsInd=[0 NP];

% initialize Fit array
FitPos=zeros(1,NP);

% Calculate x-squared;
x_sq=x.*x;

%calculate forwards/backwards arrays
RCumx_sq=cumsum(x_sq(1:end));
LCumx_sq=fliplr(cumsum(fliplr(x_sq(1:end))));
RCumx=cumsum(x(1:end));
LCumx=fliplr(cumsum(fliplr(x(1:end))));
R_N=1:1:length(x);
L_N=length(x):-1:1;

% calculate DeltaChi2 for each step location
Metric = zeros(1,NP); % the last point is not used except for plotting
Metric(1:end-1) = UpdateMetric( RCumx, LCumx, R_N, L_N );

NoFittedSteps=1;    % current number of steps found =1

FitSquidual=(std(x))^2*NP; % Chi2 = Variance*NP = total noise

% need to get the first step running
%tempDeltaQ = - MinDeltaQ - 1; %#ok<NASGU>

% Calculate starting Schwarz Information Criterion
SIC_new = (NoFittedSteps+2)*log(NP)+NP*log(FitSquidual/NP);
SIC_old = SIC_new;

%% Main stepfitter loop
 while(SIC_new <= SIC_old) 
    
    SIC_old = SIC_new;
    
    % Save current params in case this iteration is the last one and the
    % proposed step gets rejected in the end
    RCumx_old = RCumx;
    R_N_old = R_N;
    LCumx_old = LCumx;
    L_N_old = L_N;
    StepsInd_old = StepsInd;
    FitPos_old = FitPos;
    
    [~, StepInd]=min(Metric);
    
    % Find the correct place to insert StepInd into StepsInd
    StepLocs=find((StepsInd-StepInd)>0);
    StepLoc=StepLocs(1);
    
    PrevInd=StepsInd(StepLoc-1);
    NextInd=StepsInd(StepLoc);
    
    % only accept steps longer than MinPointsInStep
    if ( (StepInd-PrevInd)>=MinPointsInStep && ...
            (NextInd-StepInd)>=MinPointsInStep)
        
        StepsInd=[StepsInd(1:StepLoc-1) StepInd StepsInd(StepLoc:end)]; 
        
        % assign mean of surrounding periods as step levels
        FitPos(PrevInd+1:StepInd)=mean(x(PrevInd+1:StepInd));
        FitPos(StepInd+1:NextInd)=mean(x(StepInd+1:NextInd));

        Metric(StepInd)=0; % do not allow duplicate steps at the same point
        
        NoFittedSteps=NoFittedSteps+1;  % increase valid step tally
        
        %correcting the previous transition points
        if(PrevInd==0)%means no previous dwell
            
            LCumx(PrevInd+1:StepInd)= ...
                fliplr(cumsum(fliplr(x(PrevInd+1:StepInd))));
            L_N(PrevInd+1:StepInd)=StepInd-PrevInd:-1:1;
            LCumx_sq(PrevInd+1:StepInd)= ...
                fliplr(cumsum(fliplr(x_sq(PrevInd+1:StepInd))));
           
            Metric(PrevInd+1:StepInd-1) = UpdateMetric(...
                RCumx(PrevInd+1:StepInd), LCumx(PrevInd+1:StepInd), ...
                R_N(PrevInd+1:StepInd), L_N(PrevInd+1:StepInd) );
            
        else % There is a previous dwell; correct previous step's position
            Prev2Ind=StepsInd(StepLoc-2);
            tempRCumx=cumsum(x(Prev2Ind+1:StepInd));
            tempR_N=1:1:StepInd-Prev2Ind;
            tempLCumx=fliplr(cumsum(fliplr(x(Prev2Ind+1:StepInd))));
            tempL_N=StepInd-Prev2Ind:-1:1;
                    
            tempMetric = UpdateMetric( tempRCumx, tempLCumx, ...
                                                    tempR_N, tempL_N );
            
            [~, PrevInd] = min(tempMetric);
            PrevInd = PrevInd + Prev2Ind;
                        
      
            %update forwards arrays
            RCumx(Prev2Ind+1:PrevInd)=cumsum(x(Prev2Ind+1:PrevInd));
            RCumx(PrevInd+1:StepInd)=cumsum(x(PrevInd+1:StepInd));
            R_N(Prev2Ind+1:PrevInd)=1:1:PrevInd-Prev2Ind;
            R_N(PrevInd+1:StepInd)=1:1:StepInd-PrevInd;
            
            LCumx(Prev2Ind+1:PrevInd)= ...
                fliplr(cumsum(fliplr(x(Prev2Ind+1:PrevInd))));
            LCumx(PrevInd+1:StepInd)= ...
                fliplr(cumsum(fliplr(x(PrevInd+1:StepInd))));
            L_N(Prev2Ind+1:PrevInd)=PrevInd-Prev2Ind:-1:1;
            L_N(PrevInd+1:StepInd)=StepInd-PrevInd:-1:1;
            
            RCumx_sq(Prev2Ind+1:PrevInd)= ...
                cumsum(x_sq(Prev2Ind+1:PrevInd));
            RCumx_sq(PrevInd+1:StepInd)= ...
                cumsum(x_sq(PrevInd+1:StepInd));
            LCumx_sq(Prev2Ind+1:PrevInd)= ...
                fliplr(cumsum(fliplr(x_sq(Prev2Ind+1:PrevInd))));
            LCumx_sq(PrevInd+1:StepInd)= ...
                fliplr(cumsum(fliplr(x_sq(PrevInd+1:StepInd))));

        
            % update DeltaChi2 array
            Metric(PrevInd+1:StepInd-1) = UpdateMetric( ...
                RCumx(PrevInd+1:StepInd), LCumx(PrevInd+1:StepInd), ...
                R_N(PrevInd+1:StepInd), L_N(PrevInd+1:StepInd) );
            Metric(Prev2Ind+1:PrevInd-1) = UpdateMetric( ...
                RCumx(Prev2Ind+1:PrevInd), LCumx(Prev2Ind+1:PrevInd), ...
                R_N(Prev2Ind+1:PrevInd), L_N(Prev2Ind+1:PrevInd) );
            
            Metric(PrevInd)=0;
            StepsInd(StepLoc-1)=PrevInd;
            
            % assign mean of surrounding periods as step levels
            FitPos(Prev2Ind+1:PrevInd)=mean(x(Prev2Ind+1:PrevInd));
            FitPos(PrevInd+1:StepInd)=mean(x(PrevInd+1:StepInd));
            
        end
        
        %correcting the next transition points
        if(NextInd==NP) %means no next dwell
            RCumx(StepInd+1:NextInd)=cumsum(x(StepInd+1:NextInd));
            R_N(StepInd+1:NextInd)=1:1:NextInd-StepInd;
            RCumx_sq(StepInd+1:NextInd)=cumsum(x_sq(StepInd+1:NextInd));
         
            Metric(StepInd+1:NextInd-1) = UpdateMetric( ...
                RCumx(StepInd+1:NextInd), LCumx(StepInd+1:NextInd), ...
                R_N(StepInd+1:NextInd), L_N(StepInd+1:NextInd) );
            
            if(StepInd~=NP-1)
                Metric(StepInd+1) = 0;
            end
            Metric(NextInd-1) = 0;
            
        else % There is a subsequent dwell; correct next step's position
            Next2Ind=StepsInd(StepLoc+2); % after inserting the current step
            
            tempRCumx=cumsum(x(StepInd+1:Next2Ind));
            tempR_N=1:1:Next2Ind-StepInd;
            tempLCumx=fliplr(cumsum(fliplr(x(StepInd+1:Next2Ind))));
            tempL_N=Next2Ind-StepInd:-1:1;
            
            tempMetric = UpdateMetric( tempRCumx, tempLCumx, ...
                                                 tempR_N, tempL_N );
            
            [~, NextInd] = min(tempMetric);
            NextInd = NextInd + StepInd; 
            %+1 because MetricPart excludes the 1st and last transition already
            
            %update backward arrays
            RCumx(StepInd+1:NextInd)=cumsum(x(StepInd+1:NextInd));
            RCumx(NextInd+1:Next2Ind)=cumsum(x(NextInd+1:Next2Ind));
            R_N(StepInd+1:NextInd)=1:1:NextInd-StepInd;
            R_N(NextInd+1:Next2Ind)=1:1:Next2Ind-NextInd;
            
            LCumx(StepInd+1:NextInd)= ...
                fliplr(cumsum(fliplr(x(StepInd+1:NextInd))));
            LCumx(NextInd+1:Next2Ind)= ...
                fliplr(cumsum(fliplr(x(NextInd+1:Next2Ind))));
            L_N(StepInd+1:NextInd)=NextInd-StepInd:-1:1;
            L_N(NextInd+1:Next2Ind)=Next2Ind-NextInd:-1:1;
            
            RCumx_sq(StepInd+1:NextInd)= ...
                cumsum(x_sq(StepInd+1:NextInd));
            RCumx_sq(NextInd+1:Next2Ind)= ...
                cumsum(x_sq(NextInd+1:Next2Ind));
            LCumx_sq(StepInd+1:NextInd)= ...
                fliplr(cumsum(fliplr(x_sq(StepInd+1:NextInd))));
            LCumx_sq(NextInd+1:Next2Ind)= ...
                fliplr(cumsum(fliplr(x_sq(NextInd+1:Next2Ind))));

            % update DeltaChi2 array
            Metric(StepInd+1:NextInd-1) = UpdateMetric( ...
                RCumx(StepInd+1:NextInd), LCumx(StepInd+1:NextInd), ...
                R_N(StepInd+1:NextInd), L_N(StepInd+1:NextInd) );
            Metric(NextInd+1:Next2Ind-1) = UpdateMetric( ...
                RCumx(NextInd+1:Next2Ind), LCumx(NextInd+1:Next2Ind), ...
                R_N(NextInd+1:Next2Ind), L_N(NextInd+1:Next2Ind) ); 
            
            Metric(NextInd)=0;
            StepsInd(StepLoc+1)=NextInd;
            
            % assign mean of surrounding periods as step levels
            FitPos(StepInd+1:NextInd)=mean(x(StepInd+1:NextInd));
            FitPos(NextInd+1:Next2Ind)=mean(x(NextInd+1:Next2Ind));
            
        end
        
        if(IfCorrect)
            CorrectionFactor=R_N(StepsInd(2:end))./(R_N(StepsInd(2:end))-1);
        else
            CorrectionFactor = 1;
        end
        
         % calculate variance of data from fit
        FitSquidual=sum((RCumx_sq(StepsInd(2:end))-(RCumx( ...
            StepsInd(2:end))).^2./R_N(StepsInd(2:end))).*CorrectionFactor);
        
    else  % Step is shorter than minimum allowed length; ignore it
        Metric(StepInd)=0;
        disp('a short-dwell step found and ignored')
    end
    
    % Update the Schwarz Information Criterion
    SIC_new = (NoFittedSteps+2)*log(NP)+NP*log(FitSquidual/NP);
    
 end
 
%% Prepare stepfitter results for output
 
 disp( ...
     strcat('Final Schwarz Information Criterion: ', num2str(SIC_old)));

%calculate Step Statistics
StepStatistics.StepsChi2=-RCumx_old(StepsInd_old(2:end-1)).* ...
    RCumx_old(StepsInd_old(2:end-1))./R_N_old(StepsInd_old(2:end-1))- ...
    LCumx_old(StepsInd_old(2:end-1)+1).*LCumx_old(StepsInd_old(2:end-1)+1)./ ...
    L_N_old(StepsInd_old(2:end-1)+1)+(RCumx_old(StepsInd_old(2:end-1))+ ...
    LCumx_old(StepsInd_old(2:end-1)+1)).*(RCumx_old(StepsInd_old(2:end-1))+ ...
    LCumx_old(StepsInd_old(2:end-1)+1))./(R_N_old(StepsInd_old(2:end-1))+ ...
    L_N_old(StepsInd_old(2:end-1)+1));
StepStatistics.NumberOfStepsFound=length(StepsInd_old)-2;
StepStatistics.StepSizeStats=FitPos_old(StepsInd_old(2:end-1)+1)- ...
    FitPos_old(StepsInd_old(2:end-1));
StepStatistics.StepLengths=StepsInd_old(2:end)-StepsInd_old(1:end-1);
StepStatistics.StepFit=FitPos_old;
StepStatistics.FinalSIC = SIC_old;


%% Sub-routines
function newMetric = UpdateMetric(RCumx,LCumx,R_N,L_N)

newMetric = -RCumx(1:end-1).*RCumx(1:end-1)./R_N(1:end-1)...
            -LCumx(2:end  ).*LCumx(2:end  )./L_N(2:end  )...
            +(RCumx(1:end-1)+LCumx(2:end  )).*(RCumx(1:end-1)+ ...
            LCumx(2:end ))./(R_N(1:end-1)+L_N(2:end  ));



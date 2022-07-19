% ShihBilyardStepFinder_TDV: plug-in for TrapDataViewer
% ShihBilyardStepFinder20110805F: it only uses the correction factor for calculate total squidual, no for fitting
% the correction factor can be switched off by changing IfCorrect to 0
% ShihBilyardStepFinder20110802CF: estimate noise from power spectrum if input ExpNoise2=0
% ShihBilyardStepFinder20110801C: with the correction factor.
% No one-point step would be fitted
function StepStatistics=ShihBilyardStepFinder_TDV(x,MaxNoSteps,ExpNoise2,MinDeltaQ,MinPointsInStep,DoPlot)

%----Info------------
% Least-Squidual(TM) Step-Finding Routine
% Originally written May 2011 by Thomas Bilyard and Sheng-Min Shih, University of California, Berkeley
%-------------Step-finding input parameters--------------
% x - data
% MaxNoSteps - Maximum number of steps to find, it's a terminating condition
%              MaxNoSteps=0 means there is no limit on the number of steps
% ExpNoise2 - Expected Noise2 (squared noise). The program will stop
%             fitting after the squidual reach this value. If ExpNoise2=0,
%             the program will determine noise from the power spectrum of
%             residual (recursively)
% MinDeltaQ - Minimum delta squidual (the absolute value of it!!!), it's a terminating condition
% MinPointsInStep - reject all steps lasting less than MinPointsInStep points;
% DoPlot - 0: do not plot anything
%          1: plot average squidual per pt and expected noise2
%          2: plot fitted steps and everything described before
%          3: plot delta squidual for each pt and everything described before
%--------------Outputs in 'StepStatistics' structure-------------------
% StepsChi2 - Chi2 decrease for each step added
% NumberOfStepsFound - number of steps found
% StepSizeStats - size of each step
% StepLengths - length of each step
% StepFit - stepwise fit of data
%------------------------------------------

% Correction factor for Squidual: 1->True, 0->False
IfCorrect = 1; 

% if the correction factor is included, one-point steps can not be fitted
if (IfCorrect)
    if(MinPointsInStep<2)
        MinPointsInStep = 2;
    end
end

%Minimum Delta Squidual should be at least 0
if(MinDeltaQ<0)
    MinDeltaQ = 0;
end

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

% before correction: Metric(1:end-1)=-RCumx(1:end-1).*RCumx(1:end-1)./(R_N(1:end-1)-IfCorrect)-LCumx(2:end).*LCumx(2:end)./(L_N(2:end)-IfCorrect)+(RCumx(1:end-1)+LCumx(2:end)).*(RCumx(1:end-1)+LCumx(2:end))./(R_N(1:end-1)+L_N(2:end)-IfCorrect);
% Metric(end)=0; % add a point for plotting later

NoFittedSteps=1;    % current number of steps found =1

FitSquidual=(std(x))^2*NP; % Chi2 = Variance*NP = total noise

fSamp=1; %Important: there is not input of fSamp, just use 1 and it should affect the result

if(ExpNoise2==0)
    [Spectrum, ff] = OneSidedSpectrumSS(x, fSamp);
    tempsize = size(Spectrum);
    startingpoint = ceil(tempsize(1)*tempsize(2)/5);
    ExpFinalSquidual = mean(Spectrum(startingpoint:end)) *fSamp/2 *NP;
    disp(strcat('estimated noise2 = ',num2str(ExpFinalSquidual/NP)))
else
    ExpFinalSquidual = ExpNoise2*NP; 
end
       
if (DoPlot)
    figure(11)
    plot(NoFittedSteps,FitSquidual/NP,'ro')
    %             legend('FitChi2')
    hold on

    line([2 NoFittedSteps],[ExpFinalSquidual/NP ExpFinalSquidual/NP])
    legend('Fit Squidual','Expected Final Squidual')
    hold on

    if (DoPlot>=2)
        t=1:1:NP;
        figure(10)
        plot(t,x,'b.',t,FitPos,'r-')
        
        if (DoPlot>=3)

            figure(12)
            plot(t,Metric,'r-')

        end
        pause(0.100)
    end

end

% need to get the first step running
tempDeltaQ = - MinDeltaQ - 1;

% find all further steps
while (FitSquidual>ExpFinalSquidual & ( NoFittedSteps<MaxNoSteps | NoFittedSteps==0 ) & (-tempDeltaQ)>MinDeltaQ ) 
    
    [tempDeltaQ StepInd]=min(Metric);
    
    % Find the correct place to insert StepInd into StepsInd
    StepLocs=find((StepsInd-StepInd)>0);
    StepLoc=StepLocs(1);
    
    PrevInd=StepsInd(StepLoc-1);
    NextInd=StepsInd(StepLoc);
    
%     StepSize=-mean(x(PrevInd+1:StepInd))+mean(x(StepInd+1:NextInd));   %StepSize is not currently used
    
    % only accept steps longer than MinPointsInStep
    if ( (StepInd-PrevInd)>=MinPointsInStep & (NextInd-StepInd)>=MinPointsInStep & (-tempDeltaQ)>MinDeltaQ )
        
        StepsInd=[StepsInd(1:StepLoc-1) StepInd StepsInd(StepLoc:end)]; % it was: StepsInd=sort([StepsInd StepInd]);
        
        % assign mean of surrounding periods as step levels
        FitPos(PrevInd+1:StepInd)=mean(x(PrevInd+1:StepInd));
        FitPos(StepInd+1:NextInd)=mean(x(StepInd+1:NextInd));

        Metric(StepInd)=0;
        
        NoFittedSteps=NoFittedSteps+1;  % increase valid step tally
        
        %correcting the previous transition points
        if(PrevInd==0)%means no previous dwell
            
            LCumx(PrevInd+1:StepInd)=fliplr(cumsum(fliplr(x(PrevInd+1:StepInd))));
            L_N(PrevInd+1:StepInd)=StepInd-PrevInd:-1:1;
            LCumx_sq(PrevInd+1:StepInd)=fliplr(cumsum(fliplr(x_sq(PrevInd+1:StepInd))));
            % before correction: Metric(PrevInd+1:StepInd-1)=-RCumx(PrevInd+1:StepInd-1).*RCumx(PrevInd+1:StepInd-1)./R_N(PrevInd+1:StepInd-1)-LCumx(PrevInd+2:StepInd).*LCumx(PrevInd+2:StepInd)./L_N(PrevInd+2:StepInd)+(RCumx(PrevInd+1:StepInd-1)+LCumx(PrevInd+2:StepInd)).*(RCumx(PrevInd+1:StepInd-1)+LCumx(PrevInd+2:StepInd))./(R_N(PrevInd+1:StepInd-1)+L_N(PrevInd+2:StepInd));
            Metric(PrevInd+1:StepInd-1) = UpdateMetric( RCumx(PrevInd+1:StepInd), LCumx(PrevInd+1:StepInd), R_N(PrevInd+1:StepInd), L_N(PrevInd+1:StepInd) );
            
        else
            Prev2Ind=StepsInd(StepLoc-2);
            tempRCumx=cumsum(x(Prev2Ind+1:StepInd));
            tempR_N=1:1:StepInd-Prev2Ind;
            tempLCumx=fliplr(cumsum(fliplr(x(Prev2Ind+1:StepInd))));
            tempL_N=StepInd-Prev2Ind:-1:1;
            
            tempRCumx_sq=cumsum(x_sq(Prev2Ind+1:StepInd));
            tempLCumx_sq=fliplr(cumsum(fliplr(x_sq(Prev2Ind+1:StepInd))));
            NP2D=StepInd-Prev2Ind; %number of points in these 2 dwells
            
            tempMetric = UpdateMetric( tempRCumx, tempLCumx, tempR_N, tempL_N );
            
            [temp PrevInd] = min(tempMetric);
            PrevInd = PrevInd + Prev2Ind;
                        
%             Chi2After=VarTrans(PrevInd-Prev2Ind);
            
            %update forwards arrays
            RCumx(Prev2Ind+1:PrevInd)=cumsum(x(Prev2Ind+1:PrevInd));
            RCumx(PrevInd+1:StepInd)=cumsum(x(PrevInd+1:StepInd));
            R_N(Prev2Ind+1:PrevInd)=1:1:PrevInd-Prev2Ind;
            R_N(PrevInd+1:StepInd)=1:1:StepInd-PrevInd;
            
            LCumx(Prev2Ind+1:PrevInd)=fliplr(cumsum(fliplr(x(Prev2Ind+1:PrevInd))));
            LCumx(PrevInd+1:StepInd)=fliplr(cumsum(fliplr(x(PrevInd+1:StepInd))));
            L_N(Prev2Ind+1:PrevInd)=PrevInd-Prev2Ind:-1:1;
            L_N(PrevInd+1:StepInd)=StepInd-PrevInd:-1:1;
            
            RCumx_sq(Prev2Ind+1:PrevInd)=cumsum(x_sq(Prev2Ind+1:PrevInd));
            RCumx_sq(PrevInd+1:StepInd)=cumsum(x_sq(PrevInd+1:StepInd));
            LCumx_sq(Prev2Ind+1:PrevInd)=fliplr(cumsum(fliplr(x_sq(Prev2Ind+1:PrevInd))));
            LCumx_sq(PrevInd+1:StepInd)=fliplr(cumsum(fliplr(x_sq(PrevInd+1:StepInd))));

        
            % update DeltaChi2 array
            Metric(PrevInd+1:StepInd-1) = UpdateMetric( RCumx(PrevInd+1:StepInd), LCumx(PrevInd+1:StepInd), R_N(PrevInd+1:StepInd), L_N(PrevInd+1:StepInd) );
            Metric(Prev2Ind+1:PrevInd-1) = UpdateMetric( RCumx(Prev2Ind+1:PrevInd), LCumx(Prev2Ind+1:PrevInd), R_N(Prev2Ind+1:PrevInd), L_N(Prev2Ind+1:PrevInd) );
            
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
            % before correction: Metric(StepInd+1:NextInd-1)=-RCumx(StepInd+1:NextInd-1).*RCumx(StepInd+1:NextInd-1)./R_N(StepInd+1:NextInd-1)-LCumx(StepInd+2:NextInd).*LCumx(StepInd+2:NextInd)./L_N(StepInd+2:NextInd)+(RCumx(StepInd+1:NextInd-1)+LCumx(StepInd+2:NextInd)).*(RCumx(StepInd+1:NextInd-1)+LCumx(StepInd+2:NextInd))./(R_N(StepInd+1:NextInd-1)+L_N(StepInd+2:NextInd));
            Metric(StepInd+1:NextInd-1) = UpdateMetric( RCumx(StepInd+1:NextInd), LCumx(StepInd+1:NextInd), R_N(StepInd+1:NextInd), L_N(StepInd+1:NextInd) );
            
            if(StepInd~=NP-1)
                Metric(StepInd+1) = 0;
            end
            Metric(NextInd-1) = 0;
            
        else
            Next2Ind=StepsInd(StepLoc+2); % this is after inserting the current step
            
            tempRCumx=cumsum(x(StepInd+1:Next2Ind));
            tempR_N=1:1:Next2Ind-StepInd;
            tempLCumx=fliplr(cumsum(fliplr(x(StepInd+1:Next2Ind))));
            tempL_N=Next2Ind-StepInd:-1:1;
            
            tempRCumx_sq=cumsum(x_sq(StepInd+1:Next2Ind));
            tempLCumx_sq=fliplr(cumsum(fliplr(x_sq(StepInd+1:Next2Ind))));
            NP2D=Next2Ind-StepInd; %number of points in these 2 dwells
            
            tempMetric = UpdateMetric( tempRCumx, tempLCumx, tempR_N, tempL_N );
            
            [temp NextInd] = min(tempMetric);
            NextInd = NextInd + StepInd; %+1 because MetricPart exclude the 1st and last transition already
            
            %update backward arrays
            RCumx(StepInd+1:NextInd)=cumsum(x(StepInd+1:NextInd));
            RCumx(NextInd+1:Next2Ind)=cumsum(x(NextInd+1:Next2Ind));
            R_N(StepInd+1:NextInd)=1:1:NextInd-StepInd;
            R_N(NextInd+1:Next2Ind)=1:1:Next2Ind-NextInd;
            
            LCumx(StepInd+1:NextInd)=fliplr(cumsum(fliplr(x(StepInd+1:NextInd))));
            LCumx(NextInd+1:Next2Ind)=fliplr(cumsum(fliplr(x(NextInd+1:Next2Ind))));
            L_N(StepInd+1:NextInd)=NextInd-StepInd:-1:1;
            L_N(NextInd+1:Next2Ind)=Next2Ind-NextInd:-1:1;
            
            RCumx_sq(StepInd+1:NextInd)=cumsum(x_sq(StepInd+1:NextInd));
            RCumx_sq(NextInd+1:Next2Ind)=cumsum(x_sq(NextInd+1:Next2Ind));
            LCumx_sq(StepInd+1:NextInd)=fliplr(cumsum(fliplr(x_sq(StepInd+1:NextInd))));
            LCumx_sq(NextInd+1:Next2Ind)=fliplr(cumsum(fliplr(x_sq(NextInd+1:Next2Ind))));

            % update DeltaChi2 array
            Metric(StepInd+1:NextInd-1) = UpdateMetric( RCumx(StepInd+1:NextInd), LCumx(StepInd+1:NextInd), R_N(StepInd+1:NextInd), L_N(StepInd+1:NextInd) );
            Metric(NextInd+1:Next2Ind-1) = UpdateMetric( RCumx(NextInd+1:Next2Ind), LCumx(NextInd+1:Next2Ind), R_N(NextInd+1:Next2Ind), L_N(NextInd+1:Next2Ind) ); 
            
            Metric(NextInd)=0;
            StepsInd(StepLoc+1)=NextInd;
            
            % assign mean of surrounding periods as step levels
            FitPos(StepInd+1:NextInd)=mean(x(StepInd+1:NextInd));
            FitPos(NextInd+1:Next2Ind)=mean(x(NextInd+1:Next2Ind));
            
        end
        
        % calculate variance of data from fit
        if(IfCorrect)
            CorrectionFactor=R_N(StepsInd(2:end))./(R_N(StepsInd(2:end))-1);
        else
            CorrectionFactor = 1;
        end
        FitSquidual=sum((RCumx_sq(StepsInd(2:end))-(RCumx(StepsInd(2:end))).^2./R_N(StepsInd(2:end))).*CorrectionFactor);

        if (DoPlot)% plot current steps fitted
            figure(11)
            plot(NoFittedSteps,FitSquidual/NP,'ro')
%             legend('FitChi2')
            hold on

            line([2 NoFittedSteps],[ExpFinalSquidual/NP ExpFinalSquidual/NP])
            legend('Fit Squidua','Expected Final Squidua')
            hold on
            
            if (DoPlot>=2) % Plot the fitted noise and the expected noise
                t=1:1:NP;
                figure(10)
                plot(t,x,'b.',t,FitPos,'r-')
                

                if (DoPlot>=3)
                    
                    figure(12)
                    plot(t,Metric,'r-')
                    pause(1.00)
                end

                pause(0.50)
            end
        end

    else
        
        Metric(StepInd)=0;
        disp('a short-dwell step found and ignored')
        
    end
    
    if( FitSquidual<=ExpFinalSquidual && ExpNoise2==0 )
        xresidue = x-FitPos;
        [Spectrum, ff] = OneSidedSpectrumSS(xresidue, fSamp);
        tempsize = size(Spectrum);
        startingpoint = ceil(tempsize(1)*tempsize(2)/5);
        ExpFinalSquidual = mean(Spectrum(startingpoint:end)) *fSamp/2 *NP;
        disp(strcat('Current number of steps found = ',num2str(NoFittedSteps)))
        disp(strcat('Re-estimated noise2 = ',num2str(ExpFinalSquidual/NP)))
        
    end
    

end


%print while loop exit condition
disp('termination condition:')
if NoFittedSteps>=MaxNoSteps    
    disp('All steps found')
end   
if FitSquidual<=ExpFinalSquidual
    disp('Squared Residual reached the Expected Squared Noise')
    disp(strcat('Number of steps found = ',num2str(NoFittedSteps)))
end
if ( (-tempDeltaQ)<=MinDeltaQ ) 
    disp('Squidual reduced per step smaller than the threshold')
end

FinalSquidual=sum((x-FitPos-mean(x-FitPos)).^2)/NP;

% %plot residual
% figure
% plot(x-FitPos)
% legend('residuals')

%find current error on step position
MeanSquares=RCumx_sq(StepsInd(2:end))./R_N(StepsInd(2:end));
SquareMean=RCumx(StepsInd(2:end)).*RCumx(StepsInd(2:end))./R_N(StepsInd(2:end))./R_N(StepsInd(2:end));
s=(MeanSquares-SquareMean)./R_N(StepsInd(2:end));

%find error on step transition
% s_d=sqrt(s(2:end).*s(2:end)./R_N(StepsInd(3:end))+s(1:end-1).*s(1:end-1)./R_N(StepsInd(2:end-1)));

%calculate Step Statistics
StepStatistics.StepsChi2=-RCumx(StepsInd(2:end-1)).*RCumx(StepsInd(2:end-1))./R_N(StepsInd(2:end-1))-LCumx(StepsInd(2:end-1)+1).*LCumx(StepsInd(2:end-1)+1)./L_N(StepsInd(2:end-1)+1)+(RCumx(StepsInd(2:end-1))+LCumx(StepsInd(2:end-1)+1)).*(RCumx(StepsInd(2:end-1))+LCumx(StepsInd(2:end-1)+1))./(R_N(StepsInd(2:end-1))+L_N(StepsInd(2:end-1)+1));
StepStatistics.NumberOfStepsFound=length(StepsInd)-2;
StepStatistics.StepSizeStats=FitPos(StepsInd(2:end-1)+1)-FitPos(StepsInd(2:end-1));
StepStatistics.StepLengths=StepsInd(2:end)-StepsInd(1:end-1);
StepStatistics.StepFit=FitPos;
% StepStatistics.StepSizeError=s_d;

%%
function newMetric = UpdateMetric(RCumx,LCumx,R_N,L_N)

newMetric = -RCumx(1:end-1).*RCumx(1:end-1)./R_N(1:end-1)...
            -LCumx(2:end  ).*LCumx(2:end  )./L_N(2:end  )...
            +(RCumx(1:end-1)+LCumx(2:end  )).*(RCumx(1:end-1)+LCumx(2:end  ))./(R_N(1:end-1)+L_N(2:end  ));


%%
% Variable List:
% MaxNoSteps
% MinPointsInStep
% x
% ExpNoise
% DoPlot
% 
% 
% 
% IfCorrect           replace it later
% NP
% 
% x
% x_sq
% 
% FitPos
% 
% RCumx_sq
% LCumx_sq
% 
% RCumx
% LCumx
% 
% R_N
% L_N
% 
% Metric
% 
% 
% StepLoc             position in StepsInd
% 
% StepsInd
% StepInd
% PrevInd
% NextInd
% Prev2Ind
% Next2Ind
% 
% 
% FitSquidual
% ExpNoise       ExpAveSquidual
% ExpFinalSquidual           ExpectedFinalChi2
% 
% 
% NoFittedSteps
% MaxNoSteps
% 
% 
% 
% temperary variables
% 
% tempMinVal
% StepLocs


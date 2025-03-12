    function [bestParams, bestLoss] = TestToFindStepThresholdRemoval_op(matFile)
    %OPTIMIZESTEPPARAMS Optimizes step detection parameters based on target data.
    % Auther, Joseph Slivka, Leo Li. With Use of AI
    %   Loads an nx6 array from a .mat file, extracts measured positions (col 1)
    %   and fitted positions (col 3), then finds the best step detection parameters
    %   to minimize the standard deviation of position errors.
    %
    %   INPUTS:
    %     handles  : Struct containing GUI or experimental data handles
    %     matFile  : String, filename of the .mat file containing nx6 data array
    %
    %   OUTPUTS:
    %     bestParams : [best_step_thresh, best_avg_window, best_minsteplength]
    %     bestLoss   : Minimum standard deviation of error
    
        % Load the .mat file and extract relevant columns
    
    if nargin < 1
        [matFile, matPath] = uigetfile('*.mat');
        dataStruct = load(fullfile(matPath,matFile));
    else
        dataStruct = load(matFile);
    end

    rawData = dataStruct.data.trace; % Assume the first field is the data array
    x_measured = rawData(:,1);% Column 1: Measured positions
    x_fitted_truth  = rawData(:,3); % Column 3: Fitted positions
    x_fitted_unprocessed = SICstepFinder_op(x_measured').StepFit;
     

    
    if all(isnan(x_fitted_unprocessed)) || isempty(x_fitted_unprocessed)
        error('SICstepFinder_op returned invalid steps - check input data');
    end
    % Define search ranges for parameters
    step_thresh_values = 4:1:10;  
    avg_window_values = 18:1:28;
    minsteplength_values = 5:1:15;
        
    % step_thresh_values = 1:1:40;  
    % avg_window_values = 1:1:40;
    % minsteplength_values = 1:1:40;
    
    bestLoss   = inf;
    bestParams = [NaN, NaN, NaN];
    bestfit = [];
    bestN = 0;

    if isempty(step_thresh_values) || isempty(avg_window_values) || isempty(minsteplength_values)
       error('Parameter ranges must not be empty');
    end


    if all(isnan(x_fitted_unprocessed)) || numel(unique(x_fitted_unprocessed)) < 2
        error('Initial step detection failed - check SICstepFinder_op');
    end

    lossMatrix = zeros(length(step_thresh_values), length(avg_window_values), length(minsteplength_values));

    % Grid search loop with robust loss calculation
    for stIdx = 1:length(step_thresh_values)
        st = step_thresh_values(stIdx);
        for awIdx = 1:length(avg_window_values)
            aw = avg_window_values(awIdx);
            for mlIdx = 1:length(minsteplength_values)
                ml = minsteplength_values(mlIdx);
                x_fitted_truth_temp = horzcat(nan(1,aw),x_fitted_truth',nan(1,aw));
                x_fitted_midprocessed = StepThresholdRemoval_op(x_measured', x_fitted_unprocessed, [], st, aw, ml);
                if all(isnan(x_fitted_midprocessed))
                    fprintf('Invalid params: st=%.1f aw=%d ml=%d\n', st, aw, ml);
                    continue;
                end

                if length(x_fitted_truth_temp) == length(x_fitted_midprocessed)
                    valid_mask = ~isnan(x_fitted_truth_temp) & ~isnan(x_fitted_midprocessed);
                    
                    if sum(valid_mask) < 10  % Require min 10 valid points
                        currLoss = inf;
                    else

                        dx = x_fitted_truth_temp(valid_mask) - x_fitted_midprocessed(valid_mask);
                        N = sum(abs(diff(x_fitted_midprocessed(valid_mask))) > 0); % Number of steps
                        NP = sum(valid_mask); % Number of valid points

                        % BIC -  Bayesian Information Criterion (BIC) -
                        % Gives 4,20,7 as best answers


                        % if NP > 0 && N > 0
                        %     logLikelihood = -NP/2 * log(2*pi*var(dx)) - 1/(2*var(dx)) * sum(dx.^2);
                        %     currLoss = -logLikelihood + (N+1)/2 * log(NP);
                        % else
                        %     currLoss = inf;
                        % end


                        % Schwarz Information Criterion (SIC)
                        sigma2 = var(dx); % Variance of residuals
                        currLoss = N * log(NP) + NP * log(sigma2); % SIC formula
                    end
                else
                    currLoss = inf;
                end

                lossMatrix(stIdx, awIdx, mlIdx) = currLoss;

                % Update best parameters
                if currLoss < bestLoss
                    bestLoss = currLoss;
                    bestParams = [st, aw, ml];
                    bestfit = x_fitted_midprocessed(aw:end-aw);
                    bestN = sum(abs(x_fitted_midprocessed(2:end)- x_fitted_midprocessed(1:end-1)) > 0);
                end

                

            end
        end
    end


% Display results
fprintf('Best params: stepThresh=%.2f, avgWindow=%d, minStepLength=%d, Loss=%.4f\n, N=%.4f\n',...
bestParams(1), bestParams(2), bestParams(3), bestLoss, bestN);


figure;
plot(x_measured, 'k'); hold on;
plot(x_fitted_truth, 'b'); hold on;
plot(x_fitted_unprocessed, 'g'); hold on;
plot(bestfit, 'r');
title(sprintf('st=%.1f aw=%d ml=%d', bestParams(1), bestParams(2), bestParams(3)));
legend('Truth','Processed');


plotParameterSpace(lossMatrix, step_thresh_values, avg_window_values, minsteplength_values, bestParams);

    end

    function plotParameterSpace(lossMatrix, st_vals, aw_vals, ml_vals, bestParams)
    % Extract best parameters
    best_st = bestParams(1);
    best_aw = bestParams(2);
    best_ml = bestParams(3);
    
    % Ensure best parameters are within the defined ranges
    if ~ismember(best_st, st_vals) || ~ismember(best_aw, aw_vals) || ~ismember(best_ml, ml_vals)
        error('Best parameters are outside the defined ranges.');
    end
    
    % Create 3D plots while holding one parameter fixed at its best value
    figure;
    
    % 1. Fix step_thresh at best value, vary avg_window and minsteplength
    subplot(1,3,1);
    st_idx = find(st_vals == best_st); % Index of best step_thresh
    Z1 = squeeze(lossMatrix(st_idx, :, :)); % Slice for fixed step_thresh
    [X1, Y1] = meshgrid(aw_vals, ml_vals); % Create grid for avg_window and minsteplength
    surf(X1, Y1, Z1');
    xlabel('Avg Window');
    ylabel('Min Step Length');
    zlabel('SIC Loss');
    title(sprintf('Step Thresh = %.1f (fixed)', best_st));
    colorbar;
    
    % 2. Fix avg_window at best value, vary step_thresh and minsteplength
    subplot(1,3,2);
    aw_idx = find(aw_vals == best_aw); % Index of best avg_window
    Z2 = squeeze(lossMatrix(:, aw_idx, :)); % Slice for fixed avg_window
    [X2, Y2] = meshgrid(st_vals, ml_vals); % Create grid for step_thresh and minsteplength
    surf(X2, Y2, Z2');
    xlabel('Step Thresh');
    ylabel('Min Step Length');
    zlabel('SIC Loss');
    title(sprintf('Avg Window = %d (fixed)', best_aw));
    colorbar;
    
    % 3. Fix minsteplength at best value, vary step_thresh and avg_window
    subplot(1,3,3);
    ml_idx = find(ml_vals == best_ml); % Index of best minsteplength
    Z3 = squeeze(lossMatrix(:, :, ml_idx)); % Slice for fixed minsteplength
    [X3, Y3] = meshgrid(st_vals, aw_vals); % Create grid for step_thresh and avg_window
    surf(X3, Y3, Z3');
    xlabel('Step Thresh');
    ylabel('Avg Window');
    zlabel('SIC Loss');
    title(sprintf('Min Step Length = %d (fixed)', best_ml));
    colorbar;
    end



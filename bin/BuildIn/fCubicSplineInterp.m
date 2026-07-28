function delta = fCubicSplineInterp(x, y, vals, getfunc)
% JS 2022/09/03
% Description:
%    helper function for fCalcAlignCorrection to interpolate x,y data, using MATLAB resources:
%     https://www.mathworks.com/help/matlab/ref/boundary.html
%     https://www.mathworks.com/help/matlab/ref/scatteredinterpolant.html
%     https://www.mathworks.com/help/matlab/ref/griddata.html
% Input:
%    x,y,vals: array to trace out 3d to interpolate
%    getfunc (optional): if nonempty will skip testing and rather just
%    return the function
% Output:
%    output will be either the tested delta mean error (for x,y,vals)
%    OR will be the interpolation function to apply to correct for offsets


if length(y) ~= length(x) 
    error("sizes must be the same for interpolation")
end

% if no fourth argument to get function, just do a test
if nargin < 4
    train_fraction = 0.95; %train with 95% of the points
    num_tests = 20; %how many different train/test sets to get an average
    
    % interpolation only works within the interior. So let's make sure that
    % we include the boundary in the training data
    % Documentation: https://www.mathworks.com/help/matlab/ref/boundary.html
    must_train = unique(boundary(x,y))';
    sample = setdiff(1:length(vals), must_train);
    
    %FitError = ones(1,num_tests);
    FitErrorF = ones(1,num_tests);
    for test = 1:num_tests
    
        % and split according to how much training fraction one wants
        if length(sample) < round(length(vals)*train_fraction)-length(must_train)
            error("Too little non-border points sampled. Decrease train_fraction")
        end
        trainidx = randsample(length(sample), round(length(vals)*train_fraction)-length(must_train));
        trainidx = sort([sample(trainidx), must_train]);
        testidx = setdiff(1:length(vals),trainidx);

        % Now do the interpolation
        %delta = griddata(x(trainidx), y(trainidx), vals(trainidx), x(testidx), y(testidx), 'cubic');
        %Estimate Error
        %FitError(test) = mean(abs(delta(~isnan(delta)) - vals(testidx(~isnan(delta)))));
        
        % Or get a function F(x,y) using scatteredInterpolant. It seems like C1
        % derivative continuity (this is quadratic I believe) gives a
        % better result than the cubic griddata interpolation
        % Documentation: https://www.mathworks.com/help/matlab/ref/scatteredinterpolant.html
        F = scatteredInterpolant(x(trainidx), y(trainidx), vals(trainidx), 'natural');
        FitErrorF(test) = mean(abs(F(x(testidx), y(testidx)) - vals(testidx)));
    end
    
    %AverageFitError = mean(FitError);
    AverageFitErrorF = mean(FitErrorF);
    
    delta = mean(AverageFitErrorF);
    
    % plot for visualization
%     figure()
%     legend()
%     hold on
%     title("WARNING: Data Points Outside Interpolation Region were Found!")
%     scatter(x(trainidx), y(trainidx), 'b', 'DisplayName', 'Random Chosen Train Points')
%     scatter(x(must_train), y(must_train), 'g', 'DisplayName', 'Border (must train) Points')
%     scatter(x(testidx), y(testidx), 'k',  'DisplayName', 'Test Points')
%     scatter(x(testidx(isnan(delta))), y(testidx(isnan(delta))), 'r', 'DisplayName', 'Error Points')
    
else
    
    %delta = griddata(x, y, vals, xreturn, yreturn, 'cubic');
    Fdelta = scatteredInterpolant(x, y, vals, 'natural');
    delta = Fdelta;
    
end

end


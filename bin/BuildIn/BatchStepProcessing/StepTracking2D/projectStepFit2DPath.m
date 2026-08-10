function [onPath, offPath, closestXY] = projectStepFit2DPath(xy, PathModel)
%PROJECTSTEPFIT2DPATH Project points to a linear/quadratic fitted path.
%   onPath is signed arc-length coordinate; offPath is signed normal offset.

    if size(xy,2) ~= 2
        error('projectStepFit2DPath:Input', 'XY must be N-by-2.');
    end
    xy = double(xy);
    origin = double(PathModel.Origin);
    axisU = double(PathModel.AxisU);
    axisV = double(PathModel.AxisV);
    coefficients = double(PathModel.Coefficients);
    u0 = (xy - origin) * axisU.';
    v0 = (xy - origin) * axisV.';
    uRange = double(PathModel.URange(:).');
    if numel(uRange) ~= 2 || any(~isfinite(uRange))
        error('projectStepFit2DPath:URange', ...
            'PathModel.URange must contain two finite coordinate bounds.');
    end
    uRange = sort(uRange);
    span = max(diff(uRange), 1);
    uRange = uRange + [-0.05 0.05] * span;
    derivative = polyder(coefficients);
    n = size(xy,1);
    u = zeros(n,1);
    offPath = zeros(n,1);

    if PathModel.Degree == 1
        slope = coefficients(1);
        intercept = coefficients(2);
        u = (u0 + slope * (v0 - intercept)) / (1 + slope^2);
        u = min(max(u, uRange(1)), uRange(2));
        v = polyval(coefficients, u);
        normal = [-slope, 1] / sqrt(1 + slope^2);
        offPath = (u0-u) * normal(1) + (v0-v) * normal(2);
    else
        for i = 1:n
            objective = @(q) (q-u0(i)).^2 + (polyval(coefficients,q)-v0(i)).^2;
            u(i) = fminbnd(objective, uRange(1), uRange(2));
            v = polyval(coefficients, u(i));
            normal = [-polyval(derivative,u(i)), 1];
            normal = normal / norm(normal);
            offPath(i) = dot([u0(i)-u(i), v0(i)-v], normal);
        end
    end

    % Arc length is tabulated once and interpolated, avoiding a numerical
    % integration for every localization.
    grid = linspace(uRange(1), uRange(2), 2001);
    speed = sqrt(1 + polyval(derivative,grid).^2);
    arc = cumtrapz(grid, speed);
    referenceArc = interp1(grid, arc, 0, 'linear', 'extrap');
    onPath = interp1(grid, arc, u, 'linear', 'extrap') - referenceArc;
    v = polyval(coefficients, u);
    closestXY = origin + u * axisU + v * axisV;
end

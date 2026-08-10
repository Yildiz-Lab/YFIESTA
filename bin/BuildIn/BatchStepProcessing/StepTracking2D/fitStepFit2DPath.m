function PathModel = fitStepFit2DPath(xy, degree)
%FITSTEPFIT2DPATH Fit a linear or quadratic centreline to a 2-D trajectory.
%   PathModel = fitStepFit2DPath(XY, DEGREE), where DEGREE is 1 or 2.
%   The path is fitted as v = p(u) in a PCA-aligned coordinate system.

    if nargin < 2, degree = 1; end
    if ~ismember(degree, [1 2])
        error('fitStepFit2DPath:Degree', 'Degree must be 1 (linear) or 2 (quadratic).');
    end
    if size(xy,2) ~= 2 || any(~isfinite(xy(:)))
        error('fitStepFit2DPath:Input', 'XY must be a finite N-by-2 array.');
    end
    % fminbnd, used by quadratic closest-point projection, requires double
    % bounds. MINFLUX Results can be stored as single precision, so make
    % the entire path model double at its source.
    xy = double(xy);

    origin = mean(xy, 1);
    [~,~,V] = svd(xy - origin, 'econ');
    axisU = V(:,1).';
    netMotion = xy(end,:) - xy(1,:);
    if norm(netMotion) > eps && dot(axisU, netMotion) < 0
        axisU = -axisU;
    end
    axisV = [-axisU(2), axisU(1)];
    uv = [(xy-origin) * axisU.', (xy-origin) * axisV.'];
    coefficients = polyfit(uv(:,1), uv(:,2), degree);

    PathModel.Type = sprintf('polynomial-degree-%d', degree);
    PathModel.Degree = degree;
    PathModel.Origin = origin;
    PathModel.AxisU = axisU;
    PathModel.AxisV = axisV;
    PathModel.Coefficients = coefficients;
    PathModel.URange = [min(uv(:,1)), max(uv(:,1))];
end

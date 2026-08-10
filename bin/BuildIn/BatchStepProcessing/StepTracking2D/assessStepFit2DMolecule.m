function [passed, Quality] = assessStepFit2DMolecule(Molecule, FilterOptions)
%ASSESSSTEPFIT2DMOLECULE Apply transparent, conservative pre-fit QC rules.

    passed = false;
    Quality.Passed = false;
    Quality.Reason = '';
    Quality.PointCount = 0;
    Quality.NetDisplacement = NaN;
    Quality.PathLength = NaN;
    Quality.Directionality = NaN;
    Quality.Elongation = NaN;
    if ~isfield(Molecule, 'Results') || size(Molecule.Results,2) < 4
        Quality.Reason = 'Missing Results(:,3:4) coordinates'; return
    end
    xy = Molecule.Results(:,3:4);
    Quality.PointCount = size(xy,1);
    if any(~isfinite(xy(:)))
        Quality.Reason = 'Non-finite x/y coordinate'; return
    end
    if Quality.PointCount < FilterOptions.MinPoints
        Quality.Reason = sprintf('Too short (%d < %d points)', Quality.PointCount, FilterOptions.MinPoints); return
    end
    increments = diff(xy,1,1);
    Quality.PathLength = sum(hypot(increments(:,1), increments(:,2)));
    Quality.NetDisplacement = norm(xy(end,:) - xy(1,:));
    Quality.Directionality = Quality.NetDisplacement / max(Quality.PathLength, eps);
    centred = xy - mean(xy,1);
    [~,singularValues] = svd(centred, 'econ');
    singularValues = diag(singularValues);
    if numel(singularValues) < 2
        Quality.Elongation = Inf;
    else
        Quality.Elongation = (singularValues(1)^2) / max(singularValues(2)^2, eps);
    end
    if Quality.NetDisplacement < FilterOptions.MinNetDisplacement
        Quality.Reason = 'Net displacement below threshold'; return
    end
    if Quality.Directionality < FilterOptions.MinDirectionality
        Quality.Reason = 'Low directionality (diffusive-looking trajectory)'; return
    end
    if Quality.Elongation < FilterOptions.MinElongation
        Quality.Reason = 'Low path elongation (diffusive-looking trajectory)'; return
    end
    Quality.Passed = true;
    Quality.Reason = 'Passed';
    passed = true;
end

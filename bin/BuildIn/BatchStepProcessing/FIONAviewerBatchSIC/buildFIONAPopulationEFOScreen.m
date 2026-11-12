function Screen = buildFIONAPopulationEFOScreen(medians, multiplier, gradeCap, reference)
%BUILDFIONAPOPULATIONEFOSCREEN Compare one median EFO per eligible molecule.
%   A supplied reference can come from a manually checked population. The
%   automatic median needs at least five finite, positive molecule medians.
    if ~isempty(medians)
        validateattributes(medians, {'numeric'}, {'real','vector'});
    end
    validateattributes(multiplier, {'numeric'}, {'real','scalar','finite','>',1});
    gradeCap = validatestring(gradeCap, {'off','A','B'});
    if ~isempty(reference)
        validateattributes(reference, {'numeric'}, {'real','scalar','finite','positive'});
    end
    medians = double(medians(:));
    valid = isfinite(medians) & medians > 0;
    medians(~valid) = NaN;
    nReference = nnz(valid);
    source = 'disabled';
    threshold = NaN;
    ref = NaN;
    if ~strcmp(gradeCap,'off')
        if ~isempty(reference)
            ref = double(reference);
            source = 'supplied';
        elseif nReference >= 5
            ref = median(medians(valid));
            source = 'eligible batch median';
        else
            source = 'insufficient eligible traces';
        end
        threshold = ref * multiplier;
    end
    Screen = struct('MedianEFO',medians, 'ReferenceEFO',ref, ...
        'ReferenceTraceCount',nReference, 'ReferenceSource',source, ...
        'Multiplier',double(multiplier), 'ThresholdEFO',threshold, ...
        'GradeCap',gradeCap, 'Ratio',medians/ref, 'Flagged',valid & (medians>threshold));
end

function Result = mergeStepFit2DOptions(Default, Override)
%MERGESTEPFIT2DOPTIONS Recursively apply a partial user options structure.
    Result = Default;
    if nargin < 2 || isempty(Override), return, end
    fields = fieldnames(Override);
    for i = 1:numel(fields)
        name = fields{i};
        if isfield(Result, name) && isstruct(Result.(name)) && isstruct(Override.(name))
            Result.(name) = mergeStepFit2DOptions(Result.(name), Override.(name));
        else
            Result.(name) = Override.(name);
        end
    end
end

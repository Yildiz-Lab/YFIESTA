function requireFreshFIONAEFOOutput(outputFolder)
%REQUIREFRESHFIONAEFOOUTPUT Prevent stale S/A files on a screened refit.
    files = dir(fullfile(outputFolder,'**','*_fiona.mat'));
    if any(~startsWith({files.name},'.'))
        error('batchFIONAviewerSIC:ScreenedOutputExists', ...
            ['Population EFO screening requires a fresh output folder for a new fit. ' ...
             'Use another OutputFolder, or ResumeExisting with the same population and screening settings: %s'], ...
            outputFolder);
    end
end

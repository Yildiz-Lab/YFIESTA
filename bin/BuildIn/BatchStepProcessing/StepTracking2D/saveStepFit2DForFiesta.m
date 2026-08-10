function saveStepFit2DForFiesta(filename, Molecule, Filament)
%SAVESTEPFIT2DFORFIESTA Save processed tracks in FIESTA LoadTracks format.
%   The saved MAT-file contains variables named exactly Molecule and
%   Filament, matching fMenuData > SaveTracks and LoadTracks.  Pipeline-only
%   fields are removed from this GUI copy because LoadTracks concatenates
%   incoming structs with FIESTA's standard structure definition.

    moleculeFields = {'Name','File','Comments','Results','Type','Selected', ...
        'Visible','Drift','PixelSize','Channel','Color','PathData', ...
        'PlotHandles','TrackingResults'};
    filamentFields = {'Name','File','Comments','Results','Selected','Visible', ...
        'Drift','PixelSize','Channel','Color','PosStart','PosCenter','PosEnd', ...
        'Data','PathData','PlotHandles','TrackingResults'};

    % fDefStructure preserves unrecognised fields.  FIESTA subsequently
    % horizontally concatenates the loaded tracks with its own standard
    % structures, which requires identical field names *and field order*.
    Molecule = retainFiestaFields(Molecule, moleculeFields, 'Molecule');
    if nargin < 3 || isempty(Filament)
        % Construct a genuinely empty Filament structure with FIESTA's
        % complete schema and canonical order; a Data-only struct is not
        % sufficient because it produces a different field order on load.
        Filament = emptyFiestaStruct(filamentFields);
    else
        Filament = retainFiestaFields(Filament, filamentFields, 'Filament');
    end
    [folder,~,extension] = fileparts(filename);
    if isempty(extension), filename = [filename '.mat']; end
    if ~isempty(folder) && ~isfolder(folder)
        error('saveStepFit2DForFiesta:Folder', 'Output folder does not exist: %s', folder);
    end
    save(filename, 'Molecule', 'Filament', '-v6');
end

function Object = retainFiestaFields(Object, requiredFields, objectName)
    if ~isstruct(Object)
        error(['saveStepFit2DForFiesta:' objectName 'Type'], '%s must be a structure array.', objectName);
    end
    present = fieldnames(Object);
    missing = setdiff(requiredFields, present);
    if ~isempty(missing)
        error(['saveStepFit2DForFiesta:' objectName 'Fields'], ...
            '%s is missing required FIESTA fields: %s', objectName, strjoin(missing, ', '));
    end
    extras = setdiff(present, requiredFields);
    if ~isempty(extras)
        Object = rmfield(Object, extras);
    end
    Object = orderfields(Object, requiredFields);
    % LoadTracks concatenates along columns, so use FIESTA's conventional
    % 1-by-N orientation even if an upstream pipeline made a column array.
    Object = reshape(Object, 1, []);
end

function Object = emptyFiestaStruct(requiredFields)
    Object = cell2struct(cell(1, numel(requiredFields)), requiredFields, 2);
    Object = repmat(Object, 1, 0); % Explicitly a 1-by-0 FIESTA array.
end

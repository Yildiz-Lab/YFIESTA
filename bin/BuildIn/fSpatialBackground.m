function varargout=fSpatialBackground(action,varargin)
%FSPATIALBACKGROUND Spatial-background settings, thresholding and pattern checks.
% Actions: defaults/get/set/enabled/select/params/load/reference/gaussian/validate/check/threshold.
% NormalizedMap is a fixed pattern: 1 = median signal above camera baseline.
switch action
    case 'defaults'
        varargout{1}=struct('Enabled',false,'File','','Baseline',100,'Strength',1, ...
            'NormalizedMap',[],'SourceMedian',[],'FrameCount',0,'SourceFile','','BlurSigma',[]);
    case 'get'
        threshold=varargin{1}; channel=varargin{2};
        s=fSpatialBackground('defaults');
        if isfield(threshold,'SpatialBackground')
            entries=threshold.SpatialBackground;
            if numel(entries)>=channel
                if iscell(entries), entry=entries{channel}; else, entry=entries(channel); end
                if isstruct(entry)
                    fields=fieldnames(s);
                    for k=1:numel(fields)
                        if isfield(entry,fields{k}), s.(fields{k})=entry.(fields{k}); end
                    end
                end
            end
        end
        varargout{1}=s;
    case 'set'
        threshold=varargin{1}; channel=varargin{2};
        entries={};
        if isfield(threshold,'SpatialBackground')
            entries=threshold.SpatialBackground;
            if ~iscell(entries), entries=num2cell(entries); end
        end
        if numel(entries)<channel
            entries(end+1:channel)=repmat({fSpatialBackground('defaults')},1,channel-numel(entries));
        end
        entries{channel}=varargin{3};
        threshold.SpatialBackground=entries;
        varargout{1}=threshold;
    case 'enabled'
        varargout{1}=false;
        for channel=varargin{2}
            s=fSpatialBackground('get',varargin{1},channel);
            varargout{1}=varargout{1} || s.Enabled;
        end
    case 'select'
        entries=cell(1,numel(varargin{2}));
        for k=1:numel(entries)
            entries{k}=fSpatialBackground('get',varargin{1},varargin{2}(k));
        end
        varargout{1}=entries;
    case 'params'
        p=varargin{1}; threshold=varargin{2};
        count=1;
        if iscell(threshold.Filter), count=numel(threshold.Filter); end
        p.spatial_background=fSpatialBackground('select',threshold,1:count);
        varargout{1}=p;
    case 'load'
        file=varargin{1}; baseline=varargin{2}; strength=varargin{3};
        validateattributes(baseline,{'numeric'},{'real','scalar','finite','nonnegative'});
        validateattributes(strength,{'numeric'},{'real','scalar','finite','nonnegative'});
        [image,count]=ReadReference(file);
        varargout{1}=NormalizeReference(image,file,baseline,strength,count);
        if nargout>1, varargout{2}=image; end
    case 'reference'
        [varargout{1:nargout}]=ReadReference(varargin{:});
    case 'gaussian'
        [varargout{1:nargout}]=GenerateGaussian(varargin{:});
    case 'validate'
        s=varargin{1}; imageSize=varargin{2}; mode=varargin{3};
        if ~s.Enabled, return; end
        if strcmp(mode,'variable')
            error('FIESTA:SpatialBackgroundMode','Spatial background supports Constant or Relative Intensity. Select one of these modes.');
        end
        validateattributes(s.Baseline,{'numeric'},{'real','scalar','finite','nonnegative'});
        validateattributes(s.Strength,{'numeric'},{'real','scalar','finite','nonnegative'});
        validateattributes(s.NormalizedMap,{'numeric'},{'real','2d','finite','nonempty','nonnegative'});
        if ~isempty(imageSize) && ~isequal(size(s.NormalizedMap),double(imageSize))
            error('FIESTA:BackgroundSize','Background map size %s does not match movie size %s. Load a matching map.', ...
                mat2str(size(s.NormalizedMap)),mat2str(imageSize));
        end
    case 'check'
        varargout{1}=CheckPattern(varargin{:});
    case 'threshold'
        [varargout{1:nargout}]=ThresholdImage(varargin{:});
    otherwise
        error('FIESTA:SpatialBackgroundAction','Unknown background action: %s',action);
end

function [image,count]=ReadReference(file)
% A stack contributes one fixed background: its mean over TIFF pages.
pages=imfinfo(file); count=numel(pages); image=[];
for k=1:count
    frame=imread(file,k,'Info',pages);
    validateattributes(frame,{'numeric'},{'real','2d','finite','nonempty'});
    if isempty(image), image=zeros(size(frame)); end
    if ~isequal(size(image),size(frame))
        error('FIESTA:BackgroundSize','All background TIFF pages must have the same size.');
    end
    image=image+double(frame)/count;
end

function s=NormalizeReference(image,file,baseline,strength,count)
signal=max(image-baseline,0); reference=median(signal(:));
if reference<=0
    error('FIESTA:BackgroundBaseline','The background image median must exceed the camera baseline.');
end
s=fSpatialBackground('defaults');
s.Enabled=true; s.File=char(file); s.Baseline=baseline; s.Strength=strength;
s.NormalizedMap=signal/reference; s.SourceMedian=median(image(:)); s.FrameCount=count;
s.SourceFile=char(file);

function [s,reference,blurred]=GenerateGaussian(file,sigma,baseline,strength,imageSize)
% Save floating-point counts without intensity rescaling or integer rounding.
validateattributes(sigma,{'numeric'},{'real','scalar','finite','positive'});
validateattributes(baseline,{'numeric'},{'real','scalar','finite','nonnegative'});
validateattributes(strength,{'numeric'},{'real','scalar','finite','nonnegative'});
[reference,count]=ReadReference(file);
blurred=single(imgaussfilt(reference,sigma,'Padding','symmetric'));
[folder,name]=fileparts(file); stem=sprintf('%s_gaussianblur_%.15g',name,sigma);
output=fullfile(folder,[stem '.tif']); version=2;
while isfile(output)
    output=fullfile(folder,sprintf('%s_%d.tif',stem,version)); version=version+1;
end
s=NormalizeReference(double(blurred),output,baseline,strength,1);
if nargin<5, imageSize=[]; end
fSpatialBackground('validate',s,imageSize,'constant');
tags=struct('ImageLength',size(blurred,1),'ImageWidth',size(blurred,2), ...
    'Photometric',Tiff.Photometric.MinIsBlack,'BitsPerSample',32,'SamplesPerPixel',1, ...
    'SampleFormat',Tiff.SampleFormat.IEEEFP,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky, ...
    'Compression',Tiff.Compression.LZW,'RowsPerStrip',min(64,size(blurred,1)), ...
    'Software','YFIESTA','ImageDescription',sprintf('Gaussian sigma %.15g pixels; mean of %d reference TIFF page(s). Source: %s',sigma,count,file));
tiff=Tiff(output,'w'); cleanup=onCleanup(@()close(tiff));
setTag(tiff,tags); write(tiff,blurred); clear cleanup;
s.SourceFile=char(file); s.BlurSigma=sigma;

function report=CheckPattern(image,settings)
%FSPATIALBACKGROUNDCHECK Advisory check of low-frequency pattern agreement.
% Compare shapes after normalizing brightness, not raw reference intensities.
report=struct('Warn',false,'Correlation',NaN,'ResidualFraction',0,'Message','');
if ~settings.Enabled || settings.Strength==0, return; end
fSpatialBackground('validate',settings,size(image),'constant');
image=double(image);amplitude=median(image(:))-settings.Baseline;
if amplitude<=0
    report.Warn=true;
    report.Message='Movie median is at or below the camera baseline; spatial correction has zero amplitude. Check the baseline setting.';
    return;
end
sigma=min(30,max(1,min(size(image))/16));
observed=imgaussfilt(image,sigma,'Padding','symmetric');
observed=(observed-median(observed(:)))/amplitude;
expected=imgaussfilt(double(settings.NormalizedMap)-1,sigma,'Padding','symmetric');
expected=expected-median(expected(:));
a=observed(:);b=expected(:);
if std(a)>eps && std(b)>eps
    correlation=corrcoef(a,b);report.Correlation=correlation(1,2);
end
report.ResidualFraction=sqrt(mean((a-b).^2));
variation=sqrt(mean(b.^2));
report.Warn=report.ResidualFraction>.15 || ...
    (variation>.05 && isfinite(report.Correlation) && report.Correlation<.5);
if report.Warn
    correlationText='not defined for a flat pattern';
    if isfinite(report.Correlation), correlationText=sprintf('%.2f',report.Correlation); end
    report.Message=sprintf(['The reference background pattern may differ from this movie. ' ...
        'Pattern correlation: %s; residual variation: %.0f%% of median signal above baseline. ' ...
        'Preview the mask and consider another map or a lower correction strength. This is advisory; analysis can continue.'], ...
        correlationText,100*report.ResidualFraction);
end

function [output,details]=ThresholdImage(input,params)
%FSPATIALBACKGROUNDTHRESHOLD Shared pixelwise correction for preview/analysis.
% Raw pixels and fitting images are unchanged. A baseline is a camera offset,
% not a noise standard deviation. Normalize reference signal above it to 1.
s=params.spatial_background;
if iscell(s), s=s{1}; end
mode='constant';
if ~isfield(params,'threshold') || isempty(params.threshold), mode='variable'; end
fSpatialBackground('validate',s,size(input),mode);
input=double(input);
frameMedian=median(input(:));
amplitude=max(frameMedian-s.Baseline,0);
offset=s.Strength*amplitude.*(double(s.NormalizedMap)-1);
corrected=input-offset;
threshold=params.threshold;
if ~isreal(threshold)
    values=sort(corrected(:)); count=numel(values);
    middle=values(max(1,round(count/2)));
    noise=.7413*(values(max(1,round(3*count/4)))-values(max(1,round(count/4))));
    threshold=round(middle+imag(threshold)/100*noise);
end
comparison=corrected;
if strcmp(params.binary_image_processing,'average')
    comparison=conv2(comparison,fspecial('average',3),'same');
end
if threshold>0
    output=comparison>threshold;
elseif isfield(params,'bw_region')
    output=params.bw_region;
else
    output=comparison>0;
end
if strcmp(params.binary_image_processing,'smooth')
    output=bwmorph(output,'close'); output=bwmorph(output,'open');
end
if nargout>1
    details=struct('FrameMedian',frameMedian,'Amplitude',amplitude, ...
        'Offset',offset,'CorrectedImage',corrected,'Threshold',threshold, ...
        'ThresholdMap',threshold+offset);
end

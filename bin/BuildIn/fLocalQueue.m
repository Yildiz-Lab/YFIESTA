function varargout=fLocalQueue(action,varargin)
%FLOCALQUEUE File-backed jobs; a sole matching loaded movie can stay in memory.
switch lower(action)
    case 'snapshot',varargout{1}=Snapshot(varargin{:});
    case 'append',varargout{1}=Append(varargin{:});
    case 'load',[varargout{1:nargout}]=LoadJob(varargin{:});
    case 'job',[varargout{1:nargout}]=RunJob(varargin{:});
    case 'run',[varargout{1:nargout}]=Run(varargin{:});
end

function job=Snapshot(job,config,objects,batch)
[~,job.LocalQueueId]=fileparts(tempname);
job.StackName=One(job.StackName);job.Directory=One(job.Directory);
if ~startsWith(job.Directory,filesep)&&isempty(regexp(job.Directory,'^[A-Za-z]:','once'))
    job.Directory=fullfile(pwd,job.Directory);
end
job.Directory=[strip(job.Directory,'right',filesep) filesep];
job.LocalSourceFiles={fullfile(job.Directory,job.StackName)};job.LocalSourceTimes=job.Time;
job.LocalReconnectObjects=[];
if job.FirstTFrame==0
    job.LocalReconnectObjects=objects;job.LocalSourceFiles={};
elseif ~batch&&~isempty(job.Drift)&&iscell(config.StackName)&&numel(unique(config.StackName))>1
    job.LocalSourceFiles=cell(size(config.StackName));job.LocalSourceTimes=config.Time;
    for k=1:numel(config.StackName)
        folder=config.Directory{min(k,numel(config.Directory))};job.LocalSourceFiles{k}=fullfile(folder,config.StackName{k});
    end
end

function value=One(value)
if iscell(value),value=value{1};end
value=char(value);

function queue=Append(queue,job)
% Older saved queues may lack the new metadata fields.
queue=queue(:)';
if ~isempty(job)
    if isempty(queue),queue=job;else
        fields=union(fieldnames(queue),fieldnames(job));
        for k=1:numel(fields)
            key=fields{k};if ~isfield(queue,key),[queue.(key)]=deal([]);end
            if ~isfield(job,key),job.(key)=[];end
        end
        queue=[queue orderfields(job,queue(1))];
    end
end
for k=1:numel(queue)
    if ~isfield(queue,'LocalQueueId')||isempty(queue(k).LocalQueueId)|| ...
            (k>1&&any(strcmp({queue(1:k-1).LocalQueueId},queue(k).LocalQueueId)))
        [~,queue(k).LocalQueueId]=fileparts(tempname);
    end
end

function [stack,time,job]=LoadJob(job,varargin)
if job.FirstTFrame==0,stack={};time={};return;end
if isfield(job,'LocalSourceFiles')&&~isempty(job.LocalSourceFiles),files=job.LocalSourceFiles;
else,files={fullfile(One(job.Directory),One(job.StackName))};end
if ~isempty(varargin)&&numel(files)~=1,error('FIESTA:QueueSource','Loaded-movie reuse requires one source file.');end
stack={};time={};
for source=1:numel(files)
    if ~isempty(varargin)
        s=varargin{1};t=varargin{2};
    else
        file=files{source};if ~isfile(file),error('FIESTA:QueueSource','Queued source is unavailable: %s',file);end
        [~,~,ext]=fileparts(file);
        if strcmpi(ext,'.mat')
            data=load(file,'Stack','TimeInfo');s=data.Stack;t=data.TimeInfo;
        elseif ismember(lower(ext),{'.nd2','.nd','.zvi','.dv'})
            [s,t]=fReadND2(file,job.StackReadOptions);
        else,[s,t]=fStackRead(file,job.StackReadOptions);
        end
    end
    if ~iscell(s),s={s};end
    if ~iscell(t),t={t};end
    for c=1:numel(s)
        n=size(s{c},3);planeTime=job.Time;
        if isfield(job,'LocalSourceTimes')&&~isempty(job.LocalSourceTimes),planeTime=job.LocalSourceTimes;end
        if isempty(planeTime),planeTime=NaN;else,planeTime=planeTime(min(source,numel(planeTime)));end
        valid=numel(t)>=c&&numel(t{c})==n&&all(isfinite(t{c}))&&all(diff(t{c})>0);
        if isfinite(planeTime)&&planeTime>0
            t{c}=(0:n-1)*planeTime;
        elseif ~valid
            if ~isfinite(planeTime)||planeTime==0
                if n>1,error('FIESTA:QueueTime','Queued movie needs a nonzero frame interval in ms.');end
                planeTime=0;
            end
            t{c}=(0:n-1)*abs(planeTime);
        end
    end
    stack=[stack s];time=[time t]; %#ok<AGROW>
    clear s t data;
end
if isempty(job.Drift)
    channel=job.Channel;if numel(stack)==1,channel=1;end
else,channel=1:numel(job.Drift);
end
if any(channel<1|channel>numel(stack)),error('FIESTA:QueueChannel','Queued channel is absent from the source movie.');end
stack=stack(channel);time=time(channel);
if job.FirstTFrame>size(stack{1},3)||min(job.LastFrame,size(stack{1},3))<job.FirstTFrame
    error('FIESTA:QueueFrames','Queued frame range is outside the source movie.');
end
if ~isempty(job.Region)&&~isequal(size(job.Region),[size(stack{1},1) size(stack{1},2)])
    error('FIESTA:QueueRegion','Queued region mask does not match this movie. Configure this movie separately.');
end
job.StackName=One(job.StackName);job.Directory=[strip(One(job.Directory),'right',filesep) filesep];

function [abort,file]=RunJob(job,varargin)
% File-loaded arrays are released on return; a reused GUI stack stays loaded.
[stack,time,job]=LoadJob(job,varargin{:});objects=[];
if job.FirstTFrame==0
    if ~isfield(job,'LocalReconnectObjects')||isempty(job.LocalReconnectObjects)
        error('FIESTA:QueueReconnect','Re-add this reconnect job to snapshot its source objects.');
    end
    objects=job.LocalReconnectObjects;
end
if isfield(job,'LocalReconnectObjects'),job=rmfield(job,'LocalReconnectObjects');end
[abort,file]=fAnalyseStack(stack,time,job,0,objects);

function [files,problem]=Run(main,unload)
global Queue Stack TimeInfo
files={};problem=[];if isempty(Queue),return;end
if getappdata(0,'LocalQueueRunning'),error('FIESTA:QueueBusy','The local queue is already running.');end
setappdata(0,'LocalQueueRunning',true);cleanup=onCleanup(@Finish);
Queue=Append(Queue,[]);if nargin<2,unload=true;end
try
    fRightPanel('UpdateQueue','Local');
    while ~isempty(Queue)
        job=Queue(1);fprintf('Local queue: %s (%d job(s) remaining)\n',job.StackName,numel(Queue));
        if numel(Queue)==1&&CanReuseLoaded(job)
            fprintf('Local queue: using the loaded movie without closing or reloading it.\n');
            [abort,file]=RunJob(job,Stack,TimeInfo);
        else
            if unload&&~isempty(Stack)
                main=getappdata(0,'hMainGui');
                fMenuData('CloseStack',main);main=getappdata(0,'hMainGui');
                if ~isempty(Stack),error('FIESTA:QueueUnload','Close the current movie before running the local queue.');end
            end
            [abort,file]=RunJob(job);
        end
        if abort
            error('FIESTA:QueueStopped','Analysis stopped for %s. This job and later jobs remain queued. Partial results: %s',job.StackName,file);
        end
        files{end+1}=file; %#ok<AGROW>
        % Never remove a different job if the user edited the list mid-analysis.
        if ~isempty(Queue)&&strcmp(Queue(1).LocalQueueId,job.LocalQueueId),Queue(1)=[];end
        fRightPanel('UpdateQueue','Local');drawnow;
    end
catch err
    problem=err;
end
setappdata(main.fig,'LocalQueueResults',files);setappdata(main.fig,'LocalQueueError',problem);

function reuse=CanReuseLoaded(job)
global Stack Config
reuse=false;
% Keep complex channel combinations/reconnects on the existing file-backed path.
if numel(Stack)~=1||isempty(Stack{1})||job.FirstTFrame==0||~isempty(job.Drift),return;end
if ~all(isfield(Config,{'StackName','Directory','StackReadOptions'}))|| ...
        (iscell(Config.StackName)&&numel(Config.StackName)~=1)|| ...
        (iscell(Config.Directory)&&numel(Config.Directory)~=1),return;end
if isempty(Config.StackName)||isempty(Config.Directory)|| ...
        ~isequaln(job.StackReadOptions,Config.StackReadOptions),return;end
if isfield(job,'LocalSourceFiles')&&~isempty(job.LocalSourceFiles),files=job.LocalSourceFiles;
else,files={fullfile(One(job.Directory),One(job.StackName))};end
reuse=numel(files)==1&&strcmp(files{1},fullfile(One(Config.Directory),One(Config.StackName)));

function Finish
setappdata(0,'LocalQueueRunning',false);
try,fRightPanel('UpdateQueue','Local');catch,end

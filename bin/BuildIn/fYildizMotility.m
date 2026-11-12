function varargout=fYildizMotility(action,varargin)
%FYILDIZMOTILITY Region statistics and export without changing GUI selections.
switch lower(action)
    case 'defaults'
        varargout{1}=struct('dthresh',1,'fthresh',2,'edgeregion',2, ...
            'velocity_moverstoggle',true,'runlength_moverstoggle',true, ...
            'landingrate_moverstoggle',false,'runlength_countedges',true, ...
            'pool',true,'regions','all','exportKymographs',true, ...
            'kymoDirection',1,'kymoOverlay',true);
    case 'regions',varargout{1}=RegionNumbers(varargin{:});
    case 'analyze',varargout{1}=Analyze(varargin{:});
    case 'plot',Plot(varargin{:});
    case 'export',Export(varargin{:});
    case 'savefigures',SaveFigures(varargin{:});
    case 'kymograph',varargout{1}=Kymograph(varargin{:});
    case 'savekymographs',varargout{1}=SaveKymographs(varargin{:});
end

function ids=RegionNumbers(text,regions)
valid=arrayfun(@(r)~isempty(r.ScanData)&&isfield(r.ScanData,'X')&&numel(r.ScanData.X)>1,regions);
if strcmpi(strtrim(text),'all'),ids=find(valid);else
    ids=[];tokens=regexp(strtrim(text),'[,;\s]+','split');
    for k=1:numel(tokens)
        if isempty(regexp(tokens{k},'^\d+([:-]\d+)?$','once'))
            error('FIESTA:MotilityRegions','Use all, 1:5, or 1,3,5 for regions.');
        end
        pair=regexp(tokens{k},'\d+','match');
        first=str2double(pair{1});last=first;if numel(pair)>1,last=str2double(pair{2});end
        if first<1||last<first||last>numel(regions),error('FIESTA:MotilityRegions','Region range is outside the loaded regions.');end
        ids=[ids first:last]; %#ok<AGROW>
    end
    ids=unique(ids);
    if any(~valid(ids)),error('FIESTA:MotilityRegions','Each requested region must contain a saved scan centerline.');end
end
if isempty(ids),error('FIESTA:MotilityRegions','No scan regions are available.');end

function report=Analyze(molecules,regions,context,o)
pix=context.PixelSizeNm;validateattributes(pix,{'numeric'},{'scalar','positive','finite'});
for key={'dthresh','fthresh','edgeregion'}
    validateattributes(o.(key{1}),{'numeric'},{'scalar','nonnegative','finite'});
end
frames=context.FrameRange;validateattributes(frames,{'numeric'},{'vector','numel',2,'integer','positive','finite'});
if frames(2)<frames(1),error('FIESTA:MotilityTime','The last kymograph frame precedes the first.');end
dt=context.FrameIntervalMs/1000;
if isfield(context,'TimeMs')&&numel(context.TimeMs)>=frames(2)
    times=double(context.TimeMs(:))/1000;steps=diff(times);steps=steps(steps>0&isfinite(steps));
    if ~isempty(steps),dt=median(steps);end
    duration=times(frames(2))-times(frames(1))+dt;
else,duration=diff(frames)*dt+dt;
end
validateattributes(duration,{'numeric'},{'scalar','positive','finite'});
if o.pool
    ids=RegionNumbers(o.regions,regions);scans=cell(1,numel(ids));
    for k=1:numel(ids)
        data=regions(ids(k)).ScanData;scans{k}=fScanGeometry(data.X,data.Y,data.ScanSize);
    end
else
    ids=0;scans={context.Scan};
    if isempty(context.Scan)||~isfield(context.Scan,'InterpD'),error('FIESTA:MotilityScan','Select a scan region first.');end
    for k=1:numel(regions)
        data=regions(k).ScanData;
        if ~isempty(data)&&isequal(data.X(:),context.Scan.X(:))&&isequal(data.Y(:),context.Scan.Y(:)),ids=k;break;end
    end
end
% Same >66% membership rule as a region click; overlap goes to the best match.
owner=zeros(1,numel(molecules));fraction=zeros(size(owner));
if o.pool
    for j=1:numel(molecules)
        m=molecules(j);if ~m.Visible||m.Channel~=context.Channel||isempty(m.Results),continue;end
        xy=double(m.Results(:,3:4))/pix;
        for k=1:numel(ids)
            scan=scans{k};inside=inpolygon(xy(:,1),xy(:,2),[scan.lx fliplr(scan.ux)],[scan.ly fliplr(scan.uy)]);
            share=mean(inside);
            if share>.66&&share>fraction(j),owner(j)=k;fraction(j)=share;end
        end
    end
else,owner([molecules.Selected]==1)=1;
end
rows=cell(0,19);summary=zeros(numel(ids),8);
for k=1:numel(ids)
    scan=scans{k};lengthNm=double(scan.InterpD(end))*pix;
    validateattributes(lengthNm,{'numeric'},{'scalar','positive','finite'});
    indices=find(owner==k);firstRow=size(rows,1)+1;positions=cell(1,numel(indices));
    for n=1:numel(indices)
        j=indices(n);m=molecules(j);rr=double(m.Results);
        keep=rr(:,1)>=frames(1)&rr(:,1)<=frames(2)&all(isfinite(rr(:,[1:4 6])),2);rr=rr(keep,:);
        if isempty(rr),continue;end
        distance=max(rr(:,6))-min(rr(:,6));
        if ~isempty(m.PathData)&&size(m.PathData,1)==numel(keep)
            path=double(m.PathData(keep,4));if all(isfinite(path)),distance=max(path)-min(path);end
        end
        elapsed=rr(end,2)-rr(1,2);lifetime=rr(end,1)-rr(1,1)>=o.fthresh&&elapsed>0;
        mover=distance>=o.dthresh*pix;
        xy=rr(:,3:4)/pix;pos=Project(xy,scan)*pix;net=pos(end)-pos(1);
        inside=inpolygon(xy(:,1),xy(:,2),[scan.lx fliplr(scan.ux)],[scan.ly fliplr(scan.uy)]);
        positions{n}=[pos inside];
        velocity=NaN;run=NaN;landing=NaN;interaction=NaN;
        if lifetime&&(~o.velocity_moverstoggle||mover),velocity=distance/elapsed;end
        if lifetime&&(~o.runlength_moverstoggle||mover),run=distance;interaction=elapsed;end
        if lifetime&&(~o.landingrate_moverstoggle||mover),landing=1;end
        id=regexp(m.Name,'\d+\s*$','match','once');id=str2double(id);
        rows(end+1,:)={id,velocity,run,landing,interaction,ids(k),j,m.Name,m.Channel,m.File, ...
            rr(1,1),rr(end,1),rr(1,2),rr(end,2),lifetime,mover,false,net,n}; %#ok<AGROW>
    end
    range=firstRow:size(rows,1);direction=0;
    if ~isempty(range)
        eligible=cell2mat(rows(range,15))&cell2mat(rows(range,16));
        signed=cell2mat(rows(range,18));direction=sign(sum(sign(signed(eligible))));
        for r=range
            pos=positions{rows{r,19}};
            atEnd=(direction>0&&any(pos(:,2)&pos(:,1)>=lengthNm-o.edgeregion*pix))|| ...
                (direction<0&&any(pos(:,2)&pos(:,1)<=o.edgeregion*pix));
            excluded=o.runlength_countedges&&direction~=0&&atEnd;
            rows{r,17}=excluded;if excluded,rows(r,[2 3 5])={NaN,NaN,NaN};end
        end
    end
    landingCount=sum(cell2mat(rows(range,4)),'omitnan');rate=landingCount/(lengthNm/1000*duration);
    for r=range,if isfinite(rows{r,4}),rows{r,4}=rate;end,end
    summary(k,:)=[ids(k) lengthNm/1000 duration landingCount rate direction sum(cell2mat(rows(range,17))) numel(range)];
end
names={'Molecule','Velocity_nm_s','RunLength_nm','LandingRate_um_s','InteractionTime_s', ...
    'Region','TrajectoryIndex','Name','Channel','SourceFile','StartFrame','EndFrame', ...
    'StartTime_s','EndTime_s','PassLifetime','IsMover','EndExcluded','NetAlongScan_nm'};
traces=cell2table(rows(:,1:18),'VariableNames',names);
for j=setdiff(1:18,[8 10]),traces.(names{j})=cell2mat(rows(:,j));end
regionsTable=array2table(summary,'VariableNames',{'Region','Length_um','Observation_s','LandingCount', ...
    'LandingRate_um_s','DominantDirection','EndExcluded','TrajectoryCount'});
report=struct('Traces',traces,'Regions',regionsTable,'Options',o,'Context',context, ...
    'PooledLandingRate',sum(summary(:,4))/sum(summary(:,2).*summary(:,3)));

function position=Project(xy,scan)
center=[mean(scan.InterpX,1)' mean(scan.InterpY,1)'];d=double(scan.InterpD(:));
best=inf(size(xy,1),1);position=zeros(size(best));
for k=1:size(center,1)-1
    v=center(k+1,:)-center(k,:);v2=sum(v.^2);if v2==0,continue;end
    t=max(0,min(1,(xy-center(k,:))*v'/v2));distance=sum((xy-(center(k,:)+t.*v)).^2,2);
    use=distance<best;best(use)=distance(use);position(use)=d(k)+t(use)*(d(k+1)-d(k));
end

function Plot(report,parent,region)
traces=report.Traces;stats=report.Regions;
if nargin>2,traces=traces(traces.Region==region,:);stats=stats(stats.Region==region,:);end
data={traces.Velocity_nm_s,traces.RunLength_nm,stats.LandingRate_um_s,traces.InteractionTime_s};
titles={'Velocity (nm/s)','Run length (nm)','Landing rate (um^{-1} s^{-1})','Interaction time (s)'};
for k=1:4
    ax=subplot(2,2,k,'Parent',parent);values=data{k};values=sort(values(isfinite(values)));cla(ax);
    if isempty(values),text(ax,.5,.5,'No qualifying trajectories','Units','normalized','HorizontalAlignment','center');
    else
        [x,~,group]=unique(values);y=cumsum(accumarray(group,1))/numel(values);
        stairs(ax,[x(1);x],[0;y],'LineWidth',1.3);
        legend(ax,sprintf('Mean %.3g\nStd %.3g\nN %d',mean(values),std(values),numel(values)),'Location','best');
    end
    title(ax,titles{k});ylabel(ax,'Cumulative fraction');ylim(ax,[0 1]);grid(ax,'on');
end

function Export(report,file,append)
if nargin<3,append=false;end
WriteSheet(file,1,report.Traces,append);WriteSheet(file,'Regions',report.Regions,append);
if isfield(report,'Kymographs'),WriteSheet(file,'Kymographs',report.Kymographs,append);end
o=report.Options;keys=fieldnames(o);values=cellfun(@(key)string(o.(key)),keys,'UniformOutput',false);
settings=[{'Setting','Value'};[keys values]];
settings(end+1,:)={'Pooled landing rate (um^-1 s^-1)',report.PooledLandingRate};
settings(end+1,:)={'Kymograph frame range',mat2str(report.Context.FrameRange)};
settings(end+1,:)={'Pixel size (nm)',report.Context.PixelSizeNm};
settings(end+1,:)={'Channel',report.Context.Channel};
settings(end+1,:)={'Overlap rule','Largest in-region fraction; ties to lower region number'};
settings(end+1,:)={'Dominant direction','Majority sign of net scan displacement among qualifying movers; tie = no end exclusion'};
writecell(settings,file,'Sheet','Options','Range','A1','WriteMode','overwritesheet');

function WriteSheet(file,sheet,data,append)
headers=data.Properties.VariableNames;old=cell(0,numel(headers));
if append&&isfile(file)
    sheets=sheetnames(file);exists=isnumeric(sheet)||any(strcmp(sheets,sheet));
    if exists
        previous=readcell(file,'Sheet',sheet);
        if ~isempty(previous)
            % Numeric IDs discard repeated text headers and empty legacy rows.
            good=cellfun(@(x)isnumeric(x)&&isscalar(x)&&isfinite(x),previous(:,1));
            if size(previous,2)>=8&&isequal(previous{1,8},'Name')
                good=good|cellfun(@(x)ischar(x)&&~isempty(x)&&~strcmp(x,'Name'),previous(:,8));
            end
            previous=previous(good,:);
            old=cell(size(previous,1),numel(headers));
            columns=min(size(previous,2),numel(headers));old(:,1:columns)=previous(:,1:columns);
        end
    end
end

blank=cellfun(@(x)isa(x,'missing')||(isstring(x)&&isscalar(x)&&ismissing(x)),old);old(blank)={''};
writecell([headers;old;table2cell(data)],file,'Sheet',sheet,'Range','A1','WriteMode','overwritesheet');

function kymo=Kymograph(report,region,regions,molecules,stack,settings)
% Sample exactly the region's scan width; horizontal distance, vertical time.
c=report.Context;frames=c.FrameRange(1):c.FrameRange(2);
if isempty(stack)||frames(end)>size(stack,3)
    error('FIESTA:MotilityStack','Load the current-channel stack covering the requested kymograph frames.');
end
validateattributes(report.Options.kymoDirection,{'numeric'},{'scalar','real','finite'});
if ~ismember(report.Options.kymoDirection,[-1 1]),error('FIESTA:MotilityDirection','Choose +1 or -1 for the kymograph direction.');end
data=regions(region).ScanData;scan=fScanGeometry(data.X,data.Y,data.ScanSize);
native=report.Regions.DominantDirection(report.Regions.Region==region);
mirrored=native~=0&&native~=report.Options.kymoDirection;
distance=double(scan.InterpD(:))'*c.PixelSizeNm/1000;image=zeros(numel(frames),numel(distance));
for n=1:numel(frames)
    x=scan.InterpX;y=scan.InterpY;
    if settings.CorrectDrift&&~settings.StackCorrected&&~isempty(settings.Drift)
        t=inv(settings.Drift(:,:,min(frames(n),size(settings.Drift,3))));
        nx=x*t(1,1)+y*t(2,1)+t(3,1);y=x*t(1,2)+y*t(2,2)+t(3,2);x=nx;
    end
    if settings.Method==1
        samples=interp2(double(stack(:,:,frames(n))),x,y,'nearest');image(n,:)=max(samples,[],1,'omitnan');
    else
        samples=interp2(double(stack(:,:,frames(n))),x,y);image(n,:)=mean(samples,1,'omitnan');
    end
end
image(~isfinite(image))=0;if mirrored,image=fliplr(image);end
time=(frames-frames(1))*c.FrameIntervalMs/1000;
if numel(c.TimeMs)>=frames(end),time=(double(c.TimeMs(frames))-double(c.TimeMs(frames(1))))/1000;end
rows=report.Traces(report.Traces.Region==region,:);tracks=struct('Name',{},'Index',{},'XUm',{},'Row',{},'EndExcluded',{});
if report.Options.kymoOverlay
    for j=1:height(rows)
        m=molecules(rows.TrajectoryIndex(j));
        if ~isempty(settings.Drift)&&isfield(m,'Drift')
            m=fCorrectPos(m,settings.Drift,double(settings.CorrectDrift||settings.StackCorrected));
        end
        r=double(m.Results);xy=r(:,3:4)/c.PixelSizeNm;
        inside=inpolygon(xy(:,1),xy(:,2),[scan.lx fliplr(scan.ux)],[scan.ly fliplr(scan.uy)]);
        use=r(:,1)>=frames(1)&r(:,1)<=frames(end)&inside&all(isfinite(xy),2);
        x=Project(xy,scan)*c.PixelSizeNm/1000;if mirrored,x=distance(end)-x;end
        x(~use)=NaN;
        % An acquisition gap separates line segments rather than drawing a bridge.
        order=(1:numel(x))';gap=find(diff(r(:,1))>1)+.5;order=sort([order;gap]);
        xx=nan(size(order));yy=xx;valid=order==round(order);
        xx(valid)=x(order(valid));yy(valid)=r(order(valid),1)-frames(1)+1;
        tracks(end+1)=struct('Name',m.Name,'Index',rows.TrajectoryIndex(j),'XUm',xx,'Row',yy,'EndExcluded',rows.EndExcluded(j)); %#ok<AGROW>
    end
end
kymo=struct('Image',image,'DistanceUm',distance,'TimeS',time,'Frames',frames,'Tracks',tracks, ...
    'NativeDirection',native,'Mirrored',mirrored,'Region',region,'Limits',settings.Limits);

function info=SaveKymographs(report,regions,molecules,stack,settings,folder,base)
values=zeros(height(report.Regions),8);
for n=1:height(report.Regions)
    region=report.Regions.Region(n);k=Kymograph(report,region,regions,molecules,stack,settings);
    file=fullfile(folder,sprintf('%s_region_%03d_kymograph',base,region));
    image=uint16(k.Image);
    if report.Options.kymoOverlay
        file=[file '_overlay'];image=KymographOverlay(k); %#ok<AGROW>
    end
    imwrite(image,[file '.tif'],'Compression','lzw');
    values(n,:)=[region k.NativeDirection report.Options.kymoDirection k.Mirrored report.Options.kymoOverlay settings.Method k.Frames(1) k.Frames(end)];
end
info=array2table(values,'VariableNames',{'Region','NativeDirection','RequestedDirection','Mirrored','MoleculeOverlay','Method_1Max_2Mean','FirstFrame','LastFrame'});

function rgb=KymographOverlay(k)
% Rasterize directly at the TIFF's native resolution; no figure or extra files.
limits=k.Limits;
if numel(limits)~=2||any(~isfinite(limits))||limits(2)<=limits(1)
    limits=prctile(k.Image(:),[1 99.8]);if limits(2)<=limits(1),limits(2)=limits(1)+1;end
end
gray=max(0,min(1,(k.Image-limits(1))/diff(limits)));rgb=repmat(gray,1,1,3);
[h,w]=size(gray);colors=[0 .85 1;1 .5 0];corners=[0 0;1 0;0 1;1 1];
for category=0:1
    pixels=cell(1,numel(k.Tracks));
    for j=1:numel(k.Tracks)
        track=k.Tracks(j);if track.EndExcluded~=category,continue;end
        x=1+track.XUm/k.DistanceUm(end)*(w-1);y=track.Row;points=zeros(0,2);
        for i=find(isfinite(x)&isfinite(y))'
            if i<numel(x)&&isfinite(x(i+1))&&isfinite(y(i+1))
                steps=max(1,ceil(2*max(abs([x(i+1)-x(i) y(i+1)-y(i)]))));
                points=[points;linspace(x(i),x(i+1),steps+1)' linspace(y(i),y(i+1),steps+1)']; %#ok<AGROW>
            else,points(end+1,:)=[x(i) y(i)];end %#ok<AGROW>
        end
        pixels{j}=points;
    end
    points=vertcat(pixels{:});if isempty(points),continue;end
    low=floor(points);indices=[];coverage=[];
    for corner=1:4
        p=low+corners(corner,:);weight=prod(1-abs(points-p),2);
        valid=p(:,1)>=1&p(:,1)<=w&p(:,2)>=1&p(:,2)<=h&weight>0;
        indices=[indices;p(valid,2)+(p(valid,1)-1)*h];coverage=[coverage;weight(valid)]; %#ok<AGROW>
    end
    alpha=reshape(accumarray(indices,coverage,[h*w 1],@max,0),h,w);
    rgb=rgb+(reshape(colors(category+1,:),1,1,3)-rgb).*alpha;
end
rgb=uint8(round(rgb*255));

function SaveFigures(report,folder,base)
for k=1:height(report.Regions)+double(report.Options.pool)
    f=figure('Visible','off','Color','white','Position',[100 100 900 650],'MenuBar','none','ToolBar','none');
    cleanup=onCleanup(@()close(f));
    if k<=height(report.Regions)
        region=report.Regions.Region(k);Plot(report,f,region);label=sprintf('Region %d',region);suffix=sprintf('_region_%03d',region);
    else
        Plot(report,f);label='Pooled regions';suffix='_pooled';
    end
    sgtitle(f,label);savefig(f,fullfile(folder,[base suffix '.fig']));
    print(f,fullfile(folder,[base suffix '.png']),'-dpng','-r150');clear cleanup;
end

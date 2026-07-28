% Skeletonization and then Hough
% https://www.mathworks.com/help/images/ref/bwskel.html

% Joseph Slivka: 2022/10/20


function NewLines = ASSAMM(I)

%% Description:
%   From a z-projection, find the longest, brightest line segment regions

% Parameters:
%   Image (.tif) : a z-projection picture that correlates to the image stack

% Returns:
%   NewLines : struct with fields bx, by for x,y coordinates of segmented line
%   avg_MT_length : average length of found segmented line length
%       (end-to-end)

% Extras:
%   number of hough peaks to use (straighter MTs means more peaks can be
%   used)
peaks = 20;
%   number of best MTs to return
chk = 15;
%   degree of connection to combine tracks
degthresh = 5;
%   neighbors to check
neighbors_check = 1;
%   pixel intersection to combine tracks that overlap
overlapthresh = 7;
%   for extrapolation post line finding "density" to be considered
extrapthresh = 0.15;
%   for cumulative threshold to call MT end (for just using density by
%   extrapthresh, set cumulthresh to 1)
cumulthresh = 0.85;

I = imadjust(I);
% figure
% imshow(I);
% title("Z-Projection Image")

%% Do edge detection
BW = imbinarize(I);
% Rather than binarize, do a thresholding method
% thresh = 18000;
% BW = I; BW(I <= thresh) = 0; BW(I > thresh) = 1;
% BW = logical(BW);
% figure
% imshow(BW)
% title("Binary Image from imbinarize")

%% Now skeletonize
outBW = bwskel(BW);
%outBW = bwmorph(BW,'skel',Inf);
% out is a matrix the size of BW but with points.
% figure
% imshow(labeloverlay(I,outBW,'Transparency',0))
% title("Skeleton Fit over Binary Image")

%% Calculate Hough transform from skeleton image
[H,theta,rho] = hough(outBW);
% figure
% imshow(imadjust(rescale(H)),[],...
%        'XData',theta,...
%        'YData',rho,...
%        'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('\rho')
% axis on
% axis normal 
% hold on
% colormap(gca,hot)
% title("Hough Space representation of the skeleton binary fit")

%% Find Peaks
P = houghpeaks(H,peaks,'threshold',ceil(0.3*max(H(:)))); %Straighter movies can have a higher value for P

% Show rho vs theta in Hough transform space. It is at the intersections
% that line parameters will be found and drawn.
% x = theta(P(:,2));
% y = rho(P(:,1));
% plot(x,y,'s','color','black');

%% Get lines that are associated with these peaks and show them on image

% Can change the FillGap and MinLength parameters depending on your line
% sparsity from the skeleton image
% We do some manual post processing to FillGap so you can keep this smaller
lines = houghlines(outBW,theta,rho,P,'FillGap',15,'MinLength',8);

% % Plot the Hough Lines
% figure, imshow(I), hold on
% max_len = 0;
% linelens = zeros(length(lines),1);
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%    text(xy(2,1),xy(2,2),strcat('  ',num2str(k)), 'Color', 'white')
% end

%% Post-Sorting
% after this maybe we can post process by grabbing the z-project value
% points and getting the one with the highest average

% First find the best candidates for kymos

linemeans = zeros(length(lines),1);
linestds = zeros(length(lines),1);

for k = 1:length(lines)
% make an interpolation
dist = pdist2(lines(k).point1, lines(k).point2);
xlin = linspace(lines(k).point1(1),lines(k).point2(1),round(4*dist,0));
ylin = linspace(lines(k).point1(2),lines(k).point2(2),round(4*dist,0));

% find the unique points
pts = [round(xlin,0)', round(ylin,0)'];
pts = unique(pts(:,1:2),'rows');

% get the values
ptsI = zeros(length(pts),1);
for j = 1:size(pts,1)
    ptsI(j) = double(I(pts(j,2),pts(j,1)));
end
linestds(k) = var(ptsI);
linemeans(k) = sum(ptsI);
end

[~, sortedInds] = sort(linemeans,'descend');
toplines = sortedInds(1:chk);

% % now plot these "top contenders" in blue
% 
% for i = 1:length(toplines)
%     k = toplines(i);
%     xy = [lines(k).point1; lines(k).point2];
% %     if i == 7
%         plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','cyan');
% %     else
% %         plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
% %     end
% 
%     % Plot beginnings and ends of lines
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%    
% end
% title("Hough Lines in Image Space with ''best'' lines by greatest sum in cyan")

%% Extend Top Contenders Lines
% Then, of these best candidates, try to extend them by looking for nearby
% neighbors of similar slope. If they are within a slope tolerance, then
% combine them by another segmented line when we export to region

% repackage data into an array for easier manipulation later
AllLines = zeros(length(lines),6);

for i = 1:length(lines)
    AllLines(i,1:2) = lines(i).point1; AllLines(i,3:4) = lines(i).point2;
    AllLines(i,5) = lines(i).theta; AllLines(i,6) = lines(i).rho; 
end

NewLinesbx = zeros(length(toplines),2*length(lines));
NewLinesby = zeros(length(toplines),2*length(lines));

for tcl = 1:length(toplines) %top contender line
    
    k = toplines(tcl);
    kbg = k; ked = k; % for beginning and end differences
    % For calling a FIESTA Scan, give it the center points in an array
    bx = [AllLines(kbg,1), AllLines(ked,3)];
    by = [AllLines(kbg,2), AllLines(ked,4)];
    
    % to break the for loop over connecting to this line
    extend_beginpoint = 1;
    extend_endpoint = 1;
    already_connected = [k];
    
%     %% check for spatial overlaps
%     % make region to check for overlap
%     pts = interp_region(AllLines(k,1),AllLines(k,3),AllLines(k,2),AllLines(k,4),size(I),0);
%     for ll = 1:length(lines)
%         if ~ismember(ll,already_connected)
%             ptsprime = interp_region(AllLines(ll,1),AllLines(ll,3),AllLines(ll,2),AllLines(ll,4),size(I),0);
%             % check that a certain number of these points overlap
%             ptsintersect = intersect(pts,ptsprime,'rows');
%             if length(ptsintersect) > overlapthresh
%                 fprintf(strcat("Overlap detected for molecules ", num2str(k), " and ", num2str(ll), "\n"))
%                 fprintf(strcat(num2str(length(intersect(pts,ptsprime,'rows')))), "\n")
%                 % append the pts to the backbone in order depending on
%                 % which is first x-wise
%                 if bx(1) < AllLines(ll,1)
%                 bx(2:4) = [ptsintersect(1,1), ptsintersect(end,1), AllLines(ll,3)];
%                 by(2:4) = [ptsintersect(1,2), ptsintersect(end,2), AllLines(ll,4)];
%                 else
%                 % have to switch order here
%                 bx(1:4) = [AllLines(ll,1), ptsintersect(1,1), ptsintersect(end,1), AllLines(k,3)];
%                 by(1:4) = [AllLines(ll,2), ptsintersect(1,2), ptsintersect(end,2), AllLines(k,4)];
%                 end
%                 % Just to make sure everything is in the minimal line order
%                 U = unique([bx', by'],'rows','stable');
%                 [bx, by] = minimize_line(U(:,1)', U(:,2)');
%                 already_connected = [already_connected, ll];
%             end
%         end
%     end
    
    %% Now look to extend lines by connecting with others
    for ll = 1:length(lines)
        if extend_beginpoint
        % Look at first point of k compared to others end points
        M = sqrt((AllLines(kbg,1) - AllLines(:,3)).^2 + (AllLines(kbg,2) - AllLines(:,4)).^2);
        [~,idx] = sort(M);
        
        % Verify that the line hasn't already been added
        for neighbors = 1:neighbors_check
            nn = idx(neighbors);
            for mm = 1:length(already_connected)
                % Remove self minimum if it happens
                if ismember(nn,already_connected)
                    nn = idx(mm+neighbors);
                end
            end

            % Check that the rho, theta values are relatively similar
            if (abs(AllLines(kbg,5) - AllLines(nn,5)) < degthresh) && (abs(AllLines(kbg,6) - AllLines(nn,6)) < degthresh)
                % make region to look for skel
                
                pts = interp_region(AllLines(nn,3),AllLines(kbg,1),AllLines(nn,4),AllLines(kbg,2),size(I));

                % check that a certain number of these points are in the skeleton
                count = 0;
                for p = 1:length(pts)
                    count = count + double(outBW(pts(p,2),pts(p,1))); % You have to invert it for the image
                end
                % If there is a good density of counts. In this case, 1/9 since 1/3
                % plus three parallel lanes (width = 1)
                if count > M(nn)/4
                    %fprintf(strcat("MATCHED with k=",num2str(nn),"! \n"))
                    bx = [AllLines(nn,1), AllLines(nn,3), bx]; by = [AllLines(nn,2), AllLines(nn,4), by];
                    kbg = nn;
                    already_connected = [already_connected, kbg];
                else
                    extend_beginpoint = 0;
                end

            end
        end 
        end % of first point extension
        
        
        if extend_endpoint
        % Now same thing, look at last point of k compared to others end points
        M = sqrt((AllLines(ked,3) - AllLines(:,1)).^2 + (AllLines(ked,4) - AllLines(:,2)).^2);
        [~,idx] = sort(M);

        % Verify that the line hasn't already been added
        for neighbors = 1:neighbors_check
            nn = idx(neighbors);
            for mm = 1:length(already_connected)
                % Remove self minimum if it happens
                if ismember(nn,already_connected)
                    nn = idx(mm+neighbors);
                end
            end

            % Check that the rho, theta values are relatively similar
            if (abs(AllLines(ked,5) - AllLines(nn,5)) < degthresh) && (abs(AllLines(ked,6) - AllLines(nn,6)) < degthresh)
                % makre region to look for skel
                pts = interp_region(AllLines(ked,3),AllLines(nn,1),AllLines(ked,4),AllLines(nn,2),size(I));

                % check that a certain number of these points are in the skeleton
                count = 0;
                for p = 1:length(pts)
                    count = count + double(outBW(pts(p,2),pts(p,1))); % You have to invert it for the image
                end
                % If there is a good density of counts. In this case, 1/9 since 1/3
                % / three parallel lanes (width = 1)
                if count > M(nn)/4
                    %fprintf(strcat("MATCHED with k=",num2str(nn),"! \n"))
                    bx = [bx, AllLines(nn,1), AllLines(nn,3)]; by = [by, AllLines(nn,2), AllLines(nn,4)];
                    ked = nn;
                    already_connected = [already_connected, ked];
                else
                    extend_endpoint = 0;
                end
            end

        end
        end % of first point extension

        if ~extend_beginpoint && ~extend_endpoint
            break
        end

    end

    % Append this to the NewLines array
    NewLinesbx(tcl,1:length(bx)) = bx;
    NewLinesby(tcl,1:length(by)) = by;

end

% Find only unique ones (since contender lines could attack to others)
[~, xidx] = unique(NewLinesbx, 'rows'); 
[~, yidx] = unique(NewLinesby, 'rows'); 
idx = union(xidx, yidx); % it has to not be in both to be unique
NewLinesbx = NewLinesbx(idx,:);
NewLinesby = NewLinesby(idx,:);

% %% Idea: extend ends of lines along the slope. If there is a reasonable density of skel
% % within, then add it
% NewLinesbx(:,1:10)
% for nl = 1:size(NewLinesbx,1)
%     % FROM RIGHT EDGE
%     bx = NewLinesbx(nl,:); by = NewLinesby(nl,:);
%     bx = bx(bx > 0); by = by(by > 0); % only nonzero
%     bxend = bx(end-1:end); byend = by(end-1:end); % only the last two points to make line
%     pts = extrap_region(bxend(1),bxend(2),byend(1),byend(2),size(I));
%     % check that a certain number of these points are in the skeleton
%     ct = 0;
%     counts = zeros(1,length(pts));
%     [m,~] = find(isnan(pts)); pts(m,:) = [];
%     for p = 1:length(pts)
%         ct = ct + double(outBW(pts(p,2),pts(p,1)));
%         counts(p) = ct;
%     end
%     
%     nl
%     ct/length(counts)
%     % do some integration method to find best cut off
%     % if want to just use extrapthresh, then set 
%     if ct/length(counts) > extrapthresh
%         M = abs(counts - cumulthresh*max(counts));
%         [~,idx] = sort(M);
%         NewLinesbx(nl,length(bx)) = pts(idx(1),1);
%         NewLinesby(nl,length(by)) = pts(idx(1),2);
%     end
%     
%     % FROM LEFT EDGE
%     bx = bx(bx > 0); by = by(by > 0); % only nonzero
%     bxbegin = bx(1:2); bybegin = by(1:2); % only the last two points to make line
%     pts = extrap_region(bxbegin(2),bxbegin(1),bybegin(2),bybegin(1),size(I));
%     % check that a certain number of these points are in the skeleton
%     ct = 0;
%     counts = zeros(1,length(pts));
%     [m,~] = find(isnan(pts)); pts(m,:) = [];
%     for p = 1:length(pts)
%         ct = ct + double(outBW(pts(p,2),pts(p,1)));
%         counts(p) = ct;
%     end
%     
%     ct/length(counts)
%     % do some integration method to find best cut off
%     % if want to just use extrapthresh, then set 
%     if ct/length(counts) > extrapthresh
%         M = abs(counts - cumulthresh*max(counts));
%         [~,idx] = sort(M);
%         NewLinesbx(nl,1) = pts(idx(1),1);
%         NewLinesby(nl,1) = pts(idx(1),2);
%     end
%     
% end
% NewLinesbx(:,1:10)

%% Finally, package into a struct to actually be iterable without a bunch of
% extra nans
NewLines = struct('bx',{[]},'by',{[]});
avg_MT_length = 0;
for nl = 1:size(NewLinesbx,1)
    bx = NewLinesbx(nl,:); by = NewLinesby(nl,:);
    NewLines.bx{nl} = bx(bx>0); NewLines.by{nl} = by(by>0); 
    avg_MT_length = avg_MT_length + sqrt((NewLines.bx{nl}(1) - NewLines.bx{nl}(end))^2 + (NewLines.by{nl}(1) - NewLines.by{nl}(end))^2);
end
avg_MT_length = avg_MT_length/length(NewLines.bx)

% plot the final result
figure
imshow(labeloverlay(I,outBW,'Transparency',0))
hold on
for nl = 1:length(NewLines.bx)
    plot(NewLines.bx{nl}([1,end]),NewLines.by{nl}([1,end]),'x','LineWidth',2,'Color','cyan');
    plot(NewLines.bx{nl},NewLines.by{nl},'LineWidth',2,'Color','cyan');
end
title("Final Regions after connecting likely similar MTs")


function pts = interp_region(x1,x2,y1,y2,upperlimits,width)
    
if nargin < 6
    width = 1;
end

% make interpolation
d = sqrt((x1-x2)^2 + (y1-y2)^2);
xlin = linspace(x1,x2,round(4*d,0));
ylin = linspace(y1,y2,round(4*d,0));

% find the unique points
pts = [round(xlin,0)', round(ylin,0)'];
pts = unique(pts(:,1:2),'rows');
lp = length(pts);

%then make a region of certain width one just to get nearby weird
%rounding
for p = 1:lp
    for x = -width:width
        for y = -width:width
            pts = [pts; pts(p,1) + x, pts(p,2) + y];
        end
    end
end
pts = unique(pts(:,1:2),'rows');

% remove any possible points outside of plot region
[m1,~] = find(pts(:,1) < 1); pts(m1,:) = [];
[m1,~] = find(pts(:,2) < 1); pts(m1,:) = [];
[m1,~] = find(pts(:,1) > upperlimits(1)); pts(m1,:) = [];
[m1,~] = find(pts(:,2) > upperlimits(2)); pts(m1,:) = [];


function pts = extrap_region(x1,x2,y1,y2,upperlimits,length)

d = pdist2(x2-x1,y2-y1);
if nargin < 6
    length = 1.5*d;
end
t = 1:250;
tend = t(end) + round(t(end)*length/d);
%textrap = 251:tend;

cx = polyfit([t(1),t(end)],[x1,x2],1);   % linear parametric fit
cy = polyfit([t(1),t(end)],[y1,y2],1);   % linear parametric fit
Xextrap = polyval(cx,tend);
Yextrap = polyval(cy,tend);

pts = interp_region(x2,Xextrap,y2,Yextrap,upperlimits,1);


function [linex, liney] = minimize_line(x,y)

v = perms(1:length(x));
linex = zeros(1,length(x)); liney = zeros(1,length(x));
minD = 50000;
for j = 1:length(v)
    d = sqrt(sum((x(v(j,2:end))-x(v(j,1:end-1))).^2 + (y(v(j,2:end))-y(v(j,1:end-1))).^2));
    if d < minD
        linex = x(v(j,:));
        liney = y(v(j,:));
        minD = d;
    end 
end



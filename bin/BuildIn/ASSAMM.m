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
peaks = 15;
%   number of best MTs to return
chk = 12;
%   degree of connection to combine tracks
degthresh = 5;

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

    for ll = 1:length(lines)
        if extend_beginpoint
        % Look at first point of k compared to others end points
        M = sqrt((AllLines(kbg,1) - AllLines(:,3)).^2 + (AllLines(kbg,2) - AllLines(:,4)).^2);
        [~,idx] = sort(M);

        % Verify that the line hasn't already been added
        nn = idx(1);
        for mm = 1:length(already_connected)
            % Remove self minimum if it happens
            if ismember(nn,already_connected)
                nn = idx(mm+1);
            end
        end

        % Check that the rho, theta values are relatively similar
        if (abs(AllLines(kbg,5) - AllLines(nn,5)) < 5) && (abs(AllLines(kbg,6) - AllLines(nn,6)) < 5)
            % make interpolation
            xlin = linspace(AllLines(nn,3),AllLines(kbg,1),round(4*M(nn),0));
            ylin = linspace(AllLines(nn,4),AllLines(kbg,2),round(4*M(nn),0));

            % find the unique points
            pts = [round(xlin,0)', round(ylin,0)'];
            pts = unique(pts(:,1:2),'rows');
            lp = length(pts);

            %then make a region of certain width one just to get nearby weird
            %rounding
            width = 1;
            for p = 1:lp
                for x = -width:width
                    for y = -width:width
                        pts = [pts; pts(p,1) + x, pts(p,2) + y];
                    end
                end
            end
            pts = unique(pts(:,1:2),'rows');
            length(pts);
            % remove any possible points outside of plot region
            pts(pts(:,1) < 1) = []; pts(pts(:,2) < 1) = []; 
            pts(pts(:,1) > size(I,1)) = []; pts(pts(:,2) > size(I,2)) = [];


            % check that a certain number of these points are in the skeleton
            count = 0;
            for p = 1:length(pts)
                count = count + double(outBW(pts(p,2),pts(p,1))); % You have to invert it for the image
            end
            % If there is a good density of counts. In this case, 1/9 since 1/3
            % plus three parallel lanes (width = 1)
            if count > M(nn)/3
                %fprintf(strcat("MATCHED with k=",num2str(nn),"! \n"))
                bx = [AllLines(nn,1), AllLines(nn,3), bx]; by = [AllLines(nn,2), AllLines(nn,4), by];
                kbg = nn;
                already_connected = [already_connected, kbg];
            else
                extend_beginpoint = 0;
            end

        end
        end % of first point extension

        if extend_endpoint
        % Now same thing, look at last point of k compared to others end points
        M = sqrt((AllLines(ked,3) - AllLines(:,1)).^2 + (AllLines(ked,4) - AllLines(:,2)).^2);
        [~,idx] = sort(M);

        % Verify that the line hasn't already been added
        nn = idx(1);
        for mm = 1:length(already_connected)
            % Remove self minimum if it happens
            if ismember(nn,already_connected)
                nn = idx(mm+1);
            end
        end

        % Check that the rho, theta values are relatively similar
        if (abs(AllLines(ked,5) - AllLines(nn,5)) < degthresh) && (abs(AllLines(ked,6) - AllLines(nn,6)) < degthresh)
            % make interpolation
            xlin = linspace(AllLines(ked,3),AllLines(nn,1),round(4*M(nn),0));
            ylin = linspace(AllLines(ked,4),AllLines(nn,2),round(4*M(nn),0));

            % find the unique points
            pts = [round(xlin,0)', round(ylin,0)'];
            pts = unique(pts(:,1:2),'rows');
            lp = length(pts);

            %then make a region of certain width one just to get nearby weird
            %rounding
            width = 1;
            for p = 1:lp
                for x = -width:width
                    for y = -width:width
                        pts = [pts; pts(p,1) + x, pts(p,2) + y];
                    end
                end
            end
            pts = unique(pts(:,1:2),'rows');
            length(pts);
            % remove any possible points outside of plot region
            pts(pts(:,1) < 1) = []; pts(pts(:,2) < 1) = []; 
            pts(pts(:,1) > size(I,1)) = []; pts(pts(:,2) > size(I,2)) = [];


            % check that a certain number of these points are in the skeleton
            count = 0;
            for p = 1:length(pts)
                count = count + double(outBW(pts(p,2),pts(p,1))); % You have to invert it for the image
            end
            % If there is a good density of counts. In this case, 1/9 since 1/3
            % plus three parallel lanes (width = 1)
            if count > M(nn)/3
                %fprintf(strcat("MATCHED with k=",num2str(nn),"! \n"))
                bx = [bx, AllLines(nn,1), AllLines(nn,3)]; by = [by, AllLines(nn,2), AllLines(nn,4)];
                ked = nn;
                already_connected = [already_connected, ked];
            else
                extend_endpoint = 0;
            end
        end

        end

        if ~extend_beginpoint && ~extend_endpoint
            break
        end

    end
    
    %NewLines.bx{tcl} = bx;
    %NewLines.by{tcl} = by;
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

% Finally, package into a struct to actually be iterable without a bunch of
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

end

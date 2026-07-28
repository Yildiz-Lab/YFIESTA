

% Make sure to replace with your file location
load('C:\Users\yildi\Downloads\241109-092420_minflux_large_centered.mat')

% find the final iteration
idx1 = find(itr == 4);

% find channel indices
ch1 = find(thi == 0);
ch2 = find(thi == 1);

% find the intersect of final iteration and specific channel
ch1idx = intersect(idx1,ch1);
ch2idx = intersect(idx1,ch2);

% generate figure
figure()
hold on

scatter(loc(ch1idx,1),loc(ch1idx,2),10,'b','filled','DisplayName','640')
scatter(loc(ch2idx,1),loc(ch2idx,2),10,'g','filled','DisplayName','561')
xlabel('Position (nm)')
ylabel('Position (nm)')
legend()


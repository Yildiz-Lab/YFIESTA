limit = length(traces);

traces_SIC = [];
compareresults = NaN(i,5);
for (i=1:limit)
    if length(traces{i} > 0)
    traces_SIC{i} = parse_and_fit_twosides(traces{i}(:,1)',traces{i}(:,2)',traces{i}(:,6)',1,500000);
    dwells1 = add_to_list_6col_dwells(traces_SIC{i},1,[],0.002);
    dwells2 = add_to_list_6col_dwells(traces{i},1,[],0.002);
    steps1 = add_to_list_6col_steps(traces_SIC{i},1,[]);
    steps2 = add_to_list_6col_steps(traces{i},1,[]);
    if (length(dwells1) > 0)
    compareresults(i,:) = [i nanmean(dwells1(:,3)) nanmean(dwells2(:,3)) nanmean(steps1(:,1)) nanmean(steps2(:,1)) ];
    end
    end
end

limit = length(traces_SIC);

dwells_SIC = [];
steps_SIC = [];
for (i=1:limit)
    if length(traces_SIC{i}) > 0
    dwells_SIC = add_to_list_6col_dwells(traces_SIC{i},1,dwells_SIC,0.002);
    steps_SIC = add_to_list_6col_steps(traces_SIC{i},1,steps_SIC);
    end
end
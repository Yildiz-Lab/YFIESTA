% Fitting an MTIMBS curve with the hill equation

conc1 = [5; 10; 30; 100; 300];
INT_STD1 = [382, 134.8; 1078, 288.4; 2021, 731.2; 3959, 1725.8; 3634, 1672.5];
%INT_STD1 = [382, 134.8; 1078, 288.4; 2021, 731.2; 3959, 1725.8; 3959, 1672.5];

conc2 = [5; 10; 30; 100; 500; 1000];
INT_STD2 = [2071, 925.9; 2934, 1319; 5817, 2414.3; 7617, 2560; 6552, 1934.1; 6377, 1513.8];
%INT_STD2 = [2071, 925.9; 2934, 1319; 5817, 2414.3; 7617, 2560; 7617, 1934.1; 7617, 1513.8];

% MAP7 plot
figure()
% plot(conc1, INT_STD1(:,1)/max(INT_STD1(:,1)), 'o', 'DisplayName','MAP7')
errorbar(conc1, INT_STD1(:,1)/max(INT_STD1(:,1)), INT_STD1(:,2)/max(INT_STD1(:,1)), 'o', 'DisplayName','MAP7')
hold on

mdl_hill = fittype('x/(k+x)','indep','x');
mdl_hill_n = fittype('x^n/(k+x^n)','indep','x');
fittedmdl1 = fit(conc1,INT_STD1(:,1)/max(INT_STD1(:,1)),mdl_hill,'start',[10.],'weight',1./INT_STD1(:,2));
fittedmdl2 = fit(conc2,INT_STD2(:,1)/max(INT_STD2(:,1)),mdl_hill,'start',[10.],'weight',1./INT_STD2(:,2));

fittedmdln1 = fit(conc1,INT_STD1(:,1)/max(INT_STD1(:,1)),mdl_hill_n,'start',[10.,1.],'weight',1./INT_STD1(:,2));
fittedmdln2 = fit(conc2,INT_STD2(:,1)/max(INT_STD2(:,1)),mdl_hill_n,'start',[10.,1.],'weight',1./INT_STD2(:,2));

X = 0:max(conc1);
plot(X,fittedmdl1(X),'DisplayName', strcat("k = ",num2str(fittedmdl1.k), ", n = 1") )
plot(X,fittedmdln1(X),'DisplayName',strcat("k = ",num2str(fittedmdln1.k), ", n = ", num2str(fittedmdln1.n)) )
legend()


% Rigor Kinesin Plot
figure()
hold on
X = 0:max(conc2);
% plot(conc2, INT_STD2(:,1)/max(INT_STD2(:,1)), 'o', 'DisplayName','Rigor Kinesin')
errorbar(conc2, INT_STD2(:,1)/max(INT_STD2(:,1)), INT_STD2(:,2)/max(INT_STD2(:,1)), 'o', 'DisplayName','Rigor Kinesin')
plot(X,fittedmdl2(X),'DisplayName',strcat("k = ",num2str(fittedmdl2.k), ", n = 1") )
plot(X,fittedmdln2(X),'DisplayName',strcat("k = ",num2str(fittedmdln2.k), ", n = ", num2str(fittedmdln2.n)) )
legend()

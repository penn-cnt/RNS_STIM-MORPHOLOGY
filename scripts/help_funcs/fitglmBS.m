function tbl = fitglmBS(x,y,nperm,intercept)
output = submodel([x,y],intercept);
R2 = output(1);
b = output(2);
% bootstrapping
[ci,results] = bootci(nperm,@(x) submodel(x,intercept),[x,y]);
b_null = results(:,2);
pCoeff = [];
pv = mean(b_null > 0);
pv = min(pv,1-pv);
pCoeff = [pCoeff,2*pv];

% permutating
R2_null = [];
for i = 1:nperm
    inds = randperm(length(x));
    output = submodel([x,y(inds)],intercept);
    R2_null = [R2_null;output(1)];
end

pR2 = mean(R2_null > R2);

% stats

tbl = array2table(horzcat([R2,pR2;b,pCoeff],ci'),...
    'VariableNames',{'Estimate','pValue','lowerCI','upperCI'},'RowNames',{'R2','b'});

end

function output = submodel(data,intercept)
mdl = fitglm(data(:,1),data(:,2),'Intercept',intercept);
b = mdl.Coefficients{end,"Estimate"};
% p = mdl.Coefficients{:,"pValue"};
R2 = mdl.Rsquared.Ordinary;
output = [R2,b];
end

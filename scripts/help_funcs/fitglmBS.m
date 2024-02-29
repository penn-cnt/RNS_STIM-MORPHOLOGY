function tbl = fitglmBS(x,y,nperm)
output = submodel([x,y]);
R2 = output(1);
b = output(2);
% permutated
[ci,results] = bootci(nperm,@submodel,[x,y]);
R2_null = results(:,1);
b_null = results(:,2);
pR2 = mean(R2_null < 0);

% stats
pCoeff = [];
pv = mean(b_null > 0);
pv = min(pv,1-pv);
pCoeff = [pCoeff,2*pv];
tbl = array2table(horzcat([R2,pR2;b,pCoeff],ci'),...
    'VariableNames',{'Estimate','pValue','lowerCI','upperCI'},'RowNames',{'R2','b'});

end

function output = submodel(data)
mdl = fitglm(data(:,1),data(:,2),'Intercept',false);
b = mdl.Coefficients{:,"Estimate"};
% p = mdl.Coefficients{:,"pValue"};
R2 = mdl.Rsquared.Ordinary;
output = [R2,b];
end

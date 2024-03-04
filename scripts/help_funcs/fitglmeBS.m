function [results,intercept,outliers] = fitglmeBS(tbl,formula,nperm)

tmp = strsplit(formula,'~');
outcome = strip(tmp{1},' ');
o = tbl.(outcome);
naninds = isnan(o);
inds = abs(zscore(o(~naninds))) > 3;
outliers = find(~naninds);
outliers = outliers(inds);
tbl = tbl(~isnan(tbl.(outcome)),:);
tbl(abs(zscore(tbl.(outcome))) > 3,:) = [];
mdl = fitglme(tbl,formula);
b = mdl.Coefficients{end,"Estimate"};
R2 = mdl.Rsquared.Ordinary;
intercept = mdl.Coefficients{1,"Estimate"};

% bootstrapping
nsample = size(tbl,1);
seed = 0;
while true
    try
        rng(seed)
        [~,bootinds] = bootstrp(nperm,[],ones(nsample,1));
        results = [];
        parfor i = 1:nperm
                output = submodel(tbl(bootinds(:,i),:),formula);
                results = [results;output];
        end
        break
    catch
        seed = seed+1;
    end
end
R2_null = results(:,1);
b_null = results(:,2);
ci = [prctile(R2_null,2.5),prctile(R2_null,97.5);prctile(b_null,2.5),prctile(b_null,97.5)];
pR2 = mean(R2_null < R2);

% stats
pCoeff = [];
pv = mean(b_null > 0);
pv = min(pv,1-pv);
pCoeff = [pCoeff,2*pv];
results = array2table(horzcat([R2,pR2;b,pCoeff],ci),...
    'VariableNames',{'Estimate','pValue','lowerCI','upperCI'},'RowNames',{'R2','b'});
end

function output = submodel(tbl_sub,formula)
mdl = fitglme(tbl_sub,formula);
b = mdl.Coefficients{end,"Estimate"};
R2 = mdl.Rsquared.Ordinary;
output = [R2,b];
end
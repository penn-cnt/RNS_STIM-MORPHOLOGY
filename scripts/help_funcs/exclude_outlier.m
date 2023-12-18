function mask = exclude_outlier(data,thres)

% Calculate the z-score for each row
zScores = zscore(data);

% Identify and exclude rows with z-scores beyond the threshold
mask = ~any(abs(zScores) > thres, 2);

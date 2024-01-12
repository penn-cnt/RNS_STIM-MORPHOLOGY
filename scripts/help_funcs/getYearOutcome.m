function [field_data,outcome] = getYearOutcome(data,year,outcome_field,data_field)
% given years get outcome of certain year and data if outcome not available
if isnumeric(year)
    data = data(cellfun(@(x) length(x) >= year,{data.(outcome_field)}));
    outcome = cellfun(@(x) x(year), {data.(outcome_field)});
else
    outcome = cellfun(@(x) x(end), {data.(outcome_field)});
end
if ~iscell(data_field)
    field_data = {data.(data_field)};
else
    for i = 1:length(data_field)
        field_data{i} = {data.(data_field{i})};
    end
end
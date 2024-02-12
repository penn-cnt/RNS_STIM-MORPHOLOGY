function data = paddata(data,field,n_bin)
ndim = ndims(data(1).(field));
if iscell(data(1).(field))
    for n = 1:length(data)
        if ~isempty(data(n).(field))
            try
                if ndim == 2
                    data(n).(field) = cellfun(@(x) x(1:n_bin,:),data(n).(field),'UniformOutput',false);
                elseif ndim == 3
                    data(n).(field) = cellfun(@(x) x(1:n_bin,:,:),data(n).(field),'UniformOutput',false);
                end
            catch
                data(n).(field) = cellfun(@(x) padarray(x,n_bin-size(x,1),nan,'post'),data(n).(field),'UniformOutput',false);
            end
        end
    end
else
    for n = 1:length(data)
        if ~isempty(data(n).(field))
            try
                if ndim == 2
                    data(n).(field) = data(n).(field)(1:n_bin,:);
                elseif ndim == 3
                    data(n).(field) = data(n).(field)(1:n_bin,:,:);
                end

            catch
                data(n).(field) = padarray(data(n).(field),n_bin-size(data(n).(field),1),nan,'post');
            end
        end
    end
end
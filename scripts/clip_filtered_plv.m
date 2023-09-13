function [plv] = clip_filtered_plv(eeg, fs, freq_bands,options)
% LOOK AT OPTIMIZING USING https://iopscience.iop.org/article/10.1088/1741-2552/aacfe4
% LOOK AT CREATING METHOD FOR AVERAGING ACROSS TIME AND PERSERVING TRIAL
% RESOLUTION2
% INPUTS
% eeg           -   channels x time x trials
% fs            -   sampling frequency of eeg data (double hz)
% freq_bands    -   frequency bands to calculate within and cross phase
%                   locking (conditions x boundaries)
% options       -   structure containing the following fields:
%                   selection   -   data selection array for different
%                                   conditions
%                   basis       -   wavelet basis for the cwt

% OUTPUTS
% plv           - resulting plv matrix, time x channel x channel x cond
connections =  [1,2; % same
                3,4; % same
                1,3; % diff
                1,4; % diff
                2,3; % diff
                3,4];% diff
num_connections = size(connections,1);
num_freqs = size(freq_bands,1);
filters = cell(1,num_freqs);

for f = 1:length(filters)
     b = fir1(options.orders,freq_bands(f,:)/(fs/2));
     filters{f} = b;
end

% Now we have all of the filters we will need to use for calculating the
% plv's
% Next step is to iterate between the different frequency bands
plv = zeros(num_connections,num_freqs);
for conn = 1:num_connections
    for freq = 1:num_freqs
        signal1 = squeeze(eeg(:,connections(conn,1)));
        signal2 = squeeze(eeg(:,connections(conn,2)));
        % calculate the angle of each signal using the hilbert
        % transform and angle function
        fsignal1 = hilbert(filtfilt(filters{freq},1,signal1));
        fsignal2 = hilbert(filtfilt(filters{freq},1,signal2));
        % pass angles into equation for plv
        plv(conn,freq) = def_plv(fsignal1,fsignal2);
    end
end

    function val = def_plv(signal1,signal2)
        % standardize each signal
        signal1 = signal1./abs(signal1); signal2 = signal2./abs(signal2);
        % calculate output using PLV formula
        val = abs(conj(signal2')*signal1)/length(signal1);
    end

end


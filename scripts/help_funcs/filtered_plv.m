function plv = filtered_plv(eeg, fs, freq_bands,options)
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
                2,4];% diff
num_connections = size(connections,1);
num_freqs = size(freq_bands,1);
filters = {};
for f = 1:length(freq_bands)
     [b,a] = butter(options.orders,freq_bands(f,:)/(fs/2));
     filters{f,1} = b;
     filters{f,2} = a;
end
for i = 1:size(eeg,2)
    for freq = 1:length(filters)
        waves(freq,:,i) = filtfilt(filters{freq,1},filters{freq,2},eeg(:,i));
    end
end
plv = zeros(num_connections,num_freqs);
for conn = 1:num_connections
    for freq = 1:num_freqs
        % calculate the angle of each signal using the hilbert
        % transform and angle function
        sig1 = hilbert(waves(freq,:,connections(conn,1)));
        sig2 = hilbert(waves(freq,:,connections(conn,2)));
        plv(conn,freq) = def_plv(sig1,sig2);
    end
end
% plv = [freq_bands,plv'];

function val = def_plv(signal1,signal2)
    % standardize each signal
    signal1 = signal1./abs(signal1); signal2 = signal2./abs(signal2);
    % calculate output using PLV formula
    val = abs(signal1*signal2')/length(signal1);
end
end


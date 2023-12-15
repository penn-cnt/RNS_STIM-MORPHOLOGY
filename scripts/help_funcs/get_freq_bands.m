function freq_bands = get_freq_bands(freq_range,fs)
data = rand(1000,1);
[~,f] = cwt(data,'amor',fs,'FrequencyLimits',freq_range,'VoicesPerOctave',6);
f = f(end:-1:1);
freq_bands = [[freq_range(1);f(1:end-1)],f(1:end)];
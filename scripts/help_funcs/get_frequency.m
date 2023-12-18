function f = get_frequency(siglength,freq_range,fs)
fb = cwtfilterbank('SignalLength',siglength,'Wavelet','amor',...
'SamplingFrequency',fs,'FrequencyLimits',freq_range,'VoicesPerOctave',6);
f = centerFrequencies(fb);

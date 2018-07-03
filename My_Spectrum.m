function [ freq, yfp] = My_Spectrum( wave, fs, ref, fc, span)
%  [ freq, yfp] = My_Spectrum( wave, fs, ref, fc, span)
%   Param:
%       wave: Input signal in time domain
%       fs:   Sampling frequency
%       ref:  Referency Mag.
%       fc:   Center frequency
%       span: Frequency span
%   Return:
%       freq: Frequency vector
%       yfp:  Signal in frequency domain
    N = length(wave);
    w = hann(N);%hamming hann
    wave_w = wave.* w;
    yf = fft(wave_w)/sum(w);
    yf = fftshift(yf);
    freq = linspace(-fs/2,fs/2, N);
    yfp = 20*log10(abs(yf)/ref);
    if fc+span/2>fs/2
        error('Center Frequency or frequency high should be small than half Samping Frequency');
    end
    Nhigh = floor(N/2 +  N*(fc+span/2)/fs);
    Nlow = floor(N/2 +  N*(fc-span/2)/fs);
    freq = freq(Nlow+1:Nhigh);
    yfp = yfp(Nlow+1:Nhigh);
end


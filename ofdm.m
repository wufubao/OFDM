clear;
clc;
close all;
%% Config
M = 64;                 % 16-QAM 
k = log2(M);            % Number of bits per symbol
R = 9216;               % Total bits length (Can be devide by NFFT)
OverSamp = 4;           % Upsampling factor
L = 16;                 % CP length
Fs = 737.28;            % Sampling frequency[MHz]
Fc = Fs/4;              % Carrier frequency[MHz]
NFFT = 128;             % Number of FFT
SNR = 40;               % SNR[dB]

%% Transmmition
% Baseband signal generation
rng shuffle
dataIn = randi([0 1],R,1);                          % Pseudo vector
dataIn_M = reshape(dataIn,length(dataIn)/k,k);      % Serial to parallel
dataInSymbols = bi2de(dataIn_M);                    % Binary to decimal
dataInMod = qammod(dataInSymbols,M,0);              % Generate M-QAM symbols

SignalWithCyclicPrefix = [];
for id = 1:length(dataInMod)/NFFT
    dataInModId = dataInMod((id-1)*NFFT+1:id*NFFT);
    ifftSignal = ifft([dataInModId(end-NFFT/2+1:end);zeros((OverSamp-1)*NFFT,1);dataInModId(1:NFFT/2)]);
    SignalWithCyclicPrefixId = zeros(OverSamp*NFFT + L, 1);
    SignalWithCyclicPrefixId(1:L) = ifftSignal(end-L+1:end); %CP
    SignalWithCyclicPrefixId(L+1:end) = ifftSignal;   % CONTENT
    SignalWithCyclicPrefix = [SignalWithCyclicPrefix; SignalWithCyclicPrefixId];
end;

signalBB = SignalWithCyclicPrefix;
t = 0:1/Fs:length(signalBB)*1/Fs-1/Fs;
signalIF = signalBB .* exp(-1i*2*pi*Fc*t)';

%% AWGN
SignalNoise = awgn(signalIF,SNR,'measured'); % AWGN

%% Reveicer
SignalRx = SignalNoise .* exp(1i*2*pi*Fc*t)';
N = L + NFFT*OverSamp;
TOT = length(SignalRx)/N;
fftSignal = [];
for id = 1:TOT
    SignalRxId = SignalRx(N*(id-1) + 1 : N*id);
    SignalRxIdNoCP = SignalRxId(end - NFFT*OverSamp + 1:end); % remove CP
    fftSignalId = fft(SignalRxIdNoCP);
    fftSignal = [fftSignal;fftSignalId(end-NFFT/2+1:end);fftSignalId(1:NFFT/2)];
end;

dataDemod = qamdemod(fftSignal, M); % M-QAM demodulation
dataBin = de2bi(dataDemod,k);       % Decimal to binary
dataOut = dataBin(:);               % P/S

%% Plot
figure(1)
Fcenter = Fc;
Span = Fs/2;
[f,p] = My_Spectrum(SignalNoise,Fs,3,Fcenter,Span);
subplot(1,2,1);
plot(f,p)
xlim([Fcenter-Span/2, Fcenter+Span/2]);
ylim([-100,-40]);
xlabel('Frequency[MHz]');
ylabel('Mag.[dBFs]')
title('Reciever''s Spectrum');

I = real(fftSignal);
Q = imag(fftSignal);
subplot(1,2,2);
scatter(I,Q);
title(sprintf('Reciever''s Constelation'));

[numErr,BER] = biterr(dataIn,dataOut);
fprintf('The BER is %f%%, with %d errors\n',BER*100,numErr);
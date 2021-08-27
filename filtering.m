function [p_filt_1, p_filt_2] = filtering (p, Fs, f0)
% p: the pressure signal
% Fs: the sampling frequency
% f0: the centre frequency of the signal

pfft = fft(p);
N = length(p);
f = (0:N-1)/N*Fs;   % frequency array

% figure()
% plot(f, abs(squeeze(pfft)).^2)

% Select a harmonic
n_harmonic = 1;      % Fundamental
[p_filt_1, ~] = get_harmonic(f,f0,pfft,n_harmonic);

% Select a harmonic
n_harmonic = 2;      % Second harmonic
[p_filt_2, ~] = get_harmonic(f,f0,pfft,n_harmonic);

end
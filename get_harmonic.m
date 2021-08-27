
function [p_filt, pfft_filt] = get_harmonic(f,f0,pfft,n_harmonic)

% Select the n-th harmonic

f_min = f0*(n_harmonic - 1/2);
f_max = f0*(n_harmonic + 1/2);

% Only keep frequency content around the harmonic
pfft(f<f_min) = 0;
pfft(f>f_max) = 0;

% Make symmetryic around Fs/2 (real signal has symmetric real part Fourier
% transform, and antisymmetric imaginary part).
N = length(pfft);
pfft(ceil(N/2+1):N) = conj(pfft(floor(1+N/2):-1:2));

p_filt = real(ifft(pfft));      % Transform back to time domain
pfft_filt = pfft;
end
function [SIGbf, zvector, xvector] = RF2img (file_path, sampling_freq, line_start, line_finish)
% load('RESULTS/no_speckle/xWave_5deg.mat')

load(file_path);

for i = 1 : size(RF_matrix, 3)
   SIG(:, :, i) = RF_matrix(:, i:i+64, i);
end

       

%% Filter data 
filt_band = [5e6 25e6]; 

Bp = fir1(250, filt_band/(1/dt/2)); Ap = 1;

time_window = tukeywin(size(SIG,1),0.001);

SIGXf = zeros(size(SIG));

for k = 1 : size(SIG, 3)
    for h = 1 : size(SIG, 2)

        SIGXf(:,h,k) = time_window.*filtfilt(Bp,Ap,SIG(:,h,k));

    end
end
%% chanhging sampling frequency to Verasonics
Fs = 62.5e6; % Verasonics clock [Hz]
time_old = 0:dt:(size(SIGXf,1)-1)*dt;
signal_duration = (size(SIGXf,1)-1)*dt; % [s]
time_new = 0:1/Fs:signal_duration;

for i = 1 : size(SIGXf, 3)
    for j = 1 : size(SIGXf, 2)
        SIG_ds(:, j, i) = interp1(time_old, SIGXf(:, j, i), time_new);
    end
end

%     figure()
%     imagesc(SIG_ds(100:end, :, 129))
%%
param.p = transducer.pitch * dx; % pitch [m]
param.lambda = speed_of_sound / pulse.freq;
param.xAngle = pulse.angle;
param.f = 1.28;
param.Fc = pulse.freq;
param.c_us = speed_of_sound;
param.NC = pulse.num_cycles;
% recorded lines 
param.start = line_start; 
param.finish = line_finish;

switch sampling_freq 
    case 'original'
        param.Fs = 1 / dt;
        [SIGbf, zvector, xvector] = xWave_beamforming(SIGXf, param);
    case 'verasonics'
        param.Fs = Fs;
        [SIGbf, zvector, xvector] = xWave_beamforming(SIG_ds, param);           
end
    
end
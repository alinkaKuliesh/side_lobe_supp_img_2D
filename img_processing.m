clear all

resolution = 0.5;

switch resolution
    case 1
        load('RESULTS/xWave_lambda_pitch.mat')

        line_start = 28;
        line_finish = 46;
        
        file_mb = 'RESULTS/xWave_lambda_pitch.mat'; % mb = main beam
        file_cb = 'RESULTS/xWave_inv_ratio_lambda_pitch.mat'; % cb = comp beam
        
    case 0.5
        load('RESULTS/xWave_targ_38um.mat')
        line_start = 1;
        line_finish = 34;
        file_mb = 'RESULTS/xWave_targ_38um.mat'; % mb = main beam
        file_cb = 'RESULTS/xWave_inv_ratio_targ_38um.mat'; % cb = comp beam
        
    case 0.25
%         load('RESULTS/xWave_lambda_pitch.mat')
        load('RESULTS/xWave_targ_38um.mat')

        line_start = 109;
        line_finish = 179;
%         file_mb = 'RESULTS/xWave.mat'; % mb = main beam
        file_mb = 'RESULTS/xWave_targ_38um.mat';
        file_cb = 'RESULTS/xWave_inv_ratio_targ_38um.mat'; % cb = comp beam
        
end

% reconstruction of the image from single beam
[IQbf, z_vec, x_vec] = RF2img(file_mb, 'verasonics', line_start, line_finish, resolution);

% reconstruction of the image from double beam
[IQbf_inv, ~, ~] = RF2img(file_cb, 'verasonics', line_start, line_finish, resolution);

%% image
DR = 40;
Env = abs(IQbf); % real envelope
I = 20*log10(Env/max(Env,[],'all'));
I(I < -DR) = -DR;

Env_inv = abs(IQbf_inv); % real envelope
I_inv = 20*log10(Env_inv/max(Env_inv,[],'all'));
I_inv(I_inv < -DR) = -DR;

figure()
imagesc(x_vec, z_vec, I); hold on
axis image
colorbar
colormap gray

figure()
imagesc(x_vec, z_vec, I_inv); hold on
axis image
colorbar
colormap gray

%%
floorn = @(x,n) floor(x.*10^n)/10^n;
pos_z = find(floorn(z_vec, 4) == 5e-3, 1);

figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I(1:pos_z,line_start:line_finish)); 
% imagesc(I(1:pos_z,line_start:line_finish)); 
axis image
xlabel('x [mm]')
ylabel('z [mm]')
colorbar
colormap gray
title('Main beam, f = 15 MHz')

figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I_inv(1:pos_z,line_start:line_finish)); 
% imagesc(I_inv(1:pos_z,line_start:line_finish)); 
axis image
xlabel('x [mm]')
ylabel('z [mm]')
colorbar
colormap gray
title('Complementary beam, f = 7.5 MHz')

%% overlay points on the image
[xv, yv] = resolution_grid(2*transducer.pitch*dx);
xv = 2e-3 + xv;
yv = x_vec(ceil(end/2)) + yv;

figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I(1:pos_z, line_start:line_finish)); hold on
for i = 1 : length(xv)
    scatter(yv(i)*1e3, xv(i)*1e3, 'r', 'fill'); hold on
end
axis image
colorbar
colormap gray
title('Main beam, f = 15 MHz')


figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I_inv(1:pos_z, line_start:line_finish)); hold on
for i = 1 : length(xv)
    scatter(yv(i)*1e3, xv(i)*1e3, 'r', 'fill'); hold on
end
axis image
colorbar
colormap gray
title('Complementary beam, f = 7.5 MHz')

%% subtraction of IQ data (main IQ - complementary IQ) 
normalization = 1; % 0 = without normalisation of envelope data 1 = with

if normalization
    Env = Env ./ max(Env, [], 'all');
    Env_inv = Env_inv ./ max(Env_inv, [], 'all');

    Env_diff = Env - 0.3 * Env_inv;
else
    Env_diff = Env - Env_inv;
end

% I_diff = I_diff - min(I_diff, [], 2);
% I_diff = abs(I_diff);
Env_diff(Env_diff < 0) = 0;

I_diff = 20*log10(Env_diff/max(Env_diff,[],'all'));
I_diff(I_diff < -DR) = -DR;

figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I_diff(1:pos_z, line_start:line_finish)); hold on
axis image
title('IQbf\_main - IQbf\_comp (f = 7.5 MHz, NC = 2)')
colorbar
colormap gray

figure()
imagesc(x_vec(line_start:line_finish)*1e3, z_vec(1:pos_z)*1e3, I(1:pos_z, line_start:line_finish)); hold on
axis image
title('IQbf\_main')
colorbar
colormap gray

%% cross-section line of the PSF
% [~, pos_z] = min(abs(z_vec-xv(1)-0.025e-3));
pos_z = 145;

% in dB
figure()
plot(x_vec(line_start:line_finish)*1e3, I(pos_z, line_start:line_finish), 'DisplayName', 'Main beam'); hold on
plot(x_vec(line_start:line_finish)*1e3, I_inv(pos_z, line_start:line_finish), 'DisplayName', 'Complementary beam, f = 7.5 MHz'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [dB]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

figure()
plot(x_vec(line_start:line_finish)*1e3, I(pos_z, line_start:line_finish), 'DisplayName', 'Main beam'); hold on
plot(x_vec(line_start:line_finish)*1e3, I_diff(pos_z, line_start:line_finish), 'DisplayName', 'Main  - Complementary, f = 7.5 MHz'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope difference [dB]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

% in envelope
figure()
plot(x_vec(line_start:line_finish)*1e3, Env(pos_z, line_start:line_finish), 'DisplayName', 'Main beam'); hold on
plot(x_vec(line_start:line_finish)*1e3, Env_inv(pos_z, line_start:line_finish), 'DisplayName', 'Complementary beam, f = 7.5 MHz'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [~]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

figure()
plot(x_vec(line_start:line_finish)*1e3, Env(pos_z, line_start:line_finish), 'DisplayName', 'Main beam'); hold on
plot(x_vec(line_start:line_finish)*1e3, Env_diff(pos_z, line_start:line_finish), 'DisplayName', 'Complementary beam, f = 7.5 MHz'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [~]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend








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
        load('RESULTS/xWave.mat')
        line_start = 1;
        line_finish = 17;
        file_mb = 'RESULTS/xWave_16deg.mat'; % mb = main beam
        
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
x_vec = x_vec(line_start:line_finish);

%% image
interpolate = 1; % interpolate in longitudinal direction
DR = 40;
Env = abs(IQbf); % real envelope
if interpolate 
    Env = interp2(x_vec, z_vec', Env, interpn(x_vec), z_vec');
end
I = 20*log10(Env/max(Env, [], 'all'));
I(I < -DR) = -DR;


figure()
imagesc(x_vec, z_vec, I); hold on
axis image
colorbar
colormap gray

%%
floorn = @(x,n) floor(x.*10^n)/10^n;
% pos_z_start = find(floorn(z_vec, 4) == 2.5e-3, 1);
pos_z_start = 1;
pos_z_end = find(floorn(z_vec, 4) == 10e-3, 1);

figure()
% imagesc(x_vec*1e3, z_vec(pos_z_start:pos_z_end)*1e3, I(pos_z_start:pos_z_end, :)); 
imagesc(x_vec, z_vec, I(pos_z_start:pos_z_end, :)); 
axis image
xlabel('x [mm]')
ylabel('z [mm]')
colorbar
colormap gray
title(['Main beam, f = 15.625 MHz, NC = ', num2str(pulse.num_cycles)])

%% overlay points on the image
load('res_grid.mat')

xv = res_grid.xv - res_grid.margin * dx;
yv = res_grid.yv - (((transducer.num_active_elements - 1) / 2) * transducer.pitch + transducer.pitch / 2) * dx;

figure()
imagesc(x_vec*1e3, z_vec(pos_z_start:pos_z_end)*1e3, I(pos_z_start:pos_z_end, :)); hold on
for i = 1 : length(xv)
    scatter(yv(i)*1e3, xv(i)*1e3, 'r', 'fill'); hold on
end
axis image
colorbar
colormap gray
title(['Main beam, f = 15 MHz, NC = ', num2str(pulse.num_cycles)])


figure()
imagesc(x_vec*1e3, z_vec(pos_z_start:pos_z_end)*1e3, I_inv(pos_z_start:pos_z_end, :)); hold on
for i = 1 : length(xv)
    scatter(yv(i)*1e3, xv(i)*1e3, 'r', 'fill'); hold on
end
axis image
colorbar
colormap gray
title(['Complementary beam, f = 7.5 MHz, NC = ', num2str(pulse.num_cycles)])

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
% I_diff(I_diff < -DR) = -DR;
I_diff(I_diff < min(I, [], 'all')) = min(I, [], 'all');

figure()
imagesc(x_vec*1e3, z_vec(pos_z_start:pos_z_end)*1e3, I_diff(pos_z_start:pos_z_end, :)); hold on
axis image
if normalization
    title(['norm(E\_main) - 0.3*norm(E\_comp) (f = 7.5 MHz, NC = ', num2str(pulse.num_cycles), ')'])
else
    title(['E\_main - E\_comp (f = 7.5 MHz, NC = ', num2str(pulse.num_cycles-2),')'])
end
colorbar
colormap gray

figure()
% imagesc(x_vec*1e3, z_vec(pos_z_start:pos_z_end)*1e3, I(pos_z_start:pos_z_end, :)); hold on
imagesc(I(pos_z_start:pos_z_end, :));
axis image
title('IQbf\_main')
colorbar
colormap gray

%% cross-section line of the PSF
% [~, pos_z] = min(abs(z_vec-xv(end)-0.0128e-3));
pos_z = 248 + pos_z_start;

% envelope [dB]
figure()
plot(x_vec*1e3, I(pos_z, :), 'DisplayName', 'E\_main'); hold on
plot(x_vec*1e3, I_inv(pos_z, :), 'DisplayName', 'E\_comp'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [dB]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

figure()
plot(x_vec*1e3, I(pos_z,:), 'DisplayName', 'E\_main'); hold on
plot(x_vec*1e3, I_diff(pos_z, :), 'DisplayName', 'norm(E\_main) - 0.3*norm(E\_comp)'); 
% plot(x_vec*1e3, I_diff(pos_z, :), 'DisplayName', 'E\_main - E\_comp(NC = 2)'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope difference [dB]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

% envelope [~]
figure()
plot(x_vec*1e3, Env(pos_z, :), 'DisplayName', 'E\_main'); hold on
plot(x_vec*1e3, Env_inv(pos_z, :), 'DisplayName', 'E\_comp'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [~]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend

figure()
plot(x_vec*1e3, Env(pos_z, :), 'DisplayName', 'E\_main'); hold on
plot(x_vec*1e3, Env_diff(pos_z, :), 'DisplayName', 'E\_comp'); 
xlabel('longitudinal direction [mm]')
ylabel('Envelope [~]')
title (strcat('PSF at z =  ', num2str(z_vec(pos_z)*1e3), ' mm'))
legend








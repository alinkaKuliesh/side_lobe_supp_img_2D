clear all

line_start = 1;
line_finish = 23;

theta_start = 5;
theta_end = 5;
theta_step = 5;

n = (theta_end - theta_start) / theta_step + 1;

for i = 1 : 1 : n % 15 - 21 step 0.5
    theta = theta_start + (i - 1) * theta_step;
    display(theta)
    file = strcat('RESULTS/no_speckle/xWave_', num2str(theta), 'deg.mat');
    % reconstruction of the image from single beam
    if i == 1
        [IQbf(i, :, :), z_vec, x_vec] = RF2img(file, 'verasonics', line_start, line_finish);
        x_vec = x_vec(line_start:line_finish);    
    else
        [IQbf(i, :, :), ~, ~] = RF2img(file, 'verasonics', line_start, line_finish);
    end
end

%% image
interpolate = 1; % interpolate in longitudinal direction
DR = 50;

for i = 1 : size(IQbf, 1)
    
    Env = abs(squeeze(IQbf(i, :, :))); % real envelope
    
    if interpolate 
        Env = interp2(x_vec, z_vec', Env, interpn(x_vec), z_vec');
    end
    
    I = 20*log10(Env/max(Env, [], 'all'));
    I(I < -DR) = -DR;
    
    Img(i, :, :) = I;

    figure()
    imagesc(x_vec*1e3, z_vec*1e3, I);
%     imagesc(I);
    axis image
    xlabel('x [mm]')
    ylabel('z [mm]')
    if i == size(IQbf, 1)
        colorbar
    end
    colormap gray
    title(['\theta = ', num2str(theta_start + (i - 1) * theta_step)])
end
%%
idx = 528; % 454

figure()
for i = 1 : size(Img, 1)
    plot(interpn(x_vec)*1e3, squeeze(Img(i, idx, :)), 'DisplayName',...
        num2str(theta_start+(i-1)*theta_step)); hold on
end
xlabel('x [mm]')
ylabel('amplitude [dB]')
legend



%% overlay points on the image
x_gt = 8 * 51.2e-6;
z_gt = [1:10:100] * 1480 / 15.625e6; 

figure()
imagesc(x_vec*1e3, z_vec*1e3, I); hold on
for i = 1 : length(z_gt)
    scatter(x_gt*1e3, z_gt(i)*1e3, 8, 'r', 'fill'); hold on
end
axis image
colorbar
colormap gray
title(['\theta = ', num2str(theta_start + (n - 1) * 0.5)])

%% coherent compounding
IQbf_cc = squeeze(sum(IQbf, 1));

Env = abs(IQbf_cc); % real envelope

if interpolate 
    Env = interp2(x_vec, z_vec', Env, interpn(x_vec), z_vec');
end

I = 20*log10(Env/max(Env, [], 'all'));
I(I < -DR) = -DR;

figure()
imagesc(x_vec*1e3, z_vec*1e3, I); 
axis image
xlabel('x [mm]')
ylabel('z [mm]')
colorbar
colormap gray
title('CC 13 angles [15:21]')

%% cross-section
idx = 528; % 454

figure()
plot(interpn(x_vec)*1e3, squeeze(I(idx, :)), 'DisplayName', 'CC'); hold on
plot(interpn(x_vec)*1e3, squeeze(Img(1, idx, :)), 'DisplayName', '15 deg'); hold on
plot(interpn(x_vec)*1e3, squeeze(Img(end, idx, :)), 'DisplayName', '21 deg');
xlabel('x [mm]')
ylabel('amplitude [dB]')
legend

%% 3D maps of PSF
idx_start = 440;
idx_end = 468;


figure;
%surf 1
surf(squeeze(Img(end, idx_start:idx_end, :)),'FaceLighting','gouraud',...
    'MeshStyle','column',...
    'SpecularColorReflectance',0,...
    'SpecularExponent',5,...
    'SpecularStrength',0.2,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'FaceAlpha',0.2,...
    'FaceColor',[0.07 0.6 1],...
    'EdgeAlpha',0.2,...
    'DisplayName','21 deg');
hold on
%surf 2
surf(I(idx_start:idx_end, :),'SpecularExponent',1,...
    'SpecularStrength',1,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'FaceColor',[0.5 0.5 .5],...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'EdgeAlpha',1,...
    'DisplayName','CC');
legend
xlabel('x [mm]')
ylabel('z [mm]')
zlabel('amplitude [dB]')

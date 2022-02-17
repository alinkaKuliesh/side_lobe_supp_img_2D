% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% for this study the media is homogeneus
% =========================================================================

function medium = define_medium_beam(kgrid, margin, speckle_flag)
% average tissue properties
c0 = 1540; %[m/s]
rho0 = 1000; %[kg/m^3]
BonA0 = 6;
alpha_coeff_0 = 0.75; % power law absorption prefactor [dB/(MHz^y cm)]

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.02;
background_map = background_map_mean +...
    background_map_std*rand([kgrid.Nx, kgrid.Ny]);
% add speckle
if speckle_flag 
    maps.sound_speed= c0 * ones(kgrid.Nx, kgrid.Ny) .* background_map;
    maps.density = rho0 * ones(kgrid.Nx, kgrid.Ny) .* background_map;
else
    maps.sound_speed= c0 * ones(kgrid.Nx, kgrid.Ny);
    maps.density = rho0 * ones(kgrid.Nx, kgrid.Ny);    
end

maps.BonA = BonA0 * ones(kgrid.Nx, kgrid.Ny);
maps.alpha_coeff = alpha_coeff_0 * ones(kgrid.Nx, kgrid.Ny);
%% fat inclusion
c_fat = 1450; %[m/s]
rho_fat = 950; %[kg/m^3]
BonA_fat = 10;
alpha_coeff_fat = 0.6; % power law absorption prefactor [dB/(MHz^y cm)]

x_mask = margin + 250;
y_mask = kgrid.Ny / 2 - 130; 

ball = zeros(kgrid.Nx, kgrid.Ny);
radius = ceil(1e-3/kgrid.dx);
        
ball = ball + makeDisc(kgrid.Nx, kgrid.Ny, x_mask, y_mask, radius, 1);

maps.sound_speed(ball >= 1) = c_fat;
maps.density(ball >= 1) = rho_fat;
maps.BonA(ball >= 1) = BonA_fat;
maps.alpha_coeff(ball >= 1) = alpha_coeff_fat;

figure()
imagesc([0:size(maps.sound_speed, 2)-1]*kgrid.dx*1e3, [0:size(maps.sound_speed, 1)-1]*kgrid.dx*1e3, maps.sound_speed)
axis image
xlabel('x [mm]')
ylabel('z [mm]')
title('slice of SOS map')
colorbar

% define medium properties
medium.sound_speed = maps.sound_speed;
medium.density = maps.density;
medium.BonA = maps.BonA;
medium.alpha_coeff = maps.alpha_coeff; 
medium.alpha_power = 1.5;  % frequency dependance y
end
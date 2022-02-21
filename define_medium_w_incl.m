% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% for this study the media is homogeneus
% =========================================================================

function medium = define_medium_w_incl(kgrid, transducer, pulse, scatters, margin, speckle_flag)
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
%% blood inclusion
c_blood = 1584; %[m/s]
rho_blood = 1060; %[kg/m^3]
BonA_blood = 6;
alpha_coeff_blood = 0.14; % power law absorption prefactor [dB/(MHz^y cm)]

x_mask = margin + 250;
y_mask = kgrid.Ny / 2 - 130; 

ball = zeros(kgrid.Nx, kgrid.Ny);
radius = ceil(0.5e-3/kgrid.dx);
        
ball = ball + makeDisc(kgrid.Nx, kgrid.Ny, x_mask, y_mask, radius, 1);

maps.sound_speed(ball >= 1) = c_blood;
maps.density(ball >= 1) = rho_blood;
maps.BonA(ball >= 1) = BonA_blood;
maps.alpha_coeff(ball >= 1) = alpha_coeff_blood;
%% resolution grid
switch scatters
    case 'single points'
        xv = [1:10:100] * pulse.wave_length + margin * kgrid.dx ; % [m]
        if mod (transducer.pitch, 2) == 0
             yv = (((transducer.num_active_elements - 1) / 2 + 11) * transducer.pitch ...
                 + transducer.pitch / 2) * kgrid.dy;     % [m]  
        else
            yv = (((transducer.num_active_elements - 1) / 2 + 11) * transducer.pitch...
                + (transducer.pitch - 1) / 2) * kgrid.dy;  
        end   
        
    case 'resolution grid'
       % min distance between points 100 um 
       [xv, yv] = resolution_grid(2*transducer.pitch*kgrid.dy);
       xv = 4e-3 + margin * kgrid.dx + xv;
       yv = (kgrid.Ny / 2 - 100) * kgrid.dy + yv;           
         
end

x_mask = round(xv/kgrid.dx);
y_mask = round(yv/kgrid.dy); 

shape = 'disc';

switch shape
    case 'cube'
        x_mask = [x_mask - 1, x_mask, x_mask + 1];
        x_mask = repmat(x_mask, 1, 3);

        y_mask = repmat(y_mask, 1, 3);
        y_mask = [y_mask - 1, y_mask, y_mask + 1];

        for i = 1 : length(x_mask)
            maps.sound_speed(x_mask(i), y_mask(i), z_mask(k)) = 2 * c_water;
            maps.density(x_mask(i), y_mask(i), z_mask(k)) = 2 * rho_water;
        end
        
    case 'disc'
        ball_grid = zeros(kgrid.Nx, kgrid.Ny);
        radius_grid = ceil(pulse.wave_length/10/kgrid.dx);
        
        for i = 1 : length(x_mask)
            ball_grid = ball_grid + makeDisc(kgrid.Nx, kgrid.Ny, x_mask(i), y_mask(i), radius_grid, 0);
        end    
        
        maps.sound_speed(ball_grid >= 1) = 2 * c0;
        maps.density(ball_grid >= 1) = 15 * rho0;
        
    case 'none'
        for i = 1 : length(x_mask)
            maps.sound_speed(x_mask(i), y_mask(i)) = 2 * c0;
            maps.density(x_mask(i), y_mask(i)) = 15 * rho0;
        end
end


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
% medium.alpha_coeff = maps.alpha_coeff; 
% medium.alpha_power = alpha_power_water;  % frequency dependance y
end
% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% for this study the media is homogeneus
% =========================================================================

function medium = define_medium(kgrid, transducer, pulse, scatters, margin, speckle_flag)

    c_water = 1480; %[m/s]
    rho_water = 1000; %[kg/m^3]
    alpha_coeff_water = 2.17e-3; % power law absorption prefactor [dB/(MHz^y cm)]
    alpha_power_water = 2;  % frequency dependance y
    
    
    speckle_perc = 0.05; % 5% variation
    speckle_mask = (rand(kgrid.Nx, kgrid.Ny) - 0.5) * speckle_perc * 2 + 1;
    
    
    if speckle_flag
        % add speckle to the media
        maps.sound_speed= c_water * ones(kgrid.Nx, kgrid.Ny) .* speckle_mask;
        maps.density = rho_water * ones(kgrid.Nx, kgrid.Ny) .* speckle_mask;
    else
        maps.sound_speed= c_water * ones(kgrid.Nx, kgrid.Ny);
        maps.density = rho_water * ones(kgrid.Nx, kgrid.Ny);
    end
    
    maps.alpha_coeff = alpha_coeff_water * ones(kgrid.Nx, kgrid.Ny);

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
       xv = 3e-3 + margin * kgrid.dx + xv;
       if mod (transducer.pitch, 2) == 0           
           yv = (((transducer.num_active_elements - 1) / 2 + 7) * transducer.pitch + transducer.pitch / 2)...
               * kgrid.dy + yv;           
       else
           yv = (((transducer.num_active_elements - 1) / 2 + 7) * transducer.pitch + (transducer.pitch - 1) / 2)...
               * kgrid.dy + yv;      
       end   
end

    x_mask = round(xv/kgrid.dx);
    start_index_y = kgrid.Ny/2 - round(transducer.size_y/2) + 1; 
    y_mask = round(yv/kgrid.dy) + start_index_y; 

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
        ball = zeros(kgrid.Nx, kgrid.Ny);
        radius = ceil(pulse.wave_length/10/kgrid.dx);
        
        for i = 1 : length(x_mask)
            ball = ball + makeDisc(kgrid.Nx, kgrid.Ny, x_mask(i), y_mask, radius, 0);
        end    
        
        maps.sound_speed(ball >= 1) = 2 * c_water;
        maps.density(ball >= 1) = 15 * rho_water;
        
    case 'none'
        for i = 1 : length(x_mask)
            maps.sound_speed(x_mask(i), y_mask(i)) = 2 * c_water;
            maps.density(x_mask(i), y_mask(i)) = 15 * rho_water;
        end
end



    figure()
    imagesc([0:size(maps.sound_speed, 1)-1]*kgrid.dx*1e3, [0:size(maps.sound_speed, 2)-1]*kgrid.dx*1e3, maps.sound_speed)
    axis image
    xlabel('x [mm]')
    ylabel('z [mm]')
    title('slice of SOS map')
    colorbar
      
    % define medium properties
    medium.sound_speed = maps.sound_speed;
    medium.density = maps.density;
    medium.alpha_coeff = maps.alpha_coeff; 
    medium.alpha_power = alpha_power_water;  % frequency dependance y
end
function transducer = define_transducer(ultrasound_probe, grid_size)

switch ultrasound_probe
    case 'P6-3'
        transducer.num_elements = 128;  % number of transducer elements f# = 1.28
        transducer.num_active_elements = 65; % number of firing transducer elements 
        kerf = 0; % [m]
        if kerf < grid_size
            transducer.kerf = 0;
        else
            transducer.kerf = ceil(kerf/grid_size); % [voxels]
        end
        transducer.element_width = ceil(220e-6 / grid_size); % [voxels]
        transducer.pitch = transducer.kerf + transducer.element_width; % [voxels]
        transducer.element_length = ceil(1e-3 / grid_size); % [voxels]
        transducer.center_freq = 3e6; % [Hz]

    case 'L22-14v'
        transducer.num_elements = 128;  % number of transducer elements f# = 1.28
        transducer.num_active_elements = 65; % number of firing transducer elements 
        kerf = 20e-6; % [m]
        if kerf < grid_size
            transducer.kerf = 0;
        else
            transducer.kerf = ceil(kerf/grid_size); % [voxels]
        end
        transducer.element_width = ceil(80e-6 / grid_size); % [voxels]
        transducer.pitch = transducer.kerf + transducer.element_width; % [voxels]
        transducer.element_length = ceil(1e-3 / grid_size); % [voxels]
        transducer.bandwidth = 0.6; % fractional bandwidth 1 = 100 %
        transducer.center_freq = 15e6; % [Hz]
        
    case 'L22-14v_lambda/4'
        transducer.num_elements = 4 * 128;  % number of transducer elements f# = 1.28
        transducer.num_active_elements = 4 * 64 + 1; % number of firing transducer elements 
        kerf = 0; % [m]
        if kerf < grid_size
            transducer.kerf = 0;
        else
            transducer.kerf = ceil(kerf/grid_size); % [voxels]
        end
        transducer.element_width = ceil(25e-6 / grid_size); % [voxels]
        transducer.pitch = transducer.kerf + transducer.element_width; % [voxels]
        transducer.element_length = ceil(1e-3 / grid_size); % [voxels]
        transducer.bandwidth = 0.6; % fractional bandwidth 1 = 100 %
        transducer.center_freq = 15e6; % [Hz]
        
    case 'L22-14v_lambda/2'
        transducer.num_elements = 145;  % number of transducer elements 
        transducer.num_active_elements = 2 * 64 + 1; % number of firing transducer elements 
        kerf = 0; % [m]
        if kerf < grid_size
            transducer.kerf = 0;
        else
            transducer.kerf = ceil(kerf/grid_size); % [voxels]
        end
        transducer.element_width = ceil(50e-6 / grid_size); % [voxels]
        transducer.pitch = transducer.kerf + transducer.element_width; % [voxels]
        transducer.bandwidth = 0.6; % fractional bandwidth 1 = 100 %
        transducer.center_freq = 15.625e6; % [Hz]        
        
    case 'RCA_Imasonic'
        transducer.num_elements_full = 128;  % real number of transducer elements (rows/columns)
        transducer.num_elements = 64;  % number of transducer elements used for xAM scan
        kerf = 25e-6; % [m]
        if kerf < grid_size
            transducer.kerf = 0;
        else
            transducer.kerf = ceil(kerf/grid_size); % [voxels]
        end
        transducer.element_width = ceil(75e-6 / grid_size); % [voxels]
        % for RCA element length = full aperture
        transducer.element_length = transducer.num_elements_full * transducer.element_width ... 
           + (transducer.num_elements_full-1) * transducer.kerf; % [voxels]
        transducer.center_freq = 15e6; % [Hz]  
    
end

transducer.size_y = transducer.element_width * transducer.num_elements ...
 + transducer.kerf * (transducer.num_elements - 1);
transducer.size_y_active = transducer.element_width * transducer.num_active_elements ...
 + transducer.kerf * (transducer.num_active_elements - 1);

end
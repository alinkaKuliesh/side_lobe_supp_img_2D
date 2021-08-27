% =========================================================================
% DEFINE THE SOURCE
% =========================================================================

function source = define_source_both_inv(margin, kgrid, transducer, pulse, speed_of_sound, indent)
    apodization_Z = true;
    apodization_Y = true;

% DEFINE THE MASK    
    x_offset = margin;
    transducer.size_y = transducer.element_width * transducer.num_elements ...
         + transducer.kerf * (transducer.num_elements - 1);
    transducer.size_y_active = transducer.element_width * transducer.num_active_elements ...
         + transducer.kerf * (transducer.num_active_elements - 1);
    
    source.u_mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
    start_index_y = kgrid.Ny/2 - round(transducer.size_y/2) + 1 +...
        indent * (transducer.element_width + transducer.kerf);
    start_index_z = kgrid.Nz/2 - round(transducer.element_length/2) + 1;
    
    if transducer.kerf ~= 0
        pattern = [ones(transducer.element_width, 1); zeros(transducer.kerf, 1)]';
        pattern_y = [repmat(pattern, 1, transducer.num_active_elements - 1) ones(1, transducer.element_width)];
        pattern_yz = repmat(pattern_y, transducer.element_length, 1);
        source.u_mask(x_offset, start_index_y:start_index_y+transducer.size_y_active-1,...
            start_index_z:start_index_z+transducer.element_length-1) = pattern_yz';
    else
        source.u_mask(x_offset, start_index_y:start_index_y+transducer.size_y_active - 1,...
            start_index_z:start_index_z+transducer.element_length-1) = 1;    
    end
   
% DEFINE THE SIGNAL
    sampling_freq = 1/kgrid.dt;   % [Hz]
    tone_burst_freq = pulse.freq;  % [Hz]
    tone_burst_cycles = pulse.num_cycles;
    element_spacing = (transducer.element_width + transducer.kerf) * kgrid.dx;   % [m]
    element_index = 0:((transducer.num_active_elements - 1) / 2  - 1);    
    element_index = [element_index fliplr(element_index)];
    
    delays = element_spacing * element_index * sind(pulse.angle) / speed_of_sound;

    tone_burst_offset = delays / kgrid.dt;

    temporaryp = pulse.pnp.*toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles, ...
        'Envelope', 'Gaussian', 'SignalOffset', tone_burst_offset);
    
    figure()
    plot([1:length(temporaryp(1,:))]*kgrid.dt*1e6, temporaryp(1,:), 'DisplayName', '400kPa'); hold on;
    xlabel('time [us]')
    ylabel('pressure [Pa]')
    title('Pulse Shape')
    legend
    
    figure();
    stackedPlot([1:length(temporaryp(1,:))]*kgrid.dt, temporaryp);
    xlabel('sampling points');
    ylabel('Transducer Element');
    title('Transmit Pressure');
    
    
    temporaryp(1:(transducer.num_active_elements-1)/2, :) = -1 * temporaryp(1:(transducer.num_active_elements-1)/2, :);

    figure();
    stackedPlot([1:length(temporaryp(1,:))]*kgrid.dt, temporaryp);
    xlabel('sampling points');
    ylabel('Transducer Element');
    title('Transmit Pressure');
    
    winY = getWin((transducer.num_active_elements-1)/2, 'Tukey', 'Param', 0.5, 'Plot', false).'; % 0 Rectangular window; 1 Hann window
%     winY = getWin(transducer.num_active_elements/2, 'Hanning', 'Plot', true).'; % 0 Rectangular window; 1 Hann window
%     winY = kaiser(transducer.num_active_elements/2, 2).';
%     winY = chebwin(transducer.num_active_elements/2, 20).';
    winY = repmat(winY, 1, 2);
    
    if apodization_Y
        temporaryp = temporaryp.*winY';
    end
    
    % add silent element in the middle
    temporaryp = [temporaryp(1:(transducer.num_active_elements-1)/2, :);...
        zeros(1, size(temporaryp, 2));...
        temporaryp((transducer.num_active_elements-1)/2+1:end, :)];
    
    figure();
    stackedPlot([1:length(temporaryp(1,:))]*kgrid.dt, temporaryp);
    xlabel('sampling points');
    ylabel('Transducer Element');
    title('Transmit Pressure');
    
    %repeat elements along y direction
    temporaryp = repelem(temporaryp, transducer.element_width, 1);
    
    winZ = getWin(transducer.element_length, 'Tukey', 'Param', 0.5, 'Plot', false).'; % 0 Rectangular window; 1 Hann window
    
    %repeat along z direction
    num_voxels = transducer.num_active_elements * transducer.element_width; % number of voxels in transducer in longitudial direction
    
    for i = 1 : transducer.element_length
        if apodization_Z
            source.ux(num_voxels*(i-1)+1:num_voxels*(i-1)+num_voxels,:) = temporaryp * winZ(i) / 1480 / 1000; % v = p / c / rho; 
        else
            source.ux(num_voxels*(i-1)+1:num_voxels*(i-1)+num_voxels,:) = temporaryp / 1480 / 1000; 
        end
    end
    
%     voxelPlot(source.u_mask)
    
end
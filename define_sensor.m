% =========================================================================
% DEFINE THE SENSOR
% =========================================================================
function sensor = define_sensor(margin, kgrid, transducer, indent, receive_el)
x_offset = margin;

sensor.mask = zeros(kgrid.Nx, kgrid.Ny);

switch receive_el
    case 'all'
        start_index_y = kgrid.Ny/2 - round(transducer.size_y/2) + 1; 
    
        if transducer.kerf ~= 0
            pattern = [ones(transducer.element_width, 1); zeros(transducer.kerf, 1)]';
            pattern_y = [repmat(pattern, 1, transducer.num_elements - 1) ones(1, transducer.element_width)];
            sensor.mask(x_offset, start_index_y:start_index_y+transducer.size_y-1) = pattern_y;
        else
            sensor.mask(x_offset, start_index_y:start_index_y+transducer.size_y - 1) = 1;
        end
        
    case 'active'
        start_index_y = kgrid.Ny/2 - round(transducer.size_y_active/2) + 1 + ...
            indent * transducer.pitch;
    
        if transducer.kerf ~= 0
            pattern = [ones(transducer.element_width, 1); zeros(transducer.kerf, 1)]';
            pattern_y = [repmat(pattern, 1, transducer.num_active_elements - 1) ones(1, transducer.element_width)];
            sensor.mask(x_offset, start_index_y:start_index_y+transducer.size_y_active-1) = pattern_y;
        else
            sensor.mask(x_offset, start_index_y:start_index_y+transducer.size_y_active - 1) = 1;
        end
        
end
 
sensor.record={'p'}; 
    
%     voxelPlot(sensor.mask)
end

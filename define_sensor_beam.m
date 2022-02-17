% =========================================================================
% DEFINE THE SENSOR
% =========================================================================
function sensor = define_sensor_beam(kgrid, margin)
sensor.mask = zeros(kgrid.Nx, kgrid.Ny);
sensor_step = 1;

% define mask
sensor.mask(margin:sensor_step:end, 1:sensor_step:end) = 1;

sensor.record={'p_max'}; 

%     voxelPlot(sensor.mask)
end

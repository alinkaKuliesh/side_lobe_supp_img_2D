clear all;

flag = 'both';

MASK_PLANE = 'xy'; % only 'xz', we are looking at plane of PWs intersection

ultrasound_probe = 'L22-14v'; % options: L22-14v RCA_Imasonic

% simulation settings
gpu_run = false;
cpu_run_cluster = false;
record_movie = false;
maroilles_gpu_run = true;
 
if gpu_run
    addpath("/tudelft.net/staff-bulk/tnw/IST/AK/hpc/akuliesh1/Matlab_Toolboxes");
    DATA_CAST = 'gpuArray-single';
elseif cpu_run_cluster
    addpath("/tudelft.net/staff-bulk/tnw/IST/AK/hpc/akuliesh1/Matlab_Toolboxes");
    DATA_CAST = 'single';
elseif maroilles_gpu_run
    addpath(genpath("~/k-wave-toolbox-version-1.3"));
    DATA_CAST = 'gpuArray-single';
    DATA_PATH = '~';
    DEVICE_NUM = 1; % select gpu device [0 1 2] 
else
    DATA_CAST = 'single';
end

CFL = 0.3;
grid_size = 10e-6; % [m]

speed_of_sound = 1480; % [m/s] reference SoS
rho = 1000; % [kg/m^3] reference density

% pulse settings
% pulse.freqs = [18.75e6 11.25e6]; % [Hz]
pulse.center_freq = 15e6; % [Hz]
pulse.freq = 18e6; %[Hz]
pulse.num_cycles = 4; 
pulse.wave_length = speed_of_sound / pulse.center_freq;
pulse.length = pulse.num_cycles * pulse.wave_length; % [m]
pulse.pnp = 400e3; % [Pa]
pulse.angle = 16; % [deg]

% transducer settings based on US probe selection
switch ultrasound_probe
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
        transducer.element_length = ceil(1e-3 / grid_size); % [voxels]
        transducer.bandwidth = 0.6; % fractional bandwidth 1 = 100 %
        transducer.center_freq = 15e6; % [Hz]

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

image_depth = ceil(6e-3 / grid_size); % [voxels]

% define the grid: margin = minimal distance between transducer and PML
[kgrid, margin, PML] = define_grid(grid_size, pulse, transducer, image_depth);

% define the medium
medium = define_medium(kgrid);
% volumeViewer(medium.sound_speed)

% create the time array long enough for tilted PW propagation
% t_end = ((kgrid.Nx - margin) * kgrid.dx + ...
%     (transducer.element_width + transducer.kerf) * transducer.num_elements / 2 * tand(pulse.angle) * kgrid.dx ...
%      + pulse.length) / cosd(pulse.angle) / speed_of_sound; % [s]
t_end = (((kgrid.Nx - margin) + ...
    hypot((transducer.element_width + transducer.kerf) * transducer.num_elements / 2, 0.5*(kgrid.Nx - margin)))  * kgrid.dx ...
     + pulse.length) / speed_of_sound; % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% define the sensor
sensor = define_sensor(margin, kgrid, transducer);
[sensor_row, sensor_col] = find(squeeze(sensor.mask(margin,:,:)) == 1);
sensor_dim = [length(unique(sensor_row)), length(unique(sensor_col))];

% define the source INSIDE choose APODIZATION on/off
switch flag
    case 'left' % left half of the aperture
        source = define_source_left(margin, kgrid, transducer, pulse, speed_of_sound);
    case 'right' % right half of the aperture
        source = define_source_right(margin, kgrid, transducer, pulse, speed_of_sound);    
    case 'both' % cross-propagating PWs
        % indent = indentation from left edge
        indent = 0;
        source = define_source_both(margin, kgrid, transducer, pulse, speed_of_sound, indent);
end

% estimate required amount of memory
mem = memory_usage_estimation(kgrid.Nx, kgrid.Ny, kgrid.Nz, numel(source.ux), numel(find(sensor.mask == 1)));
display(['Minimum required memory in Gb ', num2str(mem.min)]);
%%
% =========================================================================
% RUN THE SIMULATION
% =========================================================================
if gpu_run
    input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE, PML.Z_SIZE],...
        'DataCast', DATA_CAST, 'Smooth', false};
    
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    
    save(strcat('kwaveSim_', num2str(grid_size*1e6), 'um_', MASK_PLANE, '_', flag, '.mat'),...
        'sensor_data', 'kgrid', 'Nj', 'j_vec', 'j_label', 'margin', 'MASK_PLANE', 'sensor_dim', 'sensor_step', '-v7.3')
    
elseif maroilles_gpu_run
    input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize', [PML.X_SIZE, PML.Y_SIZE, PML.Z_SIZE],...
        'DataCast', DATA_CAST, 'Smooth', false, 'DataPath', DATA_PATH, 'DeviceNum', DEVICE_NUM};
    
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
    
    save(strcat('kwaveSim_', num2str(grid_size*1e6), 'um_', '.mat'),...
        'sensor_data', 'transducer', 'kgrid', 'margin', 'sensor_dim', 'pulse', 'speed_of_sound', 'image_depth', '-v7.3')

elseif cpu_run_cluster
    input_args = {'PMLInside', false,'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE, PML.Z_SIZE],...
        'DataCast', DATA_CAST, 'Smooth', false};
    
    sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
    
    save(strcat('kwaveSim_', num2str(grid_size*1e6), 'um_', MASK_PLANE, '_', flag, '.mat'),...
        'sensor_data', 'kgrid', 'Nj', 'j_vec', 'j_label', 'margin', 'MASK_PLANE', 'sensor_dim', 'sensor_step', '-v7.3')
else
    if record_movie
        input_args = {'DisplayMask', source.u_mask, 'PMLInside', false, 'PlotPML', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE, PML.Z_SIZE],...
            'Smooth', false, 'DataCast', DATA_CAST, 'PlotScale', [-1, 1] *  max(source.ux*rho*speed_of_sound, [], 'all'), ...
            'RecordMovie', true, 'MovieName',strcat('3D_', num2str(grid_size*1e6), 'um_', 'PML_X_20'), 'MovieProfile', 'MPEG-4'};
    else
         input_args = {'DisplayMask', source.u_mask, 'PMLInside', false, 'PlotPML', true, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE, PML.Z_SIZE],...
            'Smooth', false, 'DataCast', DATA_CAST, 'PlotScale', [-1, 1] * max(source.ux*rho*speed_of_sound, [], 'all')};
    end
    
    sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
    
end

%%
clear all
load('RESULTS/kwaveSim_10um_.mat')

sensor_data.p = reshape(sensor_data.p, [sensor_dim, kgrid.Nt]);

RF = define_RF(sensor_data, transducer);

% plot the recorded time series
figure;
stackedPlot(kgrid.t_array*1e6, RF);
xlabel('Time [\mus]');
ylabel('Transducer Element');
title('Recorded Pressure');

figure;
plot(kgrid.t_array*1e6, RF(48,:));
xlabel('Time [\mus]');
ylabel('Transducer Element');
title('Recorded Pressure');


%%
param.p = (transducer.element_width + transducer.kerf) * kgrid.dx; % pitch [m]
param.lambda = pulse.wave_length;
param.xAngle = pulse.angle;
param.Fs = 1 / kgrid.dt;
param.f = 1.28;
param.Fc = pulse.freq;
param.c_us = speed_of_sound;
param.NC = pulse.num_cycles;

[IQbf, z_vector] = xWave_beamforming(RF, param);

figure()
plot(z_vector*1e3, abs(IQbf))


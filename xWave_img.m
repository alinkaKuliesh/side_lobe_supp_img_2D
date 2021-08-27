clear all

ultrasound_probe = 'L22-14v_lambda/2'; % options: L22-14v RCA_Imasonic P6-3 L22-14v_lambda/4

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
ratio = 2; % f1 / f2
pulse.center_freq = 15e6; % [Hz]
% pulse.freq = (1 + (ratio - 1) / (ratio + 1)) * pulse.center_freq; %[Hz]
pulse.freq = 15e6; % [Hz]
pulse.num_cycles = 4; 
pulse.wave_length = speed_of_sound / pulse.center_freq;
pulse.length = pulse.num_cycles * pulse.wave_length; % [m]
pulse.pnp = 400e3; % [Pa]
pulse.angle = 16; % [deg1]

% transducer settings based on US probe selection
transducer = define_transducer(ultrasound_probe, grid_size);

image_depth = ceil(6e-3 / grid_size); % [voxels]

% define the grid: margin = minimal distance between transducer and PML
[kgrid, margin, PML] = define_grid(grid_size, pulse, transducer, image_depth);

% define the sensor
sensor = define_sensor(margin, kgrid, transducer, 0, 'all');



% indent = indentation from left edge in transducer elements
%for indent = 0 : (transducer.num_elements - transducer.num_active_elements)
% indent = 0 corresponds to line/element 33; indent = 32 <-> 65 
% indent = 128 <-> central line 129
for indent = 0 : 33
% define the medium | in loop so speckle randomly generated every time |
% patterns: single points OR resolution grid
    speckle_flag = false;
    medium = define_medium(kgrid, transducer, pulse, 'resolution grid', margin, speckle_flag);
    
% create the time array long enough for tilted PW propagation
    t_end = (((kgrid.Nx - margin) + ...
        hypot(transducer.pitch * transducer.num_elements / 2, (kgrid.Nx - margin)))  * kgrid.dx ...
        + pulse.length) / speed_of_sound; % [s]
    kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);
    
% define the source INSIDE choose APODIZATION on/off
    source = define_source_both(margin, kgrid, transducer, pulse, speed_of_sound, indent);
    
    if gpu_run
        input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
            'DataCast', DATA_CAST, 'Smooth', false};
        sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    elseif maroilles_gpu_run
        input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize', [PML.X_SIZE, PML.Y_SIZE],...
            'DataCast', DATA_CAST, 'Smooth', false, 'DataPath', DATA_PATH, 'DeviceNum', DEVICE_NUM};
        sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
    else
        input_args = {'DisplayMask', source.u_mask, 'PMLInside', false, 'PlotPML', true, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
            'Smooth', false, 'DataCast', DATA_CAST, 'PlotScale', [-1, 1] * max(source.ux*rho*speed_of_sound, [], 'all')};
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    end
    
% transform sensor data to RF lines    
    RF = define_RF(sensor_data, transducer);
    RF_matrix(:, :, indent + 1) = RF';
    
end

dt = kgrid.dt;
dx = kgrid.dx;

save('xWave.mat', 'transducer', 'pulse', 'dt', 'speed_of_sound', 'dx', 'RF_matrix', '-v7.3');

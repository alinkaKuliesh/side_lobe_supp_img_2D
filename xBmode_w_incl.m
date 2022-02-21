clear all

ultrasound_probe = 'L22-14v'; % options: L22-14v RCA_Imasonic P6-3 L22-14v_lambda/4

% simulation settings
gpu_run = false;
cpu_run_cluster = false;
record_movie = false;
maroilles_gpu_run = true;
 
if gpu_run
    addpath("/tudelft.net/staff-bulk/tnw/IST/AK/hpc/akuliesh1/Matlab_Toolboxes");
    DATA_CAST = 'gpuArray-single';
    folderName = strcat('theta_', num2str(theta));
    mkdir(folderName)
    DATA_PATH = strcat('~/tudbulk/side_lobe_supp_img_2D/', folderName);
elseif cpu_run_cluster
    addpath("/tudelft.net/staff-bulk/tnw/IST/AK/hpc/akuliesh1/Matlab_Toolboxes");
    DATA_CAST = 'single';
elseif maroilles_gpu_run
    addpath(genpath("~/k-wave-toolbox-version-1.3"));
    DATA_CAST = 'gpuArray-single';
    DATA_PATH = '~';
    DEVICE_NUM = 2; % select gpu device [0 1 2] 
else
    DATA_CAST = 'single';
end

CFL = 0.3;
grid_size = 10e-6; % [m]

speed_of_sound = 1540; % [m/s] reference SoS
rho = 1000; % [kg/m^3] reference density

% pulse settings
pulse.center_freq = 15.625e6; % [Hz]
pulse.freq = 15.625e6; % [Hz]
pulse.num_cycles = 4; 
pulse.wave_length = speed_of_sound / pulse.center_freq;
pulse.length = pulse.num_cycles * pulse.wave_length; % [m]
pulse.pnp = 400e3; % [Pa]
pulse.angle = 15; % [deg1]

% transducer settings based on US probe selection
kerf_flag = false;
transducer = define_transducer(ultrasound_probe, grid_size, kerf_flag);

image_depth = ceil(10e-3 / grid_size); % [voxels]

% define the grid: margin = minimal distance between transducer and PML
[kgrid, margin, PML] = define_grid(grid_size, pulse, transducer, image_depth);

% define the sensor
sensor = define_sensor(margin, kgrid, transducer, 0, 'all');

% indent = indentation from left edge in transducer elements
%for indent = 0 : (transducer.num_elements - transducer.num_active_elements)
% indent = 0 corresponds to line/element 33; indent = 32 <-> 65 
% indent = 128 <-> central line 129
for indent = 0 : 10%2*(transducer.num_elements-transducer.num_active_elements)-1
    display(indent)
% define the medium | in loop so speckle randomly generated every time |
% patterns: single points OR resolution grid
    speckle_flag = true;
    medium = define_medium_w_incl(kgrid, transducer, pulse, 'resolution grid', margin, speckle_flag);
    
% create the time array long enough for tilted PW propagation
    t_end = (((kgrid.Nx - margin) + ...
        hypot(transducer.pitch * transducer.num_elements / 2, (kgrid.Nx - margin)))  * kgrid.dx ...
        + pulse.length) / speed_of_sound; % [s]
    kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);
    
% define the source INSIDE choose APODIZATION on/off
    source = define_source_both(margin, kgrid, transducer, pulse, speed_of_sound, indent);
    
    if gpu_run
        input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
            'DataCast', DATA_CAST, 'Smooth', false, 'DataPath', DATA_PATH};
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

if save_flag
    save(strcat('xWave_', num2str(pulse.angle), 'deg', '.mat'),...
        'transducer', 'pulse', 'dt', 'speed_of_sound', 'dx', 'RF_matrix', '-v7.3');
else
    save(strcat('xWave_', num2str(pulse.angle), 'deg', '.mat'),...
         'pulse', 'dt', 'RF_matrix', '-v7.3');
end



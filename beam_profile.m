clear all

theta = 15;
save_flag = false;

ultrasound_probe = 'L22-14v'; % options: L22-14v RCA_Imasonic P6-3 L22-14v_lambda/4

% simulation settings
gpu_run = false;
cpu_run_cluster = false;
record_movie = true;
maroilles_gpu_run = false;

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
DEVICE_NUM = 0; % select gpu device [0 1 2] 
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
pulse.angle = theta; % [deg1]

% transducer settings based on US probe selection
kerf_flag = 0;
transducer = define_transducer(ultrasound_probe, grid_size, kerf_flag);

image_depth = ceil(10e-3 / grid_size); % [voxels]

% define the grid: margin = minimal distance between transducer and PML
[kgrid, margin, PML] = define_grid(grid_size, pulse, transducer, image_depth);

% define the sensor
sensor =  define_sensor_beam(kgrid, margin);


speckle_flag = true;
medium = define_medium_beam(kgrid, margin, speckle_flag);

% create the time array long enough for tilted PW propagation
t_end = ((kgrid.Nx - margin) * kgrid.dx + ...
    (transducer.element_width + transducer.kerf) * transducer.num_elements / 2 * tand(pulse.angle) * kgrid.dx ...
     + pulse.length) / cosd(pulse.angle) / speed_of_sound; % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% define the source INSIDE choose APODIZATION on/off
source = define_source_beam(margin, kgrid, transducer, pulse, speed_of_sound);

if gpu_run
    input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
        'DataCast', DATA_CAST, 'Smooth', false, 'DataPath', DATA_PATH};
    sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
elseif maroilles_gpu_run
    input_args = {'PMLInside', false, 'PMLAlpha', PML.Alpha, 'PMLSize', [PML.X_SIZE, PML.Y_SIZE],...
        'DataCast', DATA_CAST, 'Smooth', false, 'DataPath', DATA_PATH, 'DeviceNum', DEVICE_NUM};
    sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
elseif record_movie
    input_args = {'DisplayMask', source.u_mask, 'PMLInside', false, 'PlotPML', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
        'Smooth', false, 'DataCast', DATA_CAST, 'PlotScale', [-1/4, 1/4] * max(source.ux*rho*speed_of_sound, [], 'all'), ...
         'RecordMovie', true, 'MovieName', 'beam_profile_w_incl', 'MovieProfile', 'MPEG-4'};
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
else
    input_args = {'DisplayMask', source.u_mask, 'PMLInside', false, 'PlotPML', false, 'PMLAlpha', PML.Alpha, 'PMLSize',[PML.X_SIZE, PML.Y_SIZE],...
        'Smooth', false, 'DataCast', DATA_CAST, 'PlotScale', [-1/4, 1/4] * max(source.ux*rho*speed_of_sound, [], 'all')};
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
end


%%
[sensor_row, sensor_col] = find(sensor.mask == 1);
sensor_dim = [length(unique(sensor_row)), length(unique(sensor_col))]; 

max_pressure_map = reshape(sensor_data.p_max, sensor_dim);
figure()
imagesc([0:size(max_pressure_map, 2)-1]*kgrid.dx*1e3, [0:size(max_pressure_map, 1)-1]*kgrid.dx*1e3, max_pressure_map)
axis image
xlabel('x [mm]')
ylabel('z [mm]')
title('Maximum pressure map')
colorbar


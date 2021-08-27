
clear all
close all hidden
clc

addpath(genpath("~/k-wave-toolbox-version-1.3"));
DATA_CAST = 'gpuArray-single';
DATA_PATH = '~';
DEVICE_NUM = 1; % select gpu device [0 1 2] 

downsampling_YES = 1;

% =========================================================================
% SET GRID PARAMETERS
% =========================================================================

dx = 5e-6;        % grid point spacing in the x direction [m]
dy = dx;        % grid point spacing in the y direction [m]




% create the computational grid
Nx = 4096*2;%size(cmap_interp,1)+1999+12;           % number of grid points in the x (row) direction
Ny = 4096*2;%size(cmap_interp,2);            % number of grid points in the y (column) direction
kgrid = makeGrid(Nx, dx, Ny, dy);

%%

grid_size_x = Nx*dx
grid_size_y = Ny*dy


SOS_assumed_by_scanner = 1500;

tone_burst_freq_HF = 5e6;  % [Hz]
wavelength = SOS_assumed_by_scanner/tone_burst_freq_HF;
nb_pts_per_wavelength = round(wavelength/dx)  % [grid points]
pitch = wavelength;%150e-6;%wavelength*0.5; % [m]
pitch_pts = round(pitch/dx)




% define a 2D binary sensor mask
sensor.mask = zeros(Nx, Ny);

% set the CFL
cfl = 0.3;

% define the properties of the PML to allow plane wave propagation
pml_size  = 2*nb_pts_per_wavelength;  % [grid points]

% define source
source_pos_x_index = 2*pml_size;                        % [grid points]

nb_El = 128;
X_El = (0:nb_El-1) * pitch;
X_El = X_El-mean(X_El); % center position of elements
El_ind = zeros(nb_El,pitch_pts);
for iEl=1:nb_El
    El_ind(iEl,:) = (1:pitch_pts)+(iEl-1)*pitch_pts + (Ny/2-(nb_El/2)*pitch_pts);
end

array_width = pitch_pts*dy*nb_El;

margin_left_between_array_and_PML_in_wavelength = (El_ind(1,1)-pml_size-1)/nb_pts_per_wavelength;
margin_right_between_array_and_PML_in_wavelength = (Ny-El_ind(nb_El,end)-pml_size)/nb_pts_per_wavelength;




sensor.record = {'u'};
sensor.mask(source_pos_x_index,El_ind(1,1):El_ind(nb_El,end)) = 1;
nb_sensor_pts = sum(sensor.mask(:))

% figure(6)
% imagesc(sensor.mask)


tone_burst_cycles_HF = 3;
source_mag_HF = 1; % [m/s]




% define the properties of the propagation medium    
medium.sound_speed = SOS_assumed_by_scanner*ones(Nx, Ny);%1570*ones(Nx, Ny); % [m/s]
medium.density     = 1000*ones(Nx, Ny);%1064*ones(Nx, Ny);  % [kg/m^3] 

%medium.sound_speed(source_pos_x_index+1:source_pos_x_index+size(cmap_interp,1),:) = cmap_interp;
%medium.density(source_pos_x_index+1:source_pos_x_index+size(rhomap_interp,1),:) = rhomap_interp;


density_target = 2000;
sound_speed_target = 3000;


x_position_center_target = [3800*2]; % depth
                        
y_position_center_target = [0];
y_position_center_target = y_position_center_target + Ny/2;
%y_position_center_target = round(y_position_center_target/2) + Ny/2;
%x_position_center_target = round(x_position_center_target/2);

short_radius_target = dx*1.42/2;
long_radius_target = dx*1.42/2;
                        
                        
for zz=1:length(x_position_center_target)
for ii=1:Nx
    for jj=1:Ny
        temp = ((ii-(x_position_center_target(zz)))*dx/long_radius_target)^2+((jj-y_position_center_target(zz))*dx/short_radius_target)^2;
        if temp <= 1
            medium.sound_speed(ii,jj)=sound_speed_target;%*(1+0.01*randn);
            medium.density(ii,jj)=density_target;%*(1+0.01*randn);
            %sensor.mask(ii,jj)=1;
        end
    end
end
end


figure(30)
subplot 121
imagesc(medium.density)
title('density')
colorbar
axis equal
axis tight
subplot 122
imagesc(medium.sound_speed)
title('sound speed (m/s)')
colorbar
axis equal
axis tight
drawnow


x_position_center_target = x_position_center_target-source_pos_x_index;
y_position_center_target = y_position_center_target - Ny/2;
coordinates_targets = [round(x_position_center_target)*dx,round(y_position_center_target)*dx];

width_mm = (0:Ny-2*pml_size-1)*dy*1e3;
width_mm = width_mm - mean(width_mm);
depth_mm = (0:Nx-pml_size-source_pos_x_index-1)*dx*1e3;
figure(31)
subplot 121
imagesc(width_mm,depth_mm,medium.density(source_pos_x_index:Nx-pml_size,pml_size+1:Ny-pml_size))
    hold on
    plot(coordinates_targets(:,2)*1e3,coordinates_targets(:,1)*1e3,'rx','markersize',6,'linewidth',2)
xlabel('probe lateral distance (mm)')
ylabel('distance from the probe (mm)')
title('mass density (kg/m^3)')
colorbar
axis equal
axis tight
set(gca,'fontsize',16)
subplot 122
imagesc(width_mm,depth_mm,medium.sound_speed(source_pos_x_index:Nx-pml_size,pml_size+1:Ny-pml_size))
    hold on
    plot(coordinates_targets(:,2)*1e3,coordinates_targets(:,1)*1e3,'rx','markersize',6,'linewidth',2)
xlabel('width (mm)')
ylabel('distance from the probe (mm)')
title('compressional sound speed (m/s)')
colorbar
axis equal
axis tight
set(gca,'fontsize',16)
drawnow

% pause

% set end time
%t_end = (sqrt((Nx-2*pml_size)^2+(Ny-2*pml_size)^2)*dx)/SOS_assumed_by_scanner*1.6;
t_end = (coordinates_targets(1) + sqrt((nb_El/2*pitch)^2+coordinates_targets(1)^2))/SOS_assumed_by_scanner + 5e-6;

% create the time array
kgrid.t_array = makeTime(kgrid, SOS_assumed_by_scanner, cfl, t_end);
sampling_freq = 1/(kgrid.t_array(2) - kgrid.t_array(1));

%time_steps_one_line = length(kgrid.t_array);


plane_wave_angle = [0]*pi/180; %[rad]
nb_Tx = length(plane_wave_angle);
%delay_firing = [1:pitch_pts]*dy*sin(plane_wave_angle)/SOS;
%delay_firing = (sqrt(focal_length_pts^2+[0:ceil(pitch_pts/2)].^2)-focal_length_pts)*dx/SOS;
%delay_firing = [fliplr(delay_firing) delay_firing(2:end)];
%delay_firing = max(abs(delay_firing))-delay_firing;


X_El = ([1:nb_El] - 1) * pitch;

delay_firing = zeros(nb_Tx,nb_El);
%----- focusing and steering -------
for ii=1:nb_Tx
delay_firing(ii,:) = X_El.*sin(plane_wave_angle(ii))/SOS_assumed_by_scanner;
delay_firing(ii,:) = delay_firing(ii,:) + abs(min(delay_firing(ii,:)));
figure(22)
plot(delay_firing(ii,:)*1e6,'linewidth',3)
hold on
xlabel('element number','Fontsize',20)
ylabel('transmission time delay [\mus]','Fontsize',20)
set(gca,'FontSize',20)
grid on
drawnow
end
legend('Tx1','Tx2','Tx3','Tx4','Tx5')








source_waveform = source_mag_HF*toneBurst(sampling_freq, tone_burst_freq_HF, tone_burst_cycles_HF);
%source_waveform = [source_waveform zeros(1,round(max(delay_firing(:))*1.5*sampling_freq))];


%apod_win = ones(pitch_pts,1);
%apod_win = hanning(pitch_pts);
%apod_win = tukeywin(nb_El,1);
apod_win = hanning(nb_El);

source.ux = zeros(pitch_pts,length(source_waveform));


for iTx=1:nb_Tx

        fprintf(['\n Tx ' int2str(iTx) ' \n'])

    % define source
    source.u_mask = zeros(Nx, Ny);
    source.u_mask(source_pos_x_index,El_ind(1,1):El_ind(end,end)) = 1;
    nb_source_pts = sum(source.u_mask(:));

    source.ux = zeros(nb_source_pts,length(source_waveform));

    for jj=1:nb_El % loop on element index
        %source.ux(1+(jj-1)*pitch_pts:jj*pitch_pts,:) = repmat(pulse_delaying(source_waveform,delay_firing(iTx,jj),sampling_freq),pitch_pts,1)*apod_win(jj);
        source.ux(1+(jj-1)*pitch_pts:jj*pitch_pts,:) = repmat(source_waveform,pitch_pts,1)*apod_win(jj);
    end

    
    
display_mask = source.u_mask;%+sensor.mask;

%max_pressure_display = source_mag_HF*1000*SOS_assumed_by_scanner/5;

% set the input arguments
%input_args = {'PMLSize', pml_size, 'PlotSim', false};
input_args = {'PlotSim', false, 'PMLSize', pml_size,...
    'DataCast', DATA_CAST,'DataPath', DATA_PATH,'DeviceNum', DEVICE_NUM};
% input_args = {'DisplayMask', display_mask,'PlotPML', false,'PlotScale', [-max_pressure_display, max_pressure_display], 'PMLSize', pml_size,...
%      'LogScale',true,'RecordMovie', true, 'MovieName', 'Array_July2015', 'MovieType', 'image', 'PlotFreq', 120, 'MovieArgs', {'fps', 30}};
% input_args = {'DisplayMask', 'off','PlotScale', 'auto', 'PMLSize', pml_size,...
%     'PMLAlpha', pml_alpha, 'PlotPML', false, 'PlotSim', false};


% run the simulation
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});

for jj=1:nb_El % loop on receiving element index
     SIG(:,jj,iTx) = mean(sensor_data.ux(1+(jj-1)*pitch_pts:jj*pitch_pts,:),1);
end

end
save('kWave_simulation_G.mat', 'SIG', '-v7.3');

%%
%clear all
load('/Users/akuliesh1/side_lobe_supp_img/RESULTS/kWave_simulation_G.mat')

if downsampling_YES
%---------- downsampling ---------------
downsampling_factor = 4;
Fs = sampling_freq/downsampling_factor;
temp = zeros(ceil(size(SIG,1)/downsampling_factor),size(SIG,2),size(SIG,3));
    for ii=1:size(SIG,3)    
        temp(:,:,ii) = downsample(squeeze(SIG(:,:,ii)),downsampling_factor);
    end
SIG = temp;
clear temp
end

clear kgrid source sensor input_args

save kwave_uniform_1PlaneWave_5um_23August2021




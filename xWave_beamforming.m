function [IQbf, z_vector, x_vector] = xWave_beamforming(SIG, param)
%DAS delay-and-sum beamforming of I/Q signals for cross-amplitude modulation imaging
%   IQbf = xAMbeamforming2D(SIG, param).
%   Demodulates a 2D RF signal to I/Q, performs DAS beamforming, and
%   returns the beamformed I/Q signal. The signal is beamformed in and
%   below the area where two cross-propagating wavefronts intersect.
%
%   Note: SIG must be real RF data. 
%   Dimensions should be: 
%   number of samples x number of elements in aperture x number of transmit events [Ns, nAp, nTx]
%
%   Note: a Lagrange interpolation method is used.
%
%   Note: an f-number should be provided for computation of receive
%   aperture as a function of depth.
%
%   Besides SIG the following parameters should be given in the struct param:
%   - p: pitch of transducer array [m]
%   - lambda: wavelength (should be equal to p) [m]
%   - xAngle: angle of tilted planewave with repect to transducer array
%   [degrees]
%   - c_us: speed of sound in medium [m/s]
%   - Fs: sampling frequency [Hz]
%   - f: f-number used for receive aperture. I propose f = 1.28
%
%   Note: if nothing else is specified it is assumed that Fs (sampling 
%   frequency) is 4 * Fc (central frequency)
%   -------------------------------------------------------------------

%% subtracting half a pulselength
half_pulse_points = round(param.NC * param.Fs / param.Fc / 2);
SIG = SIG(half_pulse_points:end, :, :);
%% input variables 
nAp = size(SIG,2);
nTx = size(SIG,3);
Da = (nAp-1) * param.p; % size of active aperture [m]
%% Set time axis
dt = 1/param.Fs; % time per sample [s]
Nt = dt * size(SIG, 1); % recording time [s]
Ts = 1 / param.Fs; % sampling period
if ~isfield(param,'Fc')
    param.Fc = param.Fs/4; % center frequenzy [Hz]
end

tax = 0:dt:Nt - 1 * dt; % time axis[s]
%% define min and max image depth and set z vector 
z_min = 0; % [m]
% z_max = ((Da - param.p) / 2 ) * cotd(param.xAngle); % [m]
z_max = 10e-3; % [m]

z_vector = z_min:param.p/4:(z_max - 10 * param.p); % depth axis [m]
x_vector = 0:param.p:Da; % x axis [m]

%% Receive aperture using f-number
nAp_receive = cell(length(z_vector),1);

for zn = 1:length(z_vector)
% compute aperture based on depth and f#
aperture = z_vector(zn) / param.f; 

[~, idx] = min(abs(x_vector-aperture)); 
aperture_x = x_vector(idx);
nEl = aperture_x / param.p; % number of elements

% ensure aperture number is unequal
if rem(nEl, 2) == 0
    nEl = nEl + 1;
end

% defining receive aperture
nAp_receive_loop = (round(nAp/2)-floor(nEl/2)):(round(nAp/2)+floor(nEl/2));

% excluding element 33 as this is always 0 anyway
% nAp_receive_loop(nAp_receive_loop==33)=[];
nAp_receive{zn}=nAp_receive_loop;
end

%% RF to I/Q demodulation
IQ=zeros(size(SIG));
% for q = 1:nTx 
for q = param.start : param.finish
    idx = 1;
    for k = 1:nAp
        SIG_loop = SIG(:,k,q);
        assert(isreal(SIG_loop),'RF must contain real RF signals.')
        % Convert to column vector (if RF is a row vector)
        wasrow = isrow(SIG_loop);
        
        if wasrow
            SIG_loop = SIG_loop(:);
        end
        
        % creating time vector
        nl = size(SIG_loop, 1);
        tax = (0:nl-1)' / param.Fs;
        
        % Demodulating RF data
        IQ(:,k,q) = double(SIG_loop).*exp(-1i * 2 * pi * param.Fc * tax);
        
        % low-pass filter is determined by the normalized cut-off frequency Wn
        Wn = min(2*param.Fc/param.Fs, 0.5);
        [b, a] = butter(5, Wn);
        IQ(:, k, q) = filtfilt(b, a, IQ(:, k, q)) * 2;
        
        % recover initial size (if originally row vector)
        if wasrow
            IQ(:, k, q) = IQ(:, k, q).';
        end
        
        idx = idx + size(SIG, 1);
    end
    
end
%% Beamforming loop (delay and sum)
SIGdas = zeros(length(z_vector), size(IQ,3));
% for q = 1:nTx
for q = param.start : param.finish
    IQ_loop = IQ(:,:,q);
    
    BFimageLine = zeros(length(z_vector), 1);

% beamforming loop
for zn = 1:length(z_vector)
    SIG_loop = zeros(length(z_vector),1);
    
    for xn = nAp_receive{zn}
        % computing tof from point A to receiver transducer element
        tTx = (Da / 2  * sind(param.xAngle) + z_vector(zn) * cosd(param.xAngle)) / param.c_us + max(z_vector(zn) - z_max,0) / param.c_us;
        tRx = norm( [ z_vector(zn) 0]-[0 ((xn-1) * param.p - Da  / 2 )]) / param.c_us;
        tTR = tTx + tRx;
        
        % Interpolation Lagrange
        [~, k] = min(abs(tax-tTR));
        eps = (tTR+Ts)/Ts - k;
        l1 = -eps * (eps - 1) * (eps - 2) / 6;
        l2 = (eps + 1) * (eps - 1) * (eps - 2) / 2;
        l3 = -(eps + 1) * eps * (eps - 2) / 2;
        l4 = (eps + 1) * eps * (eps - 1) / 6;
        
        if (k+2 > size(IQ_loop, 1))
            SIG_interp = IQ_loop(k,xn);        
        else
            SIG_interp = IQ_loop(k-1,xn) * l1 + IQ_loop(k,xn) * l2 + IQ_loop(k+1,xn) * l3 + IQ_loop(k+2,xn) * l4;
        end
               
        if ~isreal(IQ_loop) % if IQ: phase rotations
            SIG_loop(xn ) = SIG_interp * exp(2i*pi*param.Fc*tTR);
        else
            SIG_loop(xn ) = SIG_interp;
        end
        
    end
    % summing over all elements
    BFimageLine(zn) = sum(SIG_loop,1);
end
    
    SIGdas(:, q) = BFimageLine;
end
IQbf = SIGdas; % the beamformed I/Q signal


end
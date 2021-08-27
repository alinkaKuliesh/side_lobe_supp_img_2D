% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================
function [kgrid, margin, PML] = define_grid(grid_size, pulse, transducer, image_depth)
    % set the size of the perfectly matched layer (PML)
    PML.X_SIZE = 10;            % [grid points]
    PML.Y_SIZE = 10;            % [grid points]
    PML.Alpha = 12;

    % calculate the spacing between the grid points
    dx = grid_size; % [m] 
    dy = dx;                 
    
    % set min margin between the source/sensor & PML
    margin = ceil(2*pulse.wave_length/dx); % [grid points]

    Ny = transducer.num_elements * transducer.pitch + 2 * margin + 2 * PML.Y_SIZE;
    Nx = image_depth + margin + 2 * PML.X_SIZE;
    
    authorized = 4;

    Fac = factor(Ny);
    while max(Fac)>authorized || mod(Ny,2) ~= 0
        Ny = Ny+1;
        Fac = factor(Ny);
    end
    Ny = Ny-2*PML.Y_SIZE;


    Fac = factor(Nx);
    while max(Fac)>authorized || mod(Nx,2) ~= 0
        Nx = Nx+1;
        Fac = factor(Nx);
    end
    Nx = Nx-2*PML.X_SIZE;

    display([Nx, Ny])

    % create the k-space grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
end
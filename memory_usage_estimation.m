function mem = memory_usage_estimation(Nx, Ny, Nz, input, output)
A_max = 9;
B_max = 2;

mem.min = (13 * Nx * Ny * Nz + 7 * Nx /2 * Ny * Nz) * 4 / 1024^3  + ...
    (input + output) * 4 / 1024^3;

mem.max = ((13 + A_max) * Nx * Ny * Nz + (7 + B_max) * Nx /2 * Ny * Nz) * 4 / 1024^3  + ...
    (input + output) * 4 / 1024^3;

end
function [xv, yv] = resolution_grid_old(delta)
fhandle = @(x) 0.5 * x.^1.2;

yv = [0 1 2.5 4.5 7.5 11.5];
xv1 = fliplr(fhandle(yv));
xv2 = 0 .* yv;

yv = yv(1:5);
xv1 = xv1(1:5);
xv2 = xv2(1:5);

% laterally points separated proportional to delta [um]
yv = yv * delta;
% axially the min distance lambda 
scale = delta/ xv1(end); % [m]
xv1 = xv1 * scale;

figure(); 
scatter(yv*1e6, xv1*1e6, 'filled'); hold on
scatter(yv*1e6, xv2*1e6, 'filled');
xlabel('x [um]')
ylabel('z [um]')
set(gca, 'YDir','reverse')

xv = [xv1, xv2];
yv = [yv, yv];
end
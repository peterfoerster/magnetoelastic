pkg load geopdes;

geometry_file = 'magnetoelastic';

% plot geometry
% [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '.txt']);
% plot_geometry (geometry, boundaries);

% solve for mvp
[problem_data, method_data] = setup_problem (geometry_file);

tic;
[geometry, msh, space, u] = mp_solve_magnetostatics2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot absolute value of magnetic flux density (and write .dat files optionally)
vtk_pts = {linspace(0, 1, method_data.nsub(1)), linspace(0, 1, method_data.nsub(2))};
plot_curl2d_mp (u, space, geometry, vtk_pts, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))]);
shading interp; view(2);

% write .vtk files
% sp_to_vtk_mp_curl2d (u, space, geometry, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B');

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

% static case
pkg load geopdes;

geometry_file = 'magnetoelastic_v2';

% plot geometry
% [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '_mag.txt']);
% plot_geometry (geometry, boundaries);

% specify problem and material data
[problem_data, method_data] = setup_problem (geometry_file);

tic;
[geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_weak_coupling2d (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_weak_coupling2d_mag (problem_data, method_data);
% [geometry_mec, msh_mag, space_mag, A] = mp_solve_magnetostatics2d (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, u] = mp_solve_linear_elasticity2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot absolute value of magnetic flux density (and write .dat files optionally)
% figure(1);
% plot_curl2d_mp (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))]);
% shading interp; view(2);
% write .vtk files
sp_to_vtk_mp_curl2d (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B');
sp_to_vtk (u, space_mec, geometry_mec, method_data.nsub, ['u_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'u', 'value');

% plot deformed geometry, scale u
u_plot = u * 1e3;
figure(2);
geometry_def = geo_deform_mp (u_plot, space_mec, geometry_mec);
nrbplot(geometry_def(5).nurbs, method_data.nsub);
xrt = nrbeval(geometry_def(5).nurbs, {1, 1})
xrb = nrbeval(geometry_def(5).nurbs, {1, 0})
shading interp; view(2);

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

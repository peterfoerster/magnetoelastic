pkg load geopdes; clf;

geometry_file = 'magnetoelastic';

% plot geometry
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '_mech.txt']);
plot_geometry (geometry, boundaries);
return
% solve for mvp, displacement
[problem_data, method_data] = setup_problem (geometry_file);

tic;
% [geometry, msh, space, A] = mp_solve_magnetostatics2d (problem_data, method_data);
[geometry, msh, space, u] = solve_linear_elasticity2d (problem_data, method_data);
[geometry, msh, space, u_test] = solve_linear_elasticity_test (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot deformed geometry
figure(1);
geometry_def = geo_deform (u, space, geometry);
nrbplot (geometry_def.nurbs, method_data.nsub.*method_data.nquad);
shading interp; view(2);

figure(2);
geometry_def_test = geo_deform (u_test, space, geometry);
nrbplot (geometry_def_test.nurbs, method_data.nsub.*method_data.nquad);
shading interp; view(2);
return
% plot absolute value of magnetic flux density (and write .dat files optionally)
plot_curl2d_mp (A, space, geometry, method_data.nsub.*method_data.nquad, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))]);
shading interp; view(2);

% write .vtk files
% sp_to_vtk_mp_curl2d (A, space, geometry, method_data.nsub.*method_data.nquad, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B');

% signal that the program is finished
% x = linspace(1, 20, 8000);
% Y = sin(2*pi*440*x);
% sound(Y);

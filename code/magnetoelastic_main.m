pkg load geopdes;

geometry_file = 'magnetoelastic';

% plot geometry
% [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '_mec.txt']);
% plot_geometry (geometry, boundaries);

% specify problem and material data
[problem_data, method_data] = setup_problem (geometry_file);

tic;
% [geometry, msh, space, A] = mp_solve_magnetostatics2d (problem_data, method_data);
[geometry, msh, space, u] = mp_solve_linear_elasticity2d (problem_data, method_data);
% [geometry_test, msh_test, space_test, u_test] = solve_linear_elasticity_test (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot deformed geometry
figure(1);
geometry_def = geo_deform_mp (u, space, geometry);
nrbkntplot(geometry_def(5).nurbs);
% plot_geometry (geometry_def);
shading interp; view(2);

% figure(2);
% geometry_def_test = geo_deform (u_test, space_test, geometry_test);
% nrbplot (geometry_def_test.nurbs, method_data.nsub.*method_data.nquad);
% shading interp; view(2);

% plot absolute value of magnetic flux density (and write .dat files optionally)
% plot_curl2d_mp (A, space, geometry, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))]);
% shading interp; view(2);

% write .vtk files
% sp_to_vtk_mp_curl2d (A, space, geometry, method_data.nsub.*method_data.nquad, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B');

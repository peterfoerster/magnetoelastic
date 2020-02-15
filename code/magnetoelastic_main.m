% static case
pkg load geopdes;

magnetic = 0;
mechanic = 0;
geometry = 1;
msh      = 0;
geometry_file = 'magnetoelastic_v2';

% specify problem and material data
[problem_data, method_data] = setup_problem (geometry_file);

% plot geometry and mesh
if (geometry)
   [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '_mec.txt']);
   plot_geometry (geometry, boundaries);
   if (msh)
      figure;
      plot_mesh (geometry, method_data.nsub);
   end
   % gimp: pos(475,175), size(1025,800)
end
return
tic;
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d (problem_data, method_data);
[geometry_mec, msh_mec, space_mec, msh_mag, space_mag] = mp_rebuild_coupling2d (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_weak_coupling2d (problem_data, method_data);
% [geometry_mec, msh_mag, space_mag, A] = mp_solve_magnetostatics2d (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, u] = mp_solve_linear_elasticity2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% magnetic problem
if (magnetic)
   save(['A_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1)) '.mat'], 'A');
   sp_to_vtk_mp_curl2d (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B');
end

% mechanic problem
if (mechanic)
   save(['u_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1)) '.mat'], 'u');
   u_plot = u * 1e3;
   figure;
   geometry_def = geo_deform_mp (u_plot, space_mec, geometry_mec);
   % nrbplot(geometry_def(5).nurbs, method_data.nsub);
   nrbkntplot(geometry_def(5).nurbs);
   axis([0 4.125 2.5 4.5]);
   shading faceted; view(2);
   set(gcf, 'Position', get(0, 'Screensize'));
   % gimp: pos(210,180), size(1550,800)
end

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

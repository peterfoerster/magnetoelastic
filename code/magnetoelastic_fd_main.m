% frequency domain
pkg load geopdes;

magnetic = 1;
mechanic = 1;
geometry_file = 'magnetoelastic_v3';

% specify problem and material data
omega = 2*pi*50;
[problem_data, method_data] = setup_problem_fd (geometry_file, omega);

tic;
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d_fd (problem_data, method_data);
% [geometry_mec, msh_mec, space_mec, msh_mag, space_mag] = mp_rebuild_coupling2d_fd (problem_data, method_data);
% [geometry_mec, msh_mag, space_mag, A] = mp_solve_magnetoquasistatics2d (problem_data, method_data);
% [geometry_mec, msh_mag, space_mag] = mp_rebuild_magnetoquasistatics2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% magnetic problem
if (magnetic)
   save(['A_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1)) '.mat'], 'A');
   T = 0;
   T = linspace(0, 1/50, 31);
   sp_to_vtk_mp_curl2d_fd (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B', omega, T);
end

% mechanic problem
if (mechanic)
   save(['u_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1)) '.mat'], 'u');
   i_cmp = complex(0, 1);
   for it=1:length(T)
      u_plot = real(u * exp(i_cmp*omega*T(it))) * 1e3;
      figure(it);
      geometry_def = geo_deform_mp (u_plot, space_mec, geometry_mec);
      nrbplot(geometry_def(5).nurbs, method_data.nsub);
      % nrbkntplot(geometry_def(5).nurbs);
      axis([0 4 2.5 4.5]);
      shading faceted; view(2);
      set(gcf, 'Position', get(0, 'Screensize'));
      drawnow;
      % gimp: pos(210,180), size(1550,800)
      % print(gcf, ['mec_fd_le_scale=1e3_it=' num2str(it) '.png'], '-dpng', '-S2800,2100');
   end
end

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

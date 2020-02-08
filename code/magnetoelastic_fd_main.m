% frequency domain
pkg load geopdes;

geometry_file = 'magnetoelastic_v2';

% specify problem and material data
omega = 2*pi*50;
[problem_data, method_data] = setup_problem_fd (geometry_file, omega);

tic;
[geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d_fd (problem_data, method_data);
% [geometry_mec, msh_mag, space_mag, A] = mp_solve_magnetoquasistatics2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% write .vtk files
T = linspace(0, 1/50, 10);
sp_to_vtk_mp_curl2d_fd (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B', omega, T);

% plot deformed geometry
i_cmp = complex(0, 1);
for it=1:length(T)
   u_plot = real(u * exp(i_cmp*omega*T(it))) * 1e3;
   figure(it);
   geometry_def = geo_deform_mp (u_plot, space_mec, geometry_mec);
   nrbplot(geometry_def(5).nurbs, method_data.nsub);
   axis([0 4 2.5 4.5]);
   shading faceted; view(2);
   drawnow;
end

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

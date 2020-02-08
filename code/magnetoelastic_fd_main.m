% quasistatic, i.e. frequency domain
pkg load geopdes;

geometry_file = 'magnetoelastic_v2';

% specify problem and material data
omega = 2*pi*50;
[problem_data, method_data] = setup_problem_fd (geometry_file, omega);

tic;
% [geometry_mec, msh_mec, space_mec, u, msh_mag, space_mag, A] = mp_solve_coupling2d_fd (problem_data, method_data);
[geometry_mec, msh_mag, space_mag, A] = mp_solve_magnetoquasistatics2d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot real(A*exp(j*omega*t)) = a*cos(j*omega*t + phi) -> magnitude of B (only plot plate)
% figure(1);
% t = linspace(0, 1/50, 2^4);
% for it=1:length(t)
%    A_plot = real(A.*exp(i*omega*t(it)));
%    plot_curl2d_mp (A_plot, space_mag, geometry_mec, method_data.nsub);
%    drawnow;
% end

% write .vtk files
T = linspace(0, 1/50, 10);
sp_to_vtk_mp_curl2d_fd (A, space_mag, geometry_mec, method_data.nsub, ['B_degree=' num2str(method_data.degree(1)) '_nsub=' num2str(method_data.nsub(1))], 'B', omega, T);
return

for it=1:length(t)
   A_plot = real(A.*exp(i*omega*t(it)));
   sp_to_vtk_mp_curl2d_fd (A_plot, space_mag, geometry_mec, method_data.nsub, , 'B', it);
end
return
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

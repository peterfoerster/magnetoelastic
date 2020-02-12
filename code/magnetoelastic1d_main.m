% 1D case
pkg load geopdes;

% specify problem and material data
[problem_data, method_data] = setup_problem1d();

% need to modify weak form to include double integral?
tic;
[geometry, msh, space, w] = solve_weak_coupling1d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot deformation
x = linspace(0, problem_data.l, 10);
[ew, F] = sp_eval (w, space, geometry, {x, 0});
plot(x, ew(1,:));
return
% plot deformed geometry, scale u
u_plot = u * 1e3;
figure(2);
geometry_def = geo_deform_mp (u_plot, space_mec, geometry_mec);
nrbplot(geometry_def(5).nurbs, method_data.nsub);
xrt = nrbeval(geometry_def(5).nurbs, {1, 1})
xrb = nrbeval(geometry_def(5).nurbs, {1, 0})
shading interp; view(2);

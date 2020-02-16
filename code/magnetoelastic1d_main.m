% 1D case
pkg load geopdes;

% specify problem and material data
[problem_data, method_data] = setup_problem1d();

tic;
[geometry, msh1, space1, msh2, space2, w] = solve_weak_coupling1d (problem_data, method_data);
fprintf('\ntime elapsed for solution: %d min\n', toc/60);

% plot deformation
npts = 10;
[ew1, F1] = sp_eval (w(1:space1.ndof), space1, geometry, npts);
plot(linspace(0, problem_data.l, npts), ew1);

figure;
npts = 10;
[ew2, F2] = sp_eval (w(space1.ndof+1:space1.ndof+space2.ndof), space2, geometry, npts);
plot(linspace(0, problem_data.l, npts), ew2);

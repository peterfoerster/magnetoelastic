pkg load geopdes;

geometry_file = 'magnetoelastic_v1';
% specify problem and material data
[problem_data, method_data] = setup_problem (geometry_file);
degree_ref = [3 3];
nsub_ref   = [32 32];
degree     = [2 2];
N_it       = 5;

% compute errors
compute_errh1 (problem_data, method_data, degree_ref, nsub_ref, degree, N_it);
% plot errors
plot_convergence (degree_ref, nsub_ref, degree, N_it);

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

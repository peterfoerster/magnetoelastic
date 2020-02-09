pkg load geopdes;

geometry_file = 'magnetoelastic_v2';
filename      = 'ms_v2';
% specify problem and material data
[problem_data, method_data] = setup_problem (geometry_file);
degree_ref = [4 4];
nsub_ref   = [32 32];
degree     = [3 3];
N_it       = 5;

% compute errors
compute_errh1 (problem_data, method_data, degree_ref, nsub_ref, degree, N_it, filename);
% plot errors
plot_convergence (degree_ref, nsub_ref, degree, N_it, filename);

% signal that the program is finished
x = linspace(1, 20, 8000);
Y = sin(2*pi*440*x);
sound(Y);

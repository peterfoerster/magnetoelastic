% input:
% degree, N_it, nquad_offset(standard degree+1)
% output:
% u(for every iteration)

function [] = conv_study_sol (degree, N_it, nquad_offset)
geometry_file = 'photocathode_insulator';
voltage = -60e3;
[problem_data, method_data] = init_potential (geometry_file, voltage);

% convergence study
method_data.degree = [degree degree];
method_data.regularity = [degree-1 degree-1];
method_data.nquad = [degree+1+nquad_offset degree+1+nquad_offset];
for iit=0:N_it
  fprintf('\niteration with nsub = %d\n', 2^iit);
  method_data.nsub = [2^iit 2^iit];
  tic;
  [geometry, msh, space, u] = mp_solve_laplace_mod (problem_data, method_data);
  fprintf('time elapsed for field solution:%d\n', toc);
  tic;
  filename = [geometry_file '_degree=' num2str(degree) '_nsub=' num2str(2^iit) '_nquad_offset=' num2str(nquad_offset) '.mat'];
  save(filename, 'u');
end
end

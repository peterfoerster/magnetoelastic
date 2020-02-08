function [] = compute_errh1 (problem_data, method_data, degree_ref, nsub_ref, degree, N_it)
   % change filename and solver calls for other studies
   filename = 'magnetostatics2d';

   errh1 = errl2 = NaN(N_it+1,1);
   % reference solution
   method_data.degree     = degree_ref;
   method_data.regularity = degree_ref-1;
   method_data.nsub       = nsub_ref;
   method_data.nquad      = degree_ref+1;
   tic;
   [geometry, msh_ref, space_ref, u_ref] = mp_solve_magnetostatics2d (problem_data, method_data);
   fprintf('\ntime elapsed for reference solution:%d min\n', toc/60);

   % iterative solution and error computation
   for iit=0:N_it
      fprintf('\niteration with nsub = %d\n', 2^iit);
      method_data.degree     = degree;
      method_data.regularity = degree-1;
      method_data.nsub       = [2^iit 2^iit];
      method_data.nquad      = degree+1;
      [geometry, msh, space, u] = mp_solve_magnetostatics2d (problem_data, method_data);
      tic;
      [errh1(iit+1), errl2(iit+1)] = mp_errh1 (geometry, msh_ref, space_ref, u_ref, space, u);
      fprintf('\ntime elapsed for error computation:%d min\n', toc/60);
   end

   % save results
   filename = [filename '_degree_ref=' num2str(degree_ref(1)) '_nsub_ref=' num2str(nsub_ref(1)) ...
               '_degree=' num2str(degree(1)) '_N_it=' num2str(N_it) '.mat'];
   save(filename, 'errh1', 'errl2');
end

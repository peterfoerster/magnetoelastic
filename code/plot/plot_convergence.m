function [] = plot_convergence (degree_ref, nsub_ref, degree, N_it)
   % change filename corresponding with 'compute_errh1'
   filename = 'magnetostatics2d';

   filename = [filename '_degree_ref=' num2str(degree_ref(1)) '_nsub_ref=' num2str(nsub_ref(1)) ...
               '_degree=' num2str(degree(1)) '_N_it=' num2str(N_it) '.mat'];
   % 'errh1', 'errl2'
   load(filename);

   h = 1./(2.^(0:N_it));
   h_l2 = h.^(degree(1)+1);
   h_h1 = h.^degree(1);

   figure(1);
   loglog(h, errl2, h, h_l2*errl2(1));
   legend('computed', 'p+1');
   xlabel('h');
   title('relative error in L^2 norm');

   figure;
   loglog(h, errh1, h, h_h1*errh1(1));
   legend('computed', 'p');
   xlabel('h');
   title('relative error in H^1 norm');
end

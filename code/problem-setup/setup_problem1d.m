function [problem_data, method_data] = setup_problem1d ()
   % geometry
   l = 0.3;
   geometry = nrbsquare ([0 0], l, l);
   problem_data.geo_name = geometry;
   problem_data.drchlt_sides = [1 3];
   % enforced using Lagrange multipliers
   problem_data.nmnn_sides   = [3];

   % material parameters
   problem_data.b   = 0.008;
   problem_data.rho = 9250;
   problem_data.E   = 25e9;
   problem_data.I   = problem_data.b^3/12;
   problem_data.A   = problem_data.b^2;
   problem_data.B   = 1;

   problem_data.g_mag = @(x,y,ib) zeros(size(x));
   problem_data.h_mag = @(x,y,ib) zeros(size(x));

   method_data.degree     = [2 3];
   % degree-1
   method_data.regularity = method_data.degree-1;
   % to be determined by convergence study
   method_data.nsub       = [1 1];
   % degree+1
   method_data.nquad      = method_data.degree+1;
end

function [problem_data, method_data] = setup_problem1d ()
   % geometry
   l   = 0.3;
   geometry = nrbline ([0 0],[l 0]);
   problem_data.geo_name = geometry;

   % material parameters
   problem_data.b   = 0.008;
   problem_data.rho = 9250;
   problem_data.E   = 25e9;
   problem_data.I   = problem_data.b^3/12;
   problem_data.A   = problem_data.b^2;
   problem_data.B   = 1;

   % boundary conditions
   problem_data.drchlt1_ends = [true false];
   problem_data.drchlt2_ends = [true false];
   problem_data.nmnn1_ends   = [false false];
   problem_data.nmnn2_ends   = [true false];

   problem_data.f = @(x) zeros(size(x));

   method_data.degree     = 3;
   % degree-1
   method_data.regularity = method_data.degree-1;
   % to be determined by convergence study
   method_data.nsub       = 32;
   % degree+1
   method_data.nquad      = method_data.degree+1;
end

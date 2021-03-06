function [problem_data, method_data] = setup_problem1d ()
   % geometry
   l = 0.3;
   b = 0.008;
   geometry = nrbline ([0 0], [l 0]);
   problem_data.geo_name = geometry;
   problem_data.drchlt1_sides = [1];
   problem_data.drchlt2_sides = [1];
   % secondary boundary condition (rather use secondary Dirichlet)
   problem_data.nmnn2_sides   = [1];

   % material parameters
   problem_data.l   = l;
   problem_data.b   = b;
   problem_data.rho = 9250;
   problem_data.E   = 25e9;
   problem_data.I   = problem_data.b^3/12;
   problem_data.A   = problem_data.b^2;
   % coupling coefficient?
   problem_data.D   = 298e3;
   problem_data.B   = 1;

   problem_data.g = @(ib) 0;

   method_data.degree     = [2 3];
   % degree-1
   method_data.regularity = method_data.degree-1;
   method_data.nsub       = [1 2^0];
   % degree+1
   method_data.nquad      = method_data.degree+1;
end

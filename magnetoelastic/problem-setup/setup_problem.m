function [problem_data, method_data] = setup_problem (geometry_file)
   problem_data.geo_name     = [geometry_file '.txt'];
   problem_data.nmnn_sides   = [];
   problem_data.drchlt_sides = [1];

   % add E, nu, G for mechanics (differentiate between nu)

   % diagonal components of reluctivity tensor in order, as cell array of function handles
   nu11 = 1/3.77e-6;
   nu22 = 1/1.012e-6;
   problem_data.nu = {@(x,y,iptc) compute_nu(x, y, iptc, nu11), @(x,y,iptc) compute_nu(x, y, iptc, nu22)};

   % coils defined via rectangles with homogeneous current density
   coils.bll = [0.5 3.5];
   coils.bur = [2.5 4];
   coils.tll = [0.5 4.5];
   coils.tur = [2.5 5];
   % 100 windings with 5 Ampere each
   coils.current = 100*5;
   problem_data.f = @(x,y,iptc) compute_source(x, y, iptc, coils);

   problem_data.g = @(x,y,ib) zeros(size(x));
   problem_data.h = @(x,y,ib) zeros(size(x));

   method_data.degree     = [2 2];
   % degree-1
   method_data.regularity = method_data.degree - 1;
   % to be determined by convergence study
   method_data.nsub       = [32 32];
   % degree+1
   method_data.nquad      = method_data.degree + 1;
end

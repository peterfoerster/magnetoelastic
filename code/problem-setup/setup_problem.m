function [problem_data, method_data] = setup_problem (geometry_file)
   % magnetic problem
   problem_data.geo_mag = [geometry_file '_mag.txt'];
   problem_data.nmnn_sides_mag   = [];
   problem_data.drchlt_sides_mag = [1];

   % diagonal components of permeability tensor, in order, as cell array of function handles
   mu11 = 3.77e-6;
   mu22 = 1.012e-5;
   problem_data.mu = {@(x,y,iptc) compute_mu(x, y, iptc, mu11), @(x,y,iptc) compute_mu(x, y, iptc, mu22)};

   % coils defined via rectangles with homogeneous current density
   % v1
   coils.bll = [0.5 2];
   coils.bur = [2.5 2.5];
   coils.tll = [0.5 4.5];
   coils.tur = [2.5 5];

   % 100 windings with 5 Ampere each
   coils.current = 100*5;
   problem_data.f_mag = @(x,y,iptc) compute_f_mag(x, y, iptc, coils);
   problem_data.g_mag = @(x,y,ib) zeros(size(x));
   problem_data.h_mag = @(x,y,ib) zeros(size(x));

   % mechanic problem
   problem_data.geo_mec = [geometry_file '_mec.txt'];
   problem_data.nmnn_sides_mec   = [2 3];
   problem_data.drchlt_sides_mec = [1];

   % material parameters
   E1   = 2.6e10;
   E2   = 2.27e10;
   nu12 = 0.429;
   G12  = 1/11*1e11;
   tau  = 8.67e6;

   problem_data.E  = {@(x,y,iptc) compute_E(x, y, iptc, E1), @(x,y,iptc) compute_E(x, y, iptc, E2)};
   problem_data.nu = {@(x,y,iptc) compute_nu(x, y, iptc, nu12)};
   problem_data.G  = {@(x,y,iptc) compute_G(x, y, iptc, G12)};

   problem_data.f_mec = @(x,y,iptc) zeros(2, size(x,1), size(x,2));
   problem_data.g_mec = @(x,y,ib) compute_tau(x, y, ib, tau);
   problem_data.h_mec = @(x,y,ib) zeros(2, size(x,1), size(x,2));

   % coupling parameters
   e11 = 213.3;
   e21 = -17.66;
   e62 = 150;

   problem_data.f = {@(x,y,iptc) compute_f(x, y, iptc, e11/mu11), @(x,y,iptc) compute_f(x, y, iptc, e21/mu11), ...
                     @(x,y,iptc) compute_f(x, y, iptc, e62/mu22)};

   method_data.degree     = [2 2];
   % degree-1
   method_data.regularity = method_data.degree-1;
   % to be determined by convergence study
   method_data.nsub       = [128 128];
   % degree+1
   method_data.nquad      = method_data.degree+1;
end

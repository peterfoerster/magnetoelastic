function [problem_data, method_data] = setup_problem (geometry_mag, geometry_mech)
   % magnetic problem
   problem_data.geo_mag = [geometry_mag '.txt'];
   problem_data.nmnn_sides_mag   = [];
   problem_data.drchlt_sides_mag = [1];

   % diagonal components of reluctivity tensor in order, as cell array of function handles
   nu_mag11 = 1/3.77e-6;
   nu_mag22 = 1/1.012e-6;
   problem_data.nu_mag = {@(x,y,iptc) compute_nu_mag(x, y, iptc, nu_mag11), @(x,y,iptc) compute_nu_mag(x, y, iptc, nu_mag22)};

   % coils defined via rectangles with homogeneous current density
   coils.bll = [0.5 3.5];
   coils.bur = [2.5 4];
   coils.tll = [0.5 4.5];
   coils.tur = [2.5 5];
   % 100 windings with 5 Ampere each
   coils.current = 100*5;
   problem_data.f_mag = @(x,y,iptc) compute_source(x, y, iptc, coils);
   problem_data.g_mag = @(x,y,ib) zeros(size(x));
   problem_data.h_mag = @(x,y,ib) zeros(size(x));

   % mechanic problem (second identical geometry with different boundaries)
   problem_data.geo_mech = [geometry_mech '.txt'];
   problem_data.nmnn_sides_mech   = [2];
   problem_data.drchlt_sides_mech = [1];

   % material coefficients
   E1        = 2.6e10;
   E2        = 2.27e10;
   nu_mech12 = 0.429;
   G12       = 1/11*1e11;
   problem_data.E         = {@(x,y,iptc) compute_E(x, y, iptc, E1), @(x,y,iptc) compute_E(x, y, iptc, E2)};
   problem_data.nu_mech12 = {@(x,y,iptc) compute_nu_mech(x, y, iptc, nu_mech12)};
   problem_data.G12       = {@(x,y,iptc) compute_G(x, y, iptc, )}

   problem_data.f_mech = @(x,y,iptc) compute_source(x, y, iptc, coils);
   problem_data.g_mech = @(x,y,ib) zeros(size(x));
   problem_data.h_mech = @(x,y,ib) zeros(size(x));

   method_data.degree     = [2 2];
   % degree-1
   method_data.regularity = method_data.degree - 1;
   % to be determined by convergence study
   method_data.nsub       = [8 8];
   % degree+1
   method_data.nquad      = method_data.degree + 1;
end

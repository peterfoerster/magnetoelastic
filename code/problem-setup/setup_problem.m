function [problem_data, method_data] = setup_problem (geometry_file)
   % magnetic problem
   problem_data.geo_name = [geometry_file '.txt'];
   problem_data.nmnn_sides_mag   = [];
   problem_data.drchlt_sides_mag = [1];

   % diagonal components of reluctivity tensor in order, as cell array of function handles
   mu_mag11 = 3.77e-6;
   mu_mag22 = 1.012e-6;
   problem_data.nu_mag = {@(x,y,iptc) compute_nu_mag(x, y, iptc, 1/mu_mag11), @(x,y,iptc) compute_nu_mag(x, y, iptc, 1/mu_mag22)};

   % coils defined via rectangles with homogeneous current density
   coils.bll = [0.5 3.5];
   coils.bur = [2.5 4];
   coils.tll = [0.5 4.5];
   coils.tur = [2.5 5];
   % 100 windings with 5 Ampere each
   coils.current = 100*5;
   problem_data.f_mag = @(x,y,iptc) compute_f_mag(x, y, iptc, coils);
   problem_data.g_mag = @(x,y,ib) zeros(size(x));
   problem_data.h_mag = @(x,y,ib) zeros(size(x));

   % mechanic problem
   problem_data.iptc = 5;
   problem_data.nmnn_sides_mech   = [2 3 4];
   problem_data.drchlt_sides_mech = [1];

   % material parameters
   % E1        = 2.6e10;
   % E2        = 2.27e10;
   % nu_mech12 = 0.429;
   % G12       = 1/11*1e11;
   % somehow given by the magnetic problem?
   tau       = 8.67e6;

   % test case
   E1 = E2   = 1e10;
   nu_mech12 = 0.4;
   G12       = E1/(2*(1+nu_mech12));
   problem_data.lambda_lame = @(x,y) (nu_mech12*E1)/((1+nu_mech12)*(1-2*nu_mech12))*ones(size(x));
   problem_data.mu_lame     = @(x,y) G12*ones(size(x));

   problem_data.E       = {@(x,y,iptc) compute_E(x, y, iptc, E1), @(x,y,iptc) compute_E(x, y, iptc, E2)};
   problem_data.nu_mech = {@(x,y,iptc) compute_nu_mech(x, y, iptc, nu_mech12)};
   problem_data.G       = {@(x,y,iptc) compute_G(x, y, iptc, G12)};

   problem_data.f_mech = @(x,y,iptc) [zeros(size(x)), zeros(size(x))];
   problem_data.g_mech = @(x,y,ib) compute_tau(x, y, ib, tau);
   problem_data.h_mech = @(x,y,ib) [zeros(size(x)); zeros(size(x))];

   % coupling
   e11 = 213.3;
   e21 = -17.66;
   e62 = 150;
   problem_data.fc = {@(x,y,iptc) compute_fc(x, y, iptc, f11/mu_mag11), @(x,y,iptc) compute_fc(x, y, iptc, f21/mu_mag11), @(x,y,iptc) compute_fc(x, y, iptc, f62/mu_mag22)};

   method_data.degree     = [2 2];
   % degree-1
   method_data.regularity = method_data.degree - 1;
   % to be determined by convergence study
   method_data.nsub       = [16 16];
   % degree+1
   method_data.nquad      = method_data.degree + 1;
end

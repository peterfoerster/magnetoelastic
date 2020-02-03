% input:
% degree_ref, nsub_ref, nquad_offset_ref, degree, N_it, nquad_offset
% output:
% folder with files for error per element
% file with norms of values for convergence study

function [] = conv_study_err (degree_ref, nsub_ref, nquad_offset_ref, degree, N_it, nquad_offset)
geometry_file = 'photocathode_insulator';
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load ([geometry_file '.txt']);

normh1_ref = NaN(1,N_it+1);
normh1 = NaN(1,N_it+1);
norml2_ref = NaN(1,N_it+1);
norml2 = NaN(1,N_it+1);

degree_ref = [degree_ref degree_ref];
nsub_ref = [nsub_ref nsub_ref];
degree = [degree degree];
regularity = degree_ref-1;
nquad = degree_ref+1+nquad_offset_ref;
npatch = numel (geometry);
msh_ref = cell (1, npatch);
sp_ref = cell (1, npatch);

filename = [geometry_file '_degree=' num2str(degree_ref(1)) '_nsub=' num2str(nsub_ref(1)) '_nquad_offset=' num2str(nquad_offset_ref) '.mat'];
load (filename);
u_ref = u;

% recreate msh and space
for iptc=1:npatch
  [knots{iptc}, zeta{iptc}] = kntrefine (geometry(iptc).nurbs.knots, nsub_ref-1, degree_ref, regularity);
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh_ref{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));
  sp_ref{iptc} = sp_bspline (knots{iptc}, degree_ref, msh_ref{iptc});
end
msh_ref = msh_multipatch (msh_ref, boundaries);
space_ref = sp_multipatch (sp_ref, msh_ref, interfaces, boundary_interfaces);
clear sp_ref;

for iit=0:N_it
  fprintf('\niteration with nsub = %d\n', 2^iit);
  nsub = [2^iit 2^iit];
  filename = [geometry_file '_degree=' num2str(degree(1)) '_nsub=' num2str(2^iit) '_nquad_offset=' num2str(nquad_offset) '.mat'];
  load(filename);

  regularity = degree-1;
  nquad = degree+1+nquad_offset;
  msh = cell (1, npatch);
  sp = cell (1, npatch);

  % recreate msh and space
  for iptc=1:npatch
    [knots{iptc}, zeta{iptc}] = kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
    rule      = msh_gauss_nodes (nquad);
    [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
    msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));
    sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
  end
  msh = msh_multipatch (msh, boundaries);
  space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
  clear msh;
  clear sp;

  % compute norms and create error elem files simultaneously
  tic;
  filename = [geometry_file '_degree_ref=' num2str(degree_ref(1)) '_nsub_ref=' num2str(nsub_ref(1)) '_nquad_offset_ref=' num2str(nquad_offset_ref(1)) '_degree=' num2str(degree(1)) '_nsub=' num2str(2^iit) '_nquad_offset=' num2str(nquad_offset)];
  mkdir(filename);
  [normh1_ref(iit+1), normh1(iit+1), norml2_ref(iit+1), norml2(iit+1)] = mp_sp_h1_norm (space_ref, msh_ref, u_ref, u, space, geometry, filename);
  fprintf('time elapsed for norm computation:%d\n', toc);
end

% save norms and info
filename = [geometry_file '_degree_ref=' num2str(degree_ref(1)) '_nsub_ref=' num2str(nsub_ref(1)) '_nquad_offset_ref=' num2str(nquad_offset_ref(1)) '_degree=' num2str(degree(1)) '_N_it=' num2str(N_it) '_nquad=' num2str(nquad_offset) '.mat'];
save(filename, 'normh1_ref', 'normh1', 'norml2_ref', 'norml2');
end

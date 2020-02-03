% single patch

function [normh1_ref, normh1, norml2_ref, norml2] = sp_h1_norm (msh_ref, space_ref, u_ref, u, space, geometry, iptc, filename)
norml2_ref = 0; norml2 = 0; normh1s_ref = 0; normh1s = 0; errh1 = [];

% compute points for vtk
[~, F] = sp_eval (u_ref, space_ref, geometry, msh_ref.breaks);

for iel = 1:msh_ref.nel_dir(1)
  msh_col = msh_evaluate_col (msh_ref, iel);
  sp_col = sp_evaluate_col (space_ref, msh_col, 'value', true, 'gradient', true);
  % evaluate h1
  [l2_ref, l2, h1s_ref, h1s, errh1_elem] = sp_h1_norm_eval_h1 (sp_col, msh_col, u_ref, space_ref, u, space, geometry);

  % order errors for vtk
  errh1 = [errh1; errh1_elem];

  norml2_ref = norml2_ref + l2_ref.^2;
  norml2 = norml2 + l2.^2;

  normh1s_ref = normh1s_ref + h1s_ref.^2;
  normh1s = normh1s + h1s.^2;
end

% export to vtk
msh_to_vtk_cell (F, errh1, [filename '/' filename '_' num2str(iptc) '.vts'], 'err_h1');

normh1_ref = sqrt (norml2_ref + normh1s_ref);
normh1 = sqrt (norml2 + normh1s);

norml2_ref = sqrt (norml2_ref);
norml2 = sqrt (norml2);
end

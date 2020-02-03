% evaluate l2

function [l2_ref, l2, errl2_elem] = sp_l2_norm_eval_l2 (sp, msh, u_ref, space_ref, u, space, geometry)
% get quadrature points
msh_pts_x = msh.qn{1};
msh_pts_y = msh.qn{2};
msh_pts_y = reshape( msh_pts_y, 1, numel(msh_pts_y) );

valu_ref = sp_eval(u_ref, space_ref, geometry, {msh_pts_x, msh_pts_y}, 'value');
valu = sp_eval(u, space, geometry, {msh_pts_x, msh_pts_y}, 'value');

% reshape according to original
valu_ref = reshape (valu_ref, sp.ncomp, msh.nqn, msh.nel);
valu = reshape (valu, sp.ncomp, msh.nqn, msh.nel);

% compute determinant and quadrature weights
w = msh.quad_weights .* msh.jacdet;

errl2_elem = sum (reshape (sum ((valu - valu_ref).^2, 1), [msh.nqn, msh.nel]) .* w);
errl2_elem  = sqrt (errl2_elem);

l2_elem = sum (reshape (sum ((valu).^2, 1), [msh.nqn, msh.nel]) .* w);
l2  = sqrt (sum (l2_elem));
l2_elem_ref = sum (reshape (sum ((valu_ref).^2, 1), [msh.nqn, msh.nel]) .* w);
l2_ref  = sqrt (sum (l2_elem_ref));
end

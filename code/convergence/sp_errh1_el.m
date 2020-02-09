function [errh1, errl2] = sp_errh1_el (geometry, msh_ref, sp_ref, u_ref, sp, u)
   % evaluate reference
   eu_ref = sp_eval (u_ref, sp_ref, geometry, msh_ref.qn, {'value', 'gradient'});
   valu_ref  = reshape (eu_ref{1}, sp_ref.ncomp, msh_ref.nqn, msh_ref.nel);
   gradu_ref = reshape (eu_ref{2}, msh_ref.rdim, msh_ref.nqn, msh_ref.nel);
   gradu_ref = reshape (gradu_ref, sp_ref.ncomp, msh_ref.rdim, msh_ref.nqn, msh_ref.nel);

   % evaluate iterative solution
   eu = sp_eval (u, sp, geometry, msh_ref.qn, {'value', 'gradient'});
   valu  = reshape (eu{1}, sp.ncomp, msh_ref.nqn, msh_ref.nel);
   gradu = reshape (eu{2}, msh_ref.rdim, msh_ref.nqn, msh_ref.nel);
   gradu = reshape (gradu, sp.ncomp, msh_ref.rdim, msh_ref.nqn, msh_ref.nel);

   jacdet_weights = msh_ref.jacdet .* msh_ref.quad_weights;

   % change to relative error?
   tmpl2 = (valu - valu_ref);
   tmpl2 = sum(tmpl2.^2, 1);
   tmpl2 = jacdet_weights .* reshape(tmpl2, [msh_ref.nqn, msh_ref.nel]);

   tmph1 = (gradu - gradu_ref);
   tmph1 = sum(sum(tmph1.^2, 1), 2);
   tmph1 = jacdet_weights .* reshape(tmph1, [msh_ref.nqn, msh_ref.nel]);

   errl2 = sum(sum(tmpl2, 1), 2);
   errh1 = sum(sum(tmpl2 + tmph1, 1), 2);
end

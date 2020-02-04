% ||u-u_ref||_H^1 = (S |u-u_ref|^2 + |grad(u)-grad(u_ref)|^2 dx)^1/2 
function [errh1 errl2] = mp_errh1 (msh_ref, space_ref, u_ref, u, space, geometry)
   if (space_ref.npatch ~= msh_ref.npatch || space.npatch ~= msh_ref.npatch)
      error('The number of patches does not coincide');
   end

   for iptc = 1:msh_ref.npatch
      % added part for reference solution
      if (isempty(space_ref.dofs_ornt))
         u_ref_ptc = u_ref(space_ref.gnum{iptc});
      else
         u_ref_ptc = u_ref(space_ref.gnum{iptc}) .* space_ref.dofs_ornt{iptc}.';
      end
      if (isempty(space.dofs_ornt))
         u_ptc = u(space.gnum{iptc});
      else
         u_ptc = u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.';
      end
      [normh1_ref, normh1, norml2_ref, norml2] = sp_h1_norm (msh_ref.msh_patch{iptc}, space_ref.sp_patch{iptc}, u_ptc_ref, u_ptc, space.sp_patch{iptc}, geometry(iptc), iptc, filename);
   end
% continue here
normh1_ref = sqrt (sum (normh1_ref .* normh1_ref));
normh1 = sqrt (sum (normh1 .* normh1));

norml2_ref = sqrt (sum (norml2_ref.* norml2_ref));
norml2 = sqrt (sum (norml2 .* norml2));
end

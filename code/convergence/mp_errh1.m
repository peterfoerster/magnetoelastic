% ||(u-u_ref)/u_ref||_H^1 = (S |(u-u_ref)/u_ref|^2 + |(grad(u)-grad(u_ref))/grad(u_ref)|^2 dx)^1/2

function [errh1, errl2] = mp_errh1 (geometry, msh_ref, space_ref, u_ref, space, u)
   if (space_ref.npatch ~= msh_ref.npatch || space.npatch ~= msh_ref.npatch)
      error('The number of patches does not coincide.');
   end

   errh1 = errl2 = 0;
   for iptc=1:msh_ref.npatch
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
% keyboard
      [errh1_ptc, errl2_ptc] = sp_errh1 (geometry(iptc), msh_ref.msh_patch{iptc}, ...
                                         space_ref.sp_patch{iptc}, u_ref_ptc, space.sp_patch{iptc}, u_ptc);
      errh1 = errh1 + errh1_ptc;
      errl2 = errl2 + errl2_ptc;
   end
   errh1 = sqrt(errh1);
   errl2 = sqrt(errl2);
end

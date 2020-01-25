function [f] = compute_source (x, y, iptc, coils)
   % compute homogeneous current densities
   Ab_coil = (coils.bur(1)-coils.bll(1))*(coils.bur(2)-coils.bll(2));
   jb = - coils.current/Ab_coil;
   At_coil = (coils.tur(1)-coils.tll(1))*(coils.tur(2)-coils.tll(2));
   jt = coils.current/At_coil;

   % evaluate current densities depending on the patch
   switch (iptc)
      case {4}
         % upper coil
         xl = x > coils.tll(1);
         xr = x < coils.tur(1);
         yb = y > coils.tll(2);
         yt = y < coils.tur(2);
         idx_coil = (xl & xr & yb & yt);
         ft = zeros(size(x)) + jt*idx_coil;
         % lower coil
         xl = x > coils.bll(1);
         xr = x < coils.bur(1);
         yb = y > coils.bll(2);
         yt = y < coils.bur(2);
         idx_coil = (xl & xr & yb & yt);
         fb = zeros(size(x)) + jb*idx_coil;
         f = fb + ft;
      otherwise
         f = zeros(size(x));
   end%switch
end

function [tau_iptc] = compute_tau (x, y, ib, tau)
   switch (ib)
      % inhomogeneous Neumann
      case {2}
         % plate
         % tau_iptc = [tau * ones(size(x)); zeros(size(x))];
         tau_iptc = [zeros(size(x)); tau * ones(size(x))];
      otherwise
         tau_iptc = [zeros(size(x)); zeros(size(x))];
   end%switch
end

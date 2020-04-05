function [tau_iptc] = compute_tau (x, y, ib, tau)
   switch (ib)
      % inhomogeneous Neumann
      case {2}
         % plate
         tau_iptc = [reshape(tau * ones(size(x)), [1, size(x)]); zeros(1, size(x,1), size(x,2))];
      otherwise
         tau_iptc = zeros(2, size(x,1), size(x,2));
   end%switch
end

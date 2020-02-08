function [domain] = discretize_domain_v2 (plate, outer_boundary)
   % right
   p1 = nrbeval(plate(1), 0);
   p2 = nrbeval(outer_boundary(1), 1);
   domain(1) = nrbline(p1, p2);

   p1 = nrbeval(plate(1), 1);
   p2 = nrbeval(outer_boundary(2), 1);
   domain(2) = nrbline(p1, p2);

   % top
   p1 = nrbeval(plate(2), 1);
   p2 = nrbeval(outer_boundary(4), 0);
   domain(3) = nrbline(p1, p2);

   p1 = nrbeval(plate(2), 0);
   p2 = nrbeval(outer_boundary(6), 1);
   domain(4) = nrbline(p1, p2);

   % left
   p1 = nrbeval(outer_boundary(8), 1);
   p2 = nrbeval(plate(3), 1);
   domain(5) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(8), 0);
   p2 = nrbeval(plate(3), 0);
   domain(6) = nrbline(p1, p2);

   % bottom
   p1 = nrbeval(outer_boundary(10), 1);
   p2 = nrbeval(plate(4), 0);
   domain(7) = nrbline(p1, p2);

   p1 = nrbeval(outer_boundary(11), 1);
   p2 = nrbeval(plate(4), 1);
   domain(8) = nrbline(p1, p2);
end

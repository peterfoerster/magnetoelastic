function [ptcs] = create_ptcs_v2 (plate, outer_boundary, domain);
   % right
   ptcs(1) = nrbcoons(outer_boundary(12), domain(1), domain(8), outer_boundary(1));
   ptcs(2) = nrbcoons(domain(1), domain(2), plate(1), outer_boundary(2));
   ptcs(3) = nrbcoons(domain(2), outer_boundary(4), domain(3), outer_boundary(3));

   % center
   ptcs(4) = nrbcoons(plate(2), outer_boundary(5), domain(4), domain(3));
   ptcs(5) = nrbcoons(plate(4), plate(2), plate(3), plate(1));
   ptcs(6) = nrbcoons(outer_boundary(11), plate(4), domain(7), domain(8));

   % left
   ptcs(7) = nrbcoons(outer_boundary(10), domain(6), outer_boundary(9), domain(7));
   ptcs(8) = nrbcoons(domain(6), domain(5), outer_boundary(8), plate(3));
   ptcs(9) = nrbcoons(domain(5), outer_boundary(6), outer_boundary(7), domain(4));
end

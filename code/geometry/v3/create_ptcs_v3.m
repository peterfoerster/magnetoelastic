function [ptcs] = create_ptcs_v3 (plate, outer_boundary, domain);
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

   % refinement (elevate degree to retain second derivatives)
   for iptc=1:length(ptcs)
      ptcs(iptc) = nrbdegelev(ptcs(iptc), [1 1]);
   end

   ptcs(1) = nrbkntins(ptcs(1), {[1/8 1/4], [3/4 7/8]});
   ptcs(2) = nrbkntins(ptcs(2), {[1/8 1/4], [1/8 1/4 3/4 7/8]});
   ptcs(3) = nrbkntins(ptcs(3), {[1/8 1/4], [1/8 1/4]});
   ptcs(4) = nrbkntins(ptcs(4), {[1/8 1/4 3/4 7/8], [1/8 1/4]});
   ptcs(5) = nrbkntins(ptcs(5), {[1/8 1/4 3/4 7/8], [1/8 1/4 3/4 7/8]});
   ptcs(6) = nrbkntins(ptcs(6), {[1/8 1/4 3/4 7/8], [3/4 7/8]});
   ptcs(7) = nrbkntins(ptcs(7), {[3/4 7/8], [3/4 7/8]});
   ptcs(8) = nrbkntins(ptcs(8), {[3/4 7/8], [1/8 1/4 3/4 7/8]});
   ptcs(9) = nrbkntins(ptcs(9), {[3/4 7/8], [1/8 1/4]});
end

function [] = write_geometryfile (ptcs, filename)
   [interfaces, boundaries] = nrbmultipatch (ptcs);
   boundaries = [];

   boundaries(1).patches = [1 1 2 3 3 4 4 5 6 6];
   boundaries(1).faces   = [2 3 2 2 4 1 4 1 1 3];

   for ibnd = 1:length(boundaries)
      boundaries(ibnd).nsides = length(boundaries(ibnd).patches);
   end

   nrbexport (ptcs, interfaces, boundaries, filename);
end

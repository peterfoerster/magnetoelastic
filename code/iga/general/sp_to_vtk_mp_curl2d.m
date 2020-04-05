% INPUT:
%
%     A:           vector of dof weights
%     space:       object representing the space of discrete functions (see sp_multipatch)
%     geometry:    geometry structure (see geo_load)
%     npts:        number of points along each parametric direction where to evaluate
%     filename:    name of the output file.
%     fieldnames:  how to name the saved variables in the vtk file
%
% OUTPUT:
%
%    none

function sp_to_vtk_mp_curl2d (A, space, geometry, npts, filename, fieldname)

  str1 = cat (2,'<?xml version="1.0"?> \n', ...
'<VTKFile type="Collection" version="0.1"> \n', ...
'<Collection> \n');

  str2 = cat (2, '<DataSet part="%d" file="%s.vts"/> \n');

  str3 = cat (2, ...
'</Collection>\n', ...
'</VTKFile> \n');

  if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
    pvd_filename = cat (2, filename, '.pvd');
  else
    pvd_filename = filename;
    filename = filename (1:end-4);
  end

  fid = fopen (pvd_filename, 'w');
  if (fid < 0)
    error ('sp_to_vtk_mp_curl2d: could not open file %s', pvd_filename);
  end

  fprintf (fid, str1);
  ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
  if (isempty (ind)); ind = 1; end
  for iptc = 1:space.npatch
    filename_patch_without_path = cat (2, filename(ind:end), '_', num2str (iptc));
    filename_patch = cat (2, filename, '_', num2str (iptc));
    fprintf (fid, str2, iptc, filename_patch_without_path);
    if (isempty (space.dofs_ornt))
      sp_to_vtk_curl2d (A(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc), npts, ...
                           filename_patch, fieldname)
    else
      sp_to_vtk (A(space.gnum{iptc}) .* space.dofs_ornt{iptc}', space.sp_patch{iptc}, geometry(iptc), npts, ...
                           filename_patch, fieldname)
    end
  end
  fprintf (fid, str3);

  fclose (fid);

end

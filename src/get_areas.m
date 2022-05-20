function [areas] = get_areas(sinfo)
  verts = sinfo.verts;
  [~,nv] = size(verts);
  triind = sinfo.triind;
  npatches = sinfo.npatches;
  
  areas = zeros(npatches,1);
  mex_id_ = 'get_areas(i int[x], i int[x], i double[xx], i int[xx], io double[x])';
[areas] = flattrirouts(mex_id_, npatches, nv, verts, triind, areas, 1, 1, 3, nv, 3, npatches, npatches);
end

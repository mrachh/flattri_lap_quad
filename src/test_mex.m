A = load('for_Manas.mat');
sinfo = [];
sinfo.centers = A.Center';
[~,sinfo.npatches] = size(sinfo.centers);
sinfo.verts = A.P';
sinfo.normals = A.normals';
sinfo.triind = A.t';

areas = get_areas(sinfo);
sinfo.area = areas;

nttest = 100;
targs_all = A.obsPoints';
targs = targs_all(:,1:nttest);

charges = A.c;

nover = 3;
rfac=2.75;
eps = 1e-4;
grad_comp = lpeval_lap_g(targs,sinfo,charges,rfac,nover,eps);

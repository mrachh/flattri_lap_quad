addpath '../src'
% A = load('CombinedMesh.mat');
% E = load('EC.mat');
% N = length(A.Area);

% Must be followed by run of bem2_charge_engine

F = get_factorization(Center,Area,normals,EC,contrast,weight);


tic, Xsol = solve_fds(F,b,Area); tsolve = toc;

fprintf('solve time=%d\n',tsolve);


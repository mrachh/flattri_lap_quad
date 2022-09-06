function F = get_factorization(center,area,normals,qcorr,contrast,weight,opts)
%
%  This subroutine constructs a recursive skeletonization with strong
%  admissibility factorization for a laplace TE polarization problem
%  the integral equation being solved is
%  (1 + \kappa 2S')\sigma = f 
%  where \kappa is the constrast (note the flip in sign is somehow
%  compensated for the difference between fields and gradients)
%
%  The surface discretization is a first order triangulation, and
%  quadrature corrections are already scaled by the contrast
%
%  Input arguments:
%  - center(N,3): centroid of the triangles
%  - area(N,1): surface area of each triangle
%  - normals(N,3): face normals of the triangles
%  - qcorr (sparse matrix) (NxN): contrast scaled quadrature corrections, 
%                same as EC
%  - contrast(N,1): contrast at each points
%  - weight (Scalar): weight for charge conservation scaling
%  - opts: optional argument struct
%        opts.occ - set max occupancy for oct tree
%        opts.rank_or_tol - precision for constructing factorization
%
%  Output arguments:
%  - F: factorization which can be used for solving the integral equation
%      using the solve_fds routine
%  
%  Notes:
%    This subroutine implicitly scales by the square root of the weights
%    for stability reasons. Do not directly use srskelf_mv /sv directly
%    as it wouldn't work, use solve_fds instead. 
%    Both the data, and the solution need to be
%    rescaled by the square root of weights
%
%
    if(nargin == 6)
        opts = [];
    end
    occ = 256;
    rank_or_tol = 1e-4;
    if(isfield(opts,'occ'))
        occ = opts.occ;
    end
    if(isfield(opts,'rank_or_tol'))
        rank_or_tol = opts.rank_or_tol;
    end
    x = center.';
    area = area.';
    nu = normals.';
    contrast = contrast(:);
    [~,N] = size(x);
    P = zeros(N,1);
    
    rsurf = sum(area);
    wuse = weight/rsurf;
    
    
    
    zk = 0.1;
    opts_use = struct('verb',1,'symm','n','zk',zk);
    opts_use.lap_proxy = true;
    Afun_use = @(i,j) Afun_lap_te(i,j,x,nu,area,P,qcorr,contrast,wuse);
    pxyfun_use = @(x,slf,nbr,proxy,l,ctr) pxyfun_lap_neumann(x,slf,nbr,proxy,l,ctr,area);
    tic, F = srskelf_asym_new(Afun_use,x,occ,rank_or_tol,pxyfun_use,opts_use); tfac = toc;

end
function F = get_factorization(center,area,normals,qcorr,contrast,weight,opts)

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
    opts = struct('verb',1,'symm','n','zk',zk);
    Afun_use = @(i,j) Afun_lap_te(i,j,x,nu,area,P,qcorr,contrast,wuse);
    pxyfun_use = @(x,slf,nbr,proxy,l,ctr) pxyfun_lap_neumann(x,slf,nbr,proxy,l,ctr,area);
    tic, F = srskelf_asym_new(Afun_use,x,occ,rank_or_tol,pxyfun_use,opts); tfac = toc;

end
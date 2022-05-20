      implicit real *8 (a-h,o-z)
      real *8, allocatable :: centers(:,:),normals(:,:),charges(:)
      real *8, allocatable :: chargesuse(:)
      real *8, allocatable :: verts(:,:),targs(:,:)
      integer, allocatable :: triind(:,:)

      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8, allocatable :: cms(:,:),rads(:)
      integer, allocatable :: norders(:),iptype(:),ixyzs(:)

      real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: uvs_targ(:,:),gradrho(:,:)
      integer, allocatable :: ipatch_id(:)
      real *8 p0(3),p1(3),p2(3)

      call prini(6,13)
      open(unit=33,file='brain_all.dat',form='unformatted')
      read(33) ntri
      read(33) nvert
      read(33) ntarg


      allocate(centers(3,ntri),normals(3,ntri),charges(ntri))
      allocate(verts(3,nvert),targs(3,ntarg),triind(3,ntri))

      read(33) centers
      read(33) normals
      read(33) charges
      read(33) verts
      read(33) targs
      read(33) triind

      close(33)


      npatches = ntri
      norder = 1
      npols = (norder+1)*(norder+2)/2
      npts = npatches*npols

      print *, npatches,npts,norder,npols
      allocate(norders(npatches),iptype(npatches),ixyzs(npatches+1))
      allocate(srcvals(12,npts),srccoefs(9,npts))

      allocate(chargesuse(npts))

      do i=1,npatches
        ixyzs(i) = (i-1)*npols + 1
        iptype(i) = 1
        norders(i) = norder
      enddo
      ixyzs(npatches+1) = npts+1

      allocate(uvs(2,npols),umat(npols,npols),vmat(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(norder,npols,uvs,umat,vmat,wts)

      call prinf('triind=*',triind,20)
      call prin2('verts=*',verts,24)
      print *, "Here"
      do ipatch=1,npatches
        p0(1:3) = verts(1:3,triind(1,ipatch))
        p1(1:3) = verts(1:3,triind(2,ipatch))
        p2(1:3) = verts(1:3,triind(3,ipatch))
        do i=1,npols
          u = uvs(1,i)
          v = uvs(2,i)
          ipt = (ipatch-1)*npols + i
          srcvals(1:3,ipt) = p0(1:3) + u*(p1(1:3)-p0(1:3)) + &
              v*(p2(1:3)-p0(1:3)) 
          srcvals(4:6,ipt) = p1(1:3)-p0(1:3)
          srcvals(7:9,ipt) = p2(1:3)-p0(1:3)
          call cross_prod3d(srcvals(4,ipt),srcvals(7,ipt), &
            srcvals(10,ipt))
          ds = sqrt(srcvals(10,ipt)**2 + srcvals(11,ipt)**2 + &
             srcvals(12,ipt)**2)
          srcvals(10:12,ipt) = srcvals(10:12,ipt)/ds
          chargesuse(ipt) = charges(ipatch)
        enddo

        do i=1,npols
          ipt = (ipatch-1)*npols + i
          do j=1,9
            srccoefs(j,ipt) = 0
            do l=1,npols
              lpt = (ipatch-1)*npols + l
              srccoefs(j,ipt) = srccoefs(j,ipt) + &
                 umat(i,l)*srcvals(j,lpt)
            enddo
          enddo
        enddo
      enddo

      call prin2('normals1=*',srcvals(10:12,1:npts),24)
      call prin2('normals2=*',normals,24)

      allocate(cms(3,npatches),rads(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

      

      
      
      erra = 0

      errc = 0
      do i=1,npatches
        ipt = ixyzs(i)
        erra = erra + (srcvals(10,ipt) - normals(1,i))**2
        erra = erra + (srcvals(11,ipt) - normals(2,i))**2
        erra = erra + (srcvals(12,ipt) - normals(3,i))**2

        errc = errc + (cms(1,i)-centers(1,i))**2
        errc = errc + (cms(2,i)-centers(2,i))**2
        errc = errc + (cms(3,i)-centers(3,i))**2
      enddo
      erra = sqrt(erra/(npatches+0.0d0))
      errc = sqrt(errc/(npatches+0.0d0))
      call prin2('error in normals=*',erra,1)
      call prin2('error in centroids=*',errc,1)

      allocate(uvs_targ(2,ntarg),ipatch_id(ntarg))
      ipatch_id = -1
      uvs_targ = 0
      ndtarg = 3
      eps = 1.0d-4
      eps_quad = 1.0d-4
      ntarguse = ntarg
      allocate(gradrho(3,ntarg))
!      npatches = 5
!      npts = npatches*npols
      call lpcomp_grads0(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarguse,targs,ipatch_id, &
        uvs_targ,eps_quad,eps,chargesuse,gradrho)
      open(unit=44,file='grad_rho.dat')
      do i=1,ntarg
        write(44,*) gradrho(1,i),gradrho(2,i),gradrho(3,i)
      enddo
      close(44)


      stop
      end
     

      subroutine getnearquad_magnetostatics(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad, &
        wnear)
!
!  This subroutine generates the near field quadrature
!  for the representations:
!
!  wnear: 
!   \nabla_{1} S_{0}, \nabla_{2} S_{0}, \nabla_{3} S_{0},
!
!  The quadrature is computed by the following strategy
!  targets within a sphere of radius rfac0*rs
!  of a patch centroid is handled using adaptive integration
!  where rs is the radius of the bounding sphere
!  for the patch
!  
!  All other targets in the near field are handled via
!  oversampled quadrature
!
!  The recommended parameter for rfac0 is 1.25d0
!  
!  NOTES:
!    - wnear must be of size nquad,3 as 3 different layer
!      potentials are returned
!      * the first kernel is \nabla_{x} S_{0}
!      * the second kernel is \nabla_{y} S_{0}
!      * the third kernel is \nabla_{z} S_{0}
! 
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id=-1, if target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of patch if target on surface,
!        otherwise not used
!    - eps: real *8
!        precision requested
!    - dpars: real *8 (1)
!        the parameter k for the yukawa kernels
!    - iquadtype: integer
!        quadrature type
!          * iquadtype = 1, use ggq for self + adaptive integration
!            for rest
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - rfac0: integer
!        radius parameter for near field
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!
!  Output arguments
!    - wnear: complex *16(3*nquad)
!        The desired near field quadrature
!               
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches),ixyzs(npatches+1)
      integer, intent(in) :: iptype(npatches),npts
      real *8, intent(in) :: srccoefs(9,npts),srcvals(12,npts)
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg),uvs_targ(2,ntarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: eps
      integer, intent(in) :: iquadtype,nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(in) :: rfac0
      integer, intent(in) :: nquad
      real *8, intent(out) :: wnear(nquad,3)

      integer ndz,ndd,ndi
      integer ipv
      real *8 dpars
      integer ipars
      complex *16 zpars

      procedure (), pointer :: fker

      external l3d_sgradx,l3d_sgrady,l3d_sgradz

      ndd = 0
      ndz = 0
      ndi = 0
      ipv = 1
      dpars = 0
      zpars = 0
      ipars = 0

      fker => l3d_sgradx
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,1))
!      print *, "done with kernel 1"

      fker => l3d_sgrady
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,2))
!      print *, "done with kernel 2"

      fker => l3d_sgradz
      ipv = 1
      call dgetnearquad_ggq_guru(npatches,norders,ixyzs,&
       iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,&
       ipatch_id,uvs_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,ndi,&
       ipars,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear(1,3))
!      print *, "done with kernel 3"



      return
      end subroutine getnearquad_magnetostatics 
!
!
!
!
!
!
      subroutine lpcomp_grads0(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps_quad,eps,rho,gradrho)

!
!  This subroutine evaluates the potential 
!  1. grad S_{0}[rho]
!
!
!  The identity term is not included in the representation if targets are on 
!  surface
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - ipatch_id: integer(ntarg)
!        id of patch of target i, id=-1, if target is off-surface
!    - uvs_targ: real *8 (2,ntarg)
!        local uv coordinates of patch if target on surface,
!        otherwise not used
!    - eps: real *8
!        precision requested
!    - rho: real *8 (npts)
!        The charge density rho
!
!  Output arguments:
!    - gradrho: real *8 (3,npts)
!         returns the potential grad S_{0}[\rho]
!
!
!
!
      implicit none
      integer, intent(in) :: npatches,norders(npatches)
      integer, intent(in) :: iptype(npatches),ixyzs(npatches+1)
      integer, intent(in) :: npts,ndtarg,ntarg
      real *8, intent(in) :: srcvals(12,npts),srccoefs(9,npts)
      real *8, intent(in) :: targs(ndtarg)
      integer, intent(in) :: ipatch_id(ntarg)
      real *8, intent(in) :: uvs_targ(2,ntarg)
      real *8, intent(in) :: eps
      real *8, intent(in) :: rho(npts)
      real *8, intent(out) :: gradrho(3,ntarg)

      integer nptso,nnz,nquad
      integer nover,npolso,npols,norder
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      real *8, allocatable :: wnear(:,:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyzso(:),novers(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer iptype_avg,norder_avg,iquadtype,npts_over,ikerorder
      integer i,noveruse,npover
      real *8 rfac,rfac0,t1,t2,omp_get_wtime,eps_quad
      complex *16 zpars

      iptype_avg = floor(sum(iptype)/(npatches+0.0d0))
      norder_avg = floor(sum(norders)/(npatches+0.0d0))

      call get_rfacs(norder_avg,iptype_avg,rfac,rfac0)
! Get centroid and bounding sphere info

      allocate(cms(3,npatches),rads(npatches),rad_near(npatches))

      call get_centroid_rads(npatches,norders,ixyzs,iptype,npts, &
          srccoefs,cms,rads)

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npatches
        rad_near(i) = rads(i)*rfac
      enddo
!$OMP END PARALLEL DO

!
!  Find near quadrature corrections
!
!      call prin2('rad_near=*',rad_near,24)
!      call prin2('rfac=*',rfac,1)
!      call prinf('ntarg=*',ntarg,1)
      print *, "Starting findnear mem"
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      
      call findnearmem(cms,npatches,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))

      
      call findnear(cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr, & 
             col_ind)


      allocate(iquad(nnz+1)) 
      call get_iquad_rsc(npatches,ixyzs,ntarg,nnz,row_ptr,col_ind, &
              iquad)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()      
      call prin2('time in find near=*',t2-t1,1)
!      call prinf('row_ptr=*',row_ptr,20)
!      call prinf('col_ind=*',col_ind,20)

!
! Estimate oversampling required for far-field, and oversample geometry
!
      ikerorder = 0

      allocate(novers(npatches),ixyzso(npatches+1))

      zpars = 0.0d0

      print *, "starting find far order"
      call cpu_time(t1)
!$       t1 = omp_get_wtime()      
!      call get_far_order(eps_quad,npatches,norders,ixyzs,iptype,cms, &
!         rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars, &
!         nnz,row_ptr,col_ind,rfac,novers,ixyzso)
      noveruse = 3
      npover = (noveruse+1)*(noveruse+2)/2
      do i=1,npatches
        novers(i) = noveruse
        ixyzso(i) = (i-1)*npover+1
      enddo
      ixyzso(npatches+1) = npatches*npover+1



      call cpu_time(t2)
!$       t2 = omp_get_wtime()      
      call prin2('far order estimation time=*',t2-t1,1)

      npts_over = ixyzso(npatches+1)-1

      print *, "oversampling factor=",(npts_over+0.0d0)/(npts+0.0d0)

!      call prinf('novers=*',novers,10)


      allocate(srcover(12,npts_over),wover(npts_over))

      call oversample_geom(npatches,norders,ixyzs,iptype,npts, &
        srccoefs,srcvals,novers,ixyzso,npts_over,srcover)

      call get_qwts(npatches,novers,ixyzso,iptype,npts_over, &
             srcover,wover)
!
!  Compute near quadrature corrections
!

      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad,3))
      wnear = 0
      iquadtype = 1
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call getnearquad_magnetostatics(npatches,norders,ixyzs, &
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id, &
        uvs_targ,eps_quad,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0, &
        nquad,wnear)
      call cpu_time(t2)
!$      t2 = omp_get_wtime()      
      call prin2('time taken in quadrature generation=*',t2-t1,1)
!
!
!  compute layer potential
!
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call lpcomp_grads0_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,rho, &
        novers,npts_over,ixyzso,srcover,wover,gradrho)
      call cpu_time(t2)
!$       t2 = omp_get_wtime()        
      call prin2('time taken in layer pot eval=*',t2-t1,1)  

      return
      end subroutine lpcomp_grads0
!
!
!
!
!
!
      subroutine lpcomp_grads0_addsub(npatches,norders,ixyzs,&
        iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
        eps,nnz,row_ptr,col_ind,iquad,nquad,wnear,rho, &
        novers,nptso,ixyzso,srcover,whtsover,gradrho)

!
!
!  This subroutine evaluates the potential 
!  1. grad S_{0}[rho]
!
!  where the near field is precomputed and stored
!  in the row sparse compressed format.
!
!  The identity term is not included in teh representation
!
!  The fmm is used to accelerate the far-field and 
!  near-field interactions are handled via precomputed quadrature
!
!  Using add and subtract - no need to call tree and set fmm parameters
!  can directly call existing fmm library
!
!  Input arguments:
! 
!    - npatches: integer
!        number of patches
!    - norders: integer(npatches)
!        order of discretization on each patch 
!    - ixyzs: integer(npatches+1)
!        ixyzs(i) denotes the starting location in srccoefs,
!        and srcvals array corresponding to patch i
!    - iptype: integer(npatches)
!        type of patch
!        iptype = 1, triangular patch discretized using RV nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (9,npts)
!        koornwinder expansion coefficients of xyz, dxyz/du,
!        and dxyz/dv on each patch. 
!        For each point 
!          * srccoefs(1:3,i) is xyz info
!          * srccoefs(4:6,i) is dxyz/du info
!          * srccoefs(7:9,i) is dxyz/dv info
!    - srcvals: real *8 (12,npts)
!        xyz(u,v) and derivative info sampled at the 
!        discretization nodes on the surface
!          * srcvals(1:3,i) - xyz info
!          * srcvals(4:6,i) - dxyz/du info
!          * srcvals(7:9,i) - dxyz/dv info
!          * srcvals(10:12,i) - normals info
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: double precision(ndtarg,ntarg)
!        target info, the first three components
!        must be xyz coordinates
!    - eps: real *8
!        precision requested
!    - nnz: integer *8
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(npts+1)
!        row_ptr(i) is the pointer
!        to col_ind array where list of relevant source patches
!        for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear_ij array where quadrature for col_ind(i)
!        starts for a single kernel. In this case the different kernels
!        are matrix entries are located at (m-1)*nquad+iquad(i), where
!        m is the kernel number
!    - nquad: integer
!        number of near field entries corresponding to each source target
!        pair. The size of wnear is 4*nquad since there are 4 kernels
!        per source target pair
!    - wnear: real *8(3*nquad)
!        Precomputed near field quadrature
!          * the first kernel is \nabla_{x} S_{0}
!          * the second kernel is \nabla_{y} S_{0}
!          * the third kernel is \nabla_{z} S_{0}
!    - rho: real *8 (npts)
!        The charge density rho
!    - novers: integer(npatches)
!        order of discretization for oversampled sources and
!        density
!    - ixyzso: integer(npatches+1)
!        ixyzso(i) denotes the starting location in srcover,
!        corresponding to patch i
!    - nptso: integer
!        total number of oversampled points
!    - srcover: real *8 (12,nptso)
!        oversampled set of source information
!    - whtsover: real *8 (nptso)
!        smooth quadrature weights at oversampled nodes
!
!  Output arguments:
!    - gradrho: real *8 (3,ntarg)
!         returns the potential grad S_{0}[\rho]
!

      implicit none
      integer npatches,norder,npols,npts
      integer ndtarg,ntarg
      integer norders(npatches),ixyzs(npatches+1)
      integer ixyzso(npatches+1),iptype(npatches)
      real *8 srccoefs(9,npts),srcvals(12,npts),eps
      real *8 targs(ndtarg,ntarg) 
      integer nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer iquad(nnz+1)
      real *8 rho(npts)
      real *8 wnear(nquad,3)

      integer novers(npatches)
      integer nover,npolso,nptso
      real *8 srcover(12,nptso),whtsover(nptso)
      real *8 gradrho(3,ntarg)
      real *8, allocatable :: wts(:)

      real *8 rhom,rhop,rmum,uf,vf,wtmp
      real *8 u1,u2,u3,u4,w1,w2,w3,w4,w5

      real *8, allocatable :: sources(:,:),targtmp(:,:)
      real *8, allocatable :: charges0(:)
      real *8, allocatable :: sigmaover(:),abc0(:)
      real *8, allocatable :: pot_aux(:),grad_aux(:,:)
      real *8 dpottmp,dgradtmp(3)
      real *8 vtmp1(3),vtmp2(3),vtmp3(3),rncj,errncj

      integer ns,nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val,E(4)

      integer i,j,jpatch,jquadstart,jstart,count1,count2
      real *8 radexp,epsfmm

      integer ipars(2)
      real *8 dpars(1),timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      real *8, allocatable :: ctmp0(:)
      real *8 thresh,ra,erra
      real *8 rr,rmin
      real *8 over4pi
      real *8 rbl(3),rbm(3)
      integer nss,ii,l,npover
      complex *16 ima,ztmp

      integer nd,ntarg0,nmax
      integer ndd,ndz,ndi,ier

      real *8 ttot,done,pi
      data ima/(0.0d0,1.0d0)/
      data over4pi/0.07957747154594767d0/

      parameter (ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

      ifpgh = 0
      ifpghtarg = 2
      allocate(sources(3,ns),targtmp(3,ntarg))

      nmax = 0

!
!  estimate max number of sources in the near field of any target
!
!    
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)

!
!  Allocate various densities
!

      allocate(sigmaover(ns),abc0(npts))

!
!  extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)
        sources(3,i) = srcover(3,i)
      enddo
!$OMP END PARALLEL DO      



!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targtmp(1,i) = targs(1,i)
        targtmp(2,i) = targs(2,i)
        targtmp(3,i) = targs(3,i)
      enddo
!$OMP END PARALLEL DO


      nd = 1
      allocate(charges0(ns))

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npts
        abc0(i) = rho(i)
      enddo
!$OMP END PARALLEL DO

      call oversample_fun_surf(nd,npatches,norders,ixyzs,iptype,& 
          npts,abc0,novers,ixyzso,ns,sigmaover)

        
!
!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,ns
        charges0(i) = sigmaover(i)*whtsover(i)*over4pi
      enddo
!$OMP END PARALLEL DO      

      allocate(pot_aux(ntarg),grad_aux(3,ntarg))

!      print *, "before fmm"
      ier =0 
      call lfmm3d_t_c_g(eps,ns,sources,charges0,ntarg,targtmp, &
        pot_aux,grad_aux,ier)


!$OMP PARALLEL DO DEFAULT(SHARED)         
      do i=1,ntarg
        gradrho(1:3,i) = grad_aux(1:3,i)
      enddo
!$OMP END PARALLEL DO      

!      print *, "after fmm"

!
!  Add near quadrature correction
!

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l,w1,w2,w3)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixyzs(jpatch+1)-ixyzs(jpatch)
          jquadstart = iquad(j)
          jstart = ixyzs(jpatch) 
          do l=1,npols
            w1 = wnear(jquadstart+l-1,1)
            w2 = wnear(jquadstart+l-1,2)
            w3 = wnear(jquadstart+l-1,3)
            gradrho(1,i) = gradrho(1,i) + w1*abc0(jstart+l-1)
            gradrho(2,i) = gradrho(2,i) + w2*abc0(jstart+l-1)
            gradrho(3,i) = gradrho(3,i) + w3*abc0(jstart+l-1)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO     
      
!      print *, "nmax=",nmax
!      print *, "nd=",nd

      call get_fmm_thresh(12,ns,srcover,3,ntarg,targtmp,thresh)

!      print *, "Thresh=",thresh
!
! Subtract near contributions computed via fmm
!
      allocate(ctmp0(nmax),srctmp2(3,nmax))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,l,jstart,nss,dpottmp,dgradtmp)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyzso(jpatch),ixyzso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)
            srctmp2(3,nss) = srcover(3,l)

            ctmp0(nss)=charges0(l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,targtmp(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        gradrho(1:3,i) = gradrho(1:3,i) - dgradtmp(1:3)
      enddo
!$OMP END PARALLEL DO      

!      print *, "finished pot eval"

      return
      end subroutine lpcomp_grads0_addsub
!
!
!
!

      subroutine lpeval_lap_g(nt,targs,npatches,centers,nv,verts, &
         triind,area,rnormals,charges,rfac,nover,eps,grad)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nt,nv,npatches,nover
      real *8, intent(in) :: targs(3,nt),centers(3,npatches)
      real *8, intent(in) :: verts(3,nv),area(npatches)
      real *8, intent(in) :: rnormals(3,npatches),rfac,eps
      real *8, intent(in) :: charges(npatches)
      real *8, intent(out) :: grad(3,nt)
      integer, intent(in) :: triind(3,npatches)

      integer, allocatable :: row_ptr(:),col_ind(:)
      real *8, allocatable :: wnear(:,:),srcover(:,:),rad_near(:)
      real *8, allocatable :: uvs(:,:),umat(:,:),vmat(:,:),wts(:)
      real *8, allocatable :: charges_over(:),pot_aux(:)
      integer, allocatable :: ixyzso(:)
      real *8, allocatable :: ctmp0(:),srctmp2(:,:)
      real *8 p0(3),p1(3),p2(3),vtmp(3),dgradtmp(3)
      real *8 over4pi

      data over4pi/0.07957747154594767d0/

      allocate(rad_near(npatches))

!$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,npatches
        rad_near(i) = sqrt(area(i))*rfac
      enddo
!$OMP END PARALLEL DO      
      
      ndtarg = 3
      nnz = 0
      ntarg = nt




      call findnearmem(centers,npatches,rad_near,ndtarg,targs,ntarg,nnz)
      allocate(row_ptr(ntarg+1),col_ind(nnz))

      call findnear(centers,npatches,rad_near,ndtarg,targs,ntarg, &
        row_ptr,col_ind)
      

      npols = (nover+1)*(nover+2)/2
      allocate(uvs(2,npols),umat(npols,npols),vmat(npols,npols))
      allocate(wts(npols))

      call vioreanu_simplex_quad(nover,npols,uvs,umat,vmat,wts)

      npts_over = npols*npatches

      allocate(srcover(3,npts_over),charges_over(npts_over))
      do ipatch=1,npatches
        p0(1:3) = verts(1:3,triind(1,ipatch))
        p1(1:3) = verts(1:3,triind(2,ipatch))
        p2(1:3) = verts(1:3,triind(3,ipatch))
        do j=1,npols
          ipt = (ipatch-1)*npols + j
          u = uvs(1,j)
          v = uvs(2,j)
          srcover(1:3,ipt) = p0(1:3) + u*(p1(1:3) - p0(1:3)) + &
            v*(p2(1:3)-p0(1:3))
          charges_over(ipt) = charges(ipatch)*area(ipatch)*wts(j)*over4pi
        enddo
      enddo

      allocate(pot_aux(ntarg))
      pot_aux = 0
      grad = 0
      
      call cpu_time(t1)
!$      t1 = omp_get_wtime()      
      call lfmm3d_t_c_g(eps,npts_over,srcover,charges_over,ntarg,targs, &
        pot_aux,grad,ier)

      call cpu_time(t2)
!$      t2 = omp_get_wtime()      
      call prin2('fmm time=*',t2-t1,1)
      
      

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ipatch,p0,p1,p2) &
!$OMP PRIVATE(i,jtarg,vtmp)
      do jtarg=1,ntarg
        do i=row_ptr(jtarg),row_ptr(jtarg+1)-1
          ipatch = col_ind(i)
          p0(1:3) = verts(1:3,triind(1,ipatch))
          p1(1:3) = verts(1:3,triind(2,ipatch))
          p2(1:3) = verts(1:3,triind(3,ipatch))

          vtmp(1:3) = 0
          call pot_int2(p0,p1,p2,rnormals(1,ipatch),targs(1,jtarg),vtmp)
          vtmp(1:3) = vtmp(1:3)*over4pi*charges(ipatch)

          grad(1:3,jtarg) = grad(1:3,jtarg) + vtmp(1:3)
        enddo
      enddo
!$OMP END PARALLEL DO
      
      call get_fmm_thresh(3,npts_over,srcover,3,ntarg,targs,thresh)

!      print *, "Thresh=",thresh
!
! Subtract near contributions computed via fmm
!
      allocate(ixyzso(npatches+1))
!$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npatches+1
        ixyzso(i) = (i-1)*npols + 1
      enddo
!$OMP END PARALLEL DO
      
      call get_near_corr_max(ntarg,row_ptr,nnz,col_ind,npatches,ixyzso,&
        nmax)
      nd = 1
      ntarg0 = 1

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

            ctmp0(nss)=charges_over(l)
          enddo
        enddo

        dpottmp = 0
        dgradtmp = 0

        call l3ddirectcg(nd,srctmp2,ctmp0,nss,targs(1,i), &
          ntarg0,dpottmp,dgradtmp,thresh)
        grad(1:3,i) = grad(1:3,i) - dgradtmp(1:3)
      enddo
!$OMP END PARALLEL DO      

!      print *, "finished pot eval"

      return
      end subroutine lpeval_lap_g
!
!
!
!
!
!
      subroutine get_areas(npatches,nv,verts,triind,area)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: npatches,nv,triind(3,npatches)
      real *8, intent(in) :: verts(3,nv)
      real *8, intent(out) :: area(npatches)

      real *8 vtmp1(3),vtmp2(3),p0(3),p1(3),p2(3),wtmp(3)

      area = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,vtmp1,vtmp2,wtmp,p0,p1,p2)
      do i=1,npatches
        p0(1:3) = verts(1:3,triind(1,i))
        p1(1:3) = verts(1:3,triind(2,i))
        p2(1:3) = verts(1:3,triind(3,i))

        vtmp1(1:3) = p1(1:3)-p0(1:3)
        vtmp2(1:3) = p2(1:3)-p0(1:3)

        call cross_prod3d(vtmp1,vtmp2,wtmp)
        area(i) = sqrt(wtmp(1)**2 + wtmp(2)**2 + wtmp(3)**2)
      enddo
!$OMP END PARALLEL DO


      end subroutine get_areas





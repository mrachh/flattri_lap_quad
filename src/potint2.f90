
!
      subroutine pot_int2(r1,r2,r3,normal,targ,vpot)
      implicit real *8 (a-h,o-z)
      real *8 r1(3),r2(3),r3(3),normal(3),targ(3),vpot(3)
      real *8 ri(3),rs(9),betaterm(3),rr(3,6),rl(3),ru(9)
      real *8 rm(3,3),vtmp(3),rfl(3),rrl(3,3)
!
!  target depedent stuff
!
      real *8 rndot, rp(3),rmtarg(3,3),rcheck(3),rtu(3),rtp(3),rtmp(3)
      real *8 rpp(3),rpm(3),rphat(3)

      rr(1:3,1) = r2(1:3)
      rr(1:3,2) = r1(1:3)
      rr(1:3,3) = r3(1:3)
      rr(1:3,4) = r1(1:3)
      rr(1:3,5) = r3(1:3)
      rr(1:3,6) = r2(1:3)
      
      rl(1) = (r2(1) - r1(1))**2 + (r2(2)-r1(2))**2 + (r2(3)-r1(3))**2
      rl(1) = sqrt(rl(1))
      
      rl(2) = (r3(1) - r1(1))**2 + (r3(2)-r1(2))**2 + (r3(3)-r1(3))**2
      rl(2) = sqrt(rl(2))

      rl(3) = (r2(1) - r3(1))**2 + (r2(2)-r3(2))**2 + (r2(3)-r3(3))**2
      rl(3) = sqrt(rl(3))

      rrl(1:3,1) =( r2(1:3)-r1(1:3))/rl(1)
      rrl(1:3,2) =( r3(1:3)-r1(1:3))/rl(2)
      rrl(1:3,3) =( r3(1:3)-r2(1:3))/rl(3)

      call cross_prod3d(rrl(1,1),normal,ru(1))
      call cross_prod3d(rrl(1,2),-normal,ru(4))
      call cross_prod3d(rrl(1,3),normal,ru(7))

      rndot = targ(1)*normal(1) + targ(2)*normal(2) + targ(3)*normal(3)
      rp(1:3) = targ(1:3) - rndot*normal(1:3)

      rm(1:3,1) = 0.5d0*(r1(1:3) + r2(1:3))
      rm(1:3,2) = 0.5d0*(r1(1:3) + r3(1:3))
      rm(1:3,3) = 0.5d0*(r2(1:3) + r3(1:3))
      
      rmtarg(1:3,1) = rm(1:3,1) - rp(1:3)
      rmtarg(1:3,2) = rm(1:3,2) - rp(1:3)
      rmtarg(1:3,3) = rm(1:3,3) - rp(1:3)

      rf = 1.0d-6
      rfl(1:3) = rl(1:3)*rf

      rcheck(1) = rmtarg(1,1)*ru(1) + rmtarg(2,1)*ru(2) + rmtarg(3,1)*ru(3)
      rcheck(2) = rmtarg(1,2)*ru(4) + rmtarg(2,2)*ru(5) + rmtarg(3,2)*ru(6)
      rcheck(3) = rmtarg(1,3)*ru(7) + rmtarg(2,3)*ru(8) + rmtarg(3,3)*ru(9)
      
      rtu(1:3) = targ(1:3)
      if(abs(rcheck(1)).le.rfl(1)) then
        rtu(1:3) = rtu(1:3) - rfl(1)*ru(1:3)
      endif
      
      if(abs(rcheck(2)).le.rfl(2)) then
        rtu(1:3) = rtu(1:3) - rfl(2)*ru(4:6)
      endif
      
      if(abs(rcheck(3)).le.rfl(3)) then
        rtu(1:3) = rtu(1:3) - rfl(3)*ru(7:9)
      endif
      
      rndot = rtu(1)*normal(1) + rtu(2)*normal(2) + rtu(3)*normal(3)
      rtp(1:3) = rtu(1:3) - rndot*normal(1:3)


      icount = 0
      do i=1,3
        rtmp(1:3) = rtu(1:3) - rr(1:3,icount+2)
        rd = rtmp(1)*normal(1) + rtmp(2)*normal(2) + rtmp(3)*normal(3)
        d1 = normal(1)*rr(1,icount+1) + normal(2)*rr(2,icount+1) + & 
          normal(3)*rr(3,icount+1)
        d2 = normal(1)*rr(1,icount+2) + normal(2)*rr(2,icount+2) + &
         normal(3)*rr(3,icount+2)
        rpp(1:3) = rr(1:3,icount+1) - d1*normal(1:3)      
        rpm(1:3) = rr(1:3,icount+2) - d2*normal(1:3)       
        rlp = rrl(1,i)*(rpp(1)-rtp(1)) + rrl(2,i)*(rpp(2)-rtp(2)) + &
          rrl(3,i)*(rpp(3)-rtp(3))
        rlm = rrl(1,i)*(rpm(1)-rtp(1)) + rrl(2,i)*(rpm(2)-rtp(2)) + &
          rrl(3,i)*(rpm(3)-rtp(3))
        rp0 = abs(ru(1+3*(i-1))*(rpm(1)-rtp(1)) + ru(2+3*(i-1))*(rpm(2)-rtp(2)) + &
          ru(3+3*(i-1))*(rpm(3)-rtp(3)))

        pplus = sqrt(rp0*rp0 + rlp*rlp)
        pminus = sqrt(rp0*rp0 + rlm*rlm)


        rphat(1:3) = (rpm(1:3) - rtp(1:3) - rrl(1:3,i)*rlm)/rp0
        rplus = sqrt(pplus*pplus +rd*rd)
        rminus = sqrt(pminus*pminus + rd*rd)
        r0 = sqrt(rp0*rp0 + rd*rd)
        
        rdabs = abs(rd)
        d1 = atan(rp0*rlp/(r0*r0 + rdabs*rplus))
        d2 = atan(rp0*rlm/(r0*r0 + rdabs*rminus))
        d3 = log((rplus + rlp)/(rminus + rlm))

        dotphatu = rphat(1)*ru(1+3*(i-1)) + rphat(2)*ru(2+3*(i-1)) + &
          rphat(3)*ru(3+3*(i-1))

        betaterm(i) = dotphatu*(d1-d2)
        rs(1+3*(i-1)) = d3*ru(1+3*(i-1))
        rs(2+3*(i-1)) = d3*ru(2+3*(i-1))
        rs(3+3*(i-1)) = d3*ru(3+3*(i-1))
        icount = icount+2
      enddo

      beta = betaterm(1) + betaterm(2) + betaterm(3)
      ri(1) = rs(1) + rs(4) + rs(7)
      ri(2) = rs(2) + rs(5) + rs(8)
      ri(3) = rs(3) + rs(6) + rs(9)

      rsign = 1.0d0
      if(rd.le.0) rsign = -1.0d0
      
      vpot(1:3) = -normal(1:3)*rsign*beta - ri(1:3)

      
      return
      end subroutine pot_int2

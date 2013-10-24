      subroutine initcomms()

c*********************************************************************
c
c     dummy initcomms routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      return
      end

      subroutine machine(idnode,mxnode)

c*********************************************************************
c
c     dummy machine routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      implicit none

      integer idnode,mxnode

      idnode=0
      mxnode=1

      return
      end

      integer function mynode()

c*********************************************************************
c
c     dummy mynode routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      implicit none

      mynode=0

      return
      end

      integer function nodedim()

c*********************************************************************
c
c     dummy nodedim routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      implicit none

      nodedim=0

      return
      end

      integer function numnodes()

c*********************************************************************
c
c     dummy numnodes routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      implicit none

      numnodes=1

      return
      end

      subroutine csend(msgtag,buf,length,pe,idum)

c*********************************************************************
c
c     dummy csend routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c********************************************************************

      implicit none

      integer msgtag,length,pe,idum

      real(8) buf(*)

      return
      end
      
      subroutine crecv(msgtag,buf,length)

c*********************************************************************
c
c     dummy crecv routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer msgtag,length
      real(8) buf(*)

      return
      end

      subroutine gisum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dummy isum for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c***********************************************************************
      
      implicit none

      integer nnn
      integer aaa(*),bbb(*)

      return
      end

      subroutine gdsum(aaa,nnn,bbb)

c***********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c***********************************************************************

      implicit none

      integer nnn
      real(8) aaa(*),bbb(*)

      return
      end

      subroutine gimax(aaa,nnn,bbb)

c***********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c***********************************************************************
      
      implicit none

      integer nnn
      integer aaa(*),bbb(*)

      return
      end

      subroutine gstate(check)

c***********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c***********************************************************************

      implicit none

      logical check

      return
      end

      subroutine gsync()

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      return
      end

      subroutine exitcomms()

c*********************************************************************
c
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      stop
      end

      subroutine merge(idnode,mxnode,natms,nbuff,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************
      implicit none

      integer idnode,mxnode,natms,nbuff
      real(8) xxx(*),yyy(*),zzz(*),buffer(*)

      return
      end

      subroutine merge1(idnode,mxnode,natms,lstme,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,natms
      integer lstme(*)
      real(8) xxx(*),yyy(*),zzz(*),buffer(*)

      return
      end

      subroutine merge4(idnode,mxnode,ngrp,nbuff,q0,q1,q2,q3,buffer)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ngrp,nbuff
      real(8) q0(*),q1(*),q2(*),q3(*),buffer(*)

      return
      end

      subroutine shlmerge(idnode,mxnode,ntshl)
      
c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer idnode,mxnode,ntshl

      return
      end

      subroutine shmove
     x     (idnode,mxnode,natms,lashap,lishap,xxt,yyt,zzt,
     x      txx,tyy,tzz,buffer)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer idnode, mxnode, natms
      integer lishap(*),lashap(*)
      real(8) xxt(*),yyt(*),zzt(*),txx(*),tyy(*),tzz(*),buffer(*)

      return
      end

      subroutine splice
     x      (idnode,natms,listme,listot,xxx,yyy,zzz,buffer)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      implicit none

      integer idnode,natms
      integer listme(*),listot(*)
      real(8) xxx(*),yyy(*),zzz(*),buffer(*)

      return
      end

      subroutine passcon
     x     (lshmov,idnode,mxnode,natms,nscons,lashap,lishap,listme,
     x     listin,listot,listcon,lstfrz)
      
c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      use setup_module

      implicit none

      logical lshmov
      integer idnode,mxnode,natms,nscons,i,j,k
      integer lashap(mxproc),lishap(mxlshp),listme(mxatms)
      integer listin(mxatms),listot(mxatms),listcon(mxcons,3)
      integer lstfrz(mxatms)
      
      do i=1,natms
         listme(i)=0
      enddo
      
      do k=1,nscons
         
         i=listcon(k,2)
         j=listcon(k,3)
         listme(i)=listme(i)+1
         listme(j)=listme(j)+1
         
      enddo

c     keep record of all atoms subject to constraints
      
      do i=1,natms
         
         if(listme(i).gt.0)then
            listot(i)=1
         else
            listot(i)=0
         endif
         
      enddo
      
      return
      end

      subroutine passpmf
     x  (idnode,mxnode,natms,nspmf,listpm,listin,lstpmt,lstpmf,npmf)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      use error_module
      use setup_module

      implicit none

      integer idnode,mxnode,natms,nspmf,i,j,k
      integer listpm(mxpmf),listin(mxatms),lstpmt(mxpmf)
      integer lstpmf(mxspmf,mspmf),npmf(2)

      if(mxpmf.lt.natms) call error(idnode,490)

      do i=1,natms
        listpm(i)=0
      enddo
      
      do k=1,nspmf
        
        do j = 1,npmf(1)+npmf(2)

          i=lstpmf(j,k)
          listpm(i)= 1
          
        enddo

      enddo

c     keep record of all atoms subject to pmf constraints
      
      do i=1,natms
        
        if(listpm(i).gt.0)then
          lstpmt(i)=1
        else
          lstpmt(i)=0
        endif
        
      enddo
      
      return
      end

      subroutine passquat
     x  (lcnb,idnode,mxnode,natms,ngrp,nscons,ntpmls,listin,
     x  listcon,lstrgd,lstout,lstcsit,lstgtp,nummols,numgrp,numgsit)

c*********************************************************************
c     
c     dummy routine for serial DL_POLY
c     copyright daresbury laboratory
c     author - w.smith
c
c*********************************************************************

      use setup_module

      implicit none

      logical lcnb
      integer idnode,mxnode,natms,ngrp,nscons,ntpmls,jj,id,igrp1,igrp2
      integer i,j,k,jr,igrp,itmols,imols,ik,lgrp
      integer listin(mxatms),listcon(mxcons,3),lstrgd(mxgatm)
      integer lstout(mxatms),lstcsit(2*mxcons),numgsit(mxungp)
      integer lstgtp(mxgrp),nummols(mxtmls),numgrp(mxtmls)
      
c     block indices for groups
      
      igrp1 = (idnode*ngrp)/mxnode + 1
      igrp2 = ((idnode+1)*ngrp)/mxnode
      
c     locate site indices of atoms in constraints

      do i = 1,natms
        listin(i) = 0
      enddo

c     loop over molecule types

      jr = 0 
      igrp = 0
      do itmols=1,ntpmls

c     loop over molecules in system
        
        do imols=1,nummols(itmols)

c     construct rigid body site list: each processor has a different copy
          
          do lgrp=1,numgrp(itmols)
            
            igrp=igrp+1
            
            if((igrp.ge.igrp1).and.(igrp.le.igrp2)) then
                
              id = lstgtp(igrp)
              do jj = 1,numgsit(id)
                  
                jr = jr +1
                i = lstrgd(jr)
                listin(i) = jj

              enddo
            endif
          enddo
        enddo
      enddo

      lcnb = .true.
      ik = 0
      do k=1,nscons
        
        i=listcon(k,2)
        j=listcon(k,3)

        if(listin(i).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(i)
          lcnb = .false.
        endif

        if(listin(j).ne.0) then
          ik = ik + 1
          lstcsit(ik) = listin(j)
          lcnb = .false.
        endif

      enddo

c     lcnb flags bodies connected by constraints

      lcnb = (.not.lcnb)
      
      return
      end


      module setup_module
      
c***********************************************************************
c     
c     dl_poly module for defining default array sizes
c     
c copyright (c) 2033, daresbury laboratory
c     author    - w. smith    sep 2003
c all rights reserved.
c 
c redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c 
c 1. redistributions of source code must retain the above copyright notice,
c    this list of conditions and the following disclaimer.
c 2. redistributions in binary form must reproduce the above copyright
c    notice, this list of conditions and the following disclaimer in the
c    documentation and/or other materials provided with the distribution.
c 3. neither the name of the <organization> nor the names of its
c    contributors may be used to endorse or promote products derived from
c    this software without specific prior written permission.
c 
c this software is provided by the copyright holders and contributors "as is"
c and any express or implied warranties, including, but not limited to, the
c implied warranties of merchantability and fitness for a particular purpose
c are disclaimed. in no event shall the copyright owner or contributors be
c liable for any direct, indirect, incidental, special, exemplary, or
c consequential damages (including, but not limited to, procurement of
c substitute goods or services; loss of use, data, or profits; or business
c interruption) however caused and on any theory of liability, whether in
c contract, strict liability, or tort (including negligence or otherwise)
c arising in any way out of the use of this software, even if advised of the
c possibility of such damage.
c     note the following internal units apply everywhere
c     
c     unit of time      (to)    =          1 x 10**(-12) seconds
c     unit of length    (lo)    =          1 x 10**(-10) metres
c     unit of mass      (mo)    = 1.6605402  x 10**(-27) kilograms
c     unit of charge    (qo)    = 1.60217733 x 10**(-19) coulombs
c     unit of energy    (eo)    = 1.6605402  x 10**(-23) joules
c     unit of pressure  (po)    = 1.6605402  x 10**(  7) pascals
c     
c*********************************************************************
      
      use parse_module
      
      implicit none
      
c     FIXED PARAMETERS
      
c     standard pi values
      
      real(8), parameter :: pi=3.141592653589793d0
      real(8), parameter :: sqrpi=1.7724538509055159d0
      
c     conversion factor for coulombic terms in internal units
c     i.e. (unit(charge)**2/(4 pi eps0 unit(length))/unit(energy)
      
      real(8), parameter :: r4pie0=138935.4835d0
      
c     boltzmann constant in internal units
      
      real(8), parameter :: boltz=8.31451115d-1
      
c     planck's constant in internal units
      
      real(8), parameter :: hbar=6.350780719d0
      
c     conversion factor for pressure from internal units to katm
      
      real(8), parameter :: prsunt=0.163882576d0
      
c     main input channel
      
      integer, parameter :: nread=5
      
c     main output channel
      
      integer, parameter :: nrite=6
      
c     force field input channel
      
      integer, parameter :: nfield=9
      
c     configuration file input channel
      
      integer, parameter :: nconf=10
      
c     statistical data file output channel
      
      integer, parameter :: nstats=20
      
c     trajectory history file channel
      
      integer, parameter :: nhist=21
      
c     acummulators restart dump file
      
      integer, parameter :: nrest=22
      
c     tabulated potential file channel
      
      integer, parameter :: ntable=23
      
c     rdf file channel number
      
      integer, parameter :: nrdfdt=24
      
c     z density file channel number
      
      integer, parameter :: nzdndt=25
      
c     hyperdynamics reference basin file
      
      integer, parameter :: nbsn=30
      
c     neb reaction path profile file
      
      integer, parameter :: npro=31
      
c     hyperdynamics events journal file
      
      integer, parameter :: nevnt=33
      
c     hyperdynamics tracking file
      
      integer, parameter :: ntrk=32
      
c     hyperdynamics restart file
      
      integer, parameter :: nhrs=35
      
c     free energy data file
      
      integer, parameter :: nfrnwr=41
      
c     solvation data file
      
      integer, parameter :: nsolwr=43
      
c     data dumping interval in event of system crash
      
      integer, parameter :: ndump=1000
      
c     maximum number of neb calculations
      
      integer, parameter :: maxneb=10
      
c     array allocation parameters (set by subroutine parset)
      
      integer kmaxa,kmaxb,kmaxc,minnode,msatms,msbad,msgrp
      integer mspmf,msteth,mxangl,mxatms,mxbond,mxbuff,mxcell
      integer mxcons,mxdihd,mxewld,mxexcl,mxfbp,mxfld,mxgatm,mxgrid
      integer mxgrp,mxinv,mxlist,mxlshp,mxneut,mxngp,mxnstk,mxpang
      integer mxpbnd,mxpdih,mxpfbp,mxpinv,mxpmf,mxproc,mxptbp,mxpvdw
      integer mxrdf,mxzdn,mxshl,mxsite,mxspmf,mxstak,mxtang,mxtbnd
      integer mxtbp,mxtcon,mxtdih,mxteth,mxtinv,mxtmls,mxtshl,mxungp
      integer mxvdw,mxxdf,mx2tbp,mx3fbp,mxebuf,mxquat,mxshak,mxspl
      integer kmaxd,kmaxe,kmaxf,mxspme,mxftab,mxhko,mxegrd,mxhke
      integer mxmet,mxsmet,mxpmet,mxter,mxpter,mxatyp,mxxtyp
      integer mxtmls_fre,mxewld_fre,mxebuf_fre,mxatms_fre,mxatyp_exc
      integer mxtmls_exc,mxtmls_sol,mxebuf_sol,mxatms_sol
      
      save kmaxa,kmaxb,kmaxc,minnode,msatms,msbad,msgrp
      save mspmf,msteth,mxangl,mxatms,mxbond,mxbuff,mxcell
      save mxcons,mxdihd,mxewld,mxexcl,mxfbp,mxfld,mxgatm,mxgrid
      save mxgrp,mxinv,mxlist,mxlshp,mxneut,mxngp,mxnstk,mxpang
      save mxpbnd,mxpdih,mxpfbp,mxpinv,mxpmf,mxproc,mxptbp,mxpvdw
      save mxrdf,mxzdn,mxshl,mxsite,mxspmf,mxstak,mxtang,mxtbnd
      save mxtbp,mxtcon,mxtdih,mxteth,mxtinv,mxtmls,mxtshl,mxungp
      save mxvdw,mxxdf,mx2tbp,mx3fbp,mxebuf,mxquat,mxshak,mxspl
      save kmaxd,kmaxe,kmaxf,mxspme,mxftab,mxhko,mxegrd,mxhke
      save mxmet,mxsmet,mxpmet,mxter,mxpter,mxatyp,mxxtyp
      save mxtmls_fre,mxewld_fre,mxebuf_fre,mxatms_fre,mxatyp_exc
      save mxtmls_exc,mxtmls_sol,mxebuf_sol,mxatms_sol
      
      contains
      
      subroutine parset(redirect,idnode,mxnode,buffer)
      
c***********************************************************************
c     
c     dl_poly subroutine to determine required array sizes for 
c     allocation of memory manager
c     
c     copyright daresbury laboratory 1997
c     author - w.smith june 1997
c     
c***********************************************************************
      
      logical loglnk,lewald,lspme,lhke,nolink,lcshft
      logical lsolva,lfree,lfrmas,lghost,redirect
      real(8) cell,celprp,rctbp,rcfbp,volm,xhi,yhi,zhi,rcut,rvdw
      real(8) densvar,delr,cut,dens,ratio,drdf,dzdn,rcter,buffer
      real(8) zlen
      integer imcon,nhko,ilx,ily,ilz,ncells
      integer idnode,mxnode,mxn1
      
      dimension cell(9),celprp(10),buffer(10)
      
      lhke=.false.
      lspme=.false.
      lewald=.false.
      lcshft=.false.
      nolink=.false.
      redirect=.false.
      mxtmls_sol=1
      mxebuf_sol=1
      mxatms_sol=1
      mxtmls_fre=1
      mxewld_fre=1
      mxebuf_fre=1
      mxatms_fre=1
      mxatyp_exc=1
      mxtmls_exc=1

c     specify maximum and minimum nodes
      
      mxproc=mxnode
      minnode=mxnode
      
c     scan the FIELD file data
      
      call fldscan(idnode,mxn1,rctbp,rcfbp,rcter)
      
c     scan CONFIG file data
      
      call cfgscan
     x  (idnode,nconf,imcon,volm,xhi,yhi,zhi,cell,buffer)
      
c     scan CONTROL file data
      
      call conscan
     x  (redirect,lewald,lspme,lhke,nolink,lcshft,lsolva,lfree,lfrmas,
     x  lghost,idnode,imcon,nhko,rcut,rvdw,delr,densvar,drdf,dzdn,
     x  zlen,cell)
      
c     set dimension of working coordinate arrays
      
      msatms=max(1,(mxatms+minnode-1)/minnode)
      if(lsolva)mxatms_sol=mxatms
      if(lfree.or.lghost)mxatms_fre=mxatms
      
c     maximum number of molecule types
      
      mxtmls=max(mxtmls,1)
      if(lsolva)mxtmls_sol=mxtmls
      if(lfree)mxtmls_fre=mxtmls
      if(lghost)then
        
        mxtmls_exc=mxtmls
        mxtmls_fre=mxtmls
        
      endif
      
c     maximum number of specified bondlength constraints
      
      mxtcon=max(mxtcon,1)
      
c     maximum number of chemical bond potentials
      
      mxtbnd=max(mxtbnd,1)
      
c     maximum number of different bond angle potentials
      
      mxtang=max(mxtang,1)
      
c     maximum number of different torsional potentials
      
      mxtdih=max(mxtdih,1)
      
c     maximum number of different inversion potentials
      
      mxtinv=max(mxtinv,1)
      
c     maximum number of unique rigid body units
      
      mxungp=max(mxungp,1)
      
c     maximum number of tethered atom potentials
      
      mxteth=max(mxteth,1)
      
c     maximum number of core-shell units
      
      mxshl=max(mxshl,1)
      
c     set maximum number of unique atom types
      
      mxatyp=max(1,mxatyp)
      mxxtyp=(mxatyp*(mxatyp+1))/2
      if(lghost)mxatyp_exc=mxatyp
      
c     maximum number of vdw potentials
      
      mxvdw=max(mxvdw,1)+1
      
c     maximum number of metal potentials
      
      mxmet=max(mxmet,1)+1
      mxsmet=mxatyp
      
c     maximum number of tersoff potentials
      
      if(mxter.gt.0)then
        
        mxter=mxatyp
        
      endif
      
c     maximum number of three body potentials
      
      if(mxtbp.eq.0)then
        
        mx2tbp=0
        
      else
        
        mx2tbp=(mxatyp*(mxatyp+1))/2
        mxtbp=mx2tbp*mxatyp
        
      endif
      
c     maximum number of four body potentials
      
      if(mxfbp.eq.0)then
        
        mx3fbp=0
        
      else
        
        mx3fbp=(mxatyp*(mxatyp+1)*(mxatyp+2))/6
        mxfbp=mxatyp*mx3fbp
        
      endif
      
c     maximum number of angular potential parameters

      mxpang=6

c     maximum number of three body potential parameters
      
      mxptbp=mxpang+1
      
c     maximum number of four body potential parameters
      
      mxpfbp=3
      
c     maximum number of parameters for dihedrals
      
      mxpdih=5
      
c     maximum number of parameters for inversion potentials
      
      mxpinv=2
      
c     maximum number of parameters for bond potentials
      
      mxpbnd=4
      
c     maximum number of parameters for vdw potentials
      
      mxpvdw=6
      
c     maximum number of parameters for metal potentials
      
      mxpmet=7
      
c     maximum number of parameters for tersoff potential
      
      mxpter=11
      
c     maximum number of external field parameters
      
      mxfld=10
      
c     maximum number of excluded atoms per atom
      
      mxexcl=max(mxexcl,1)
      
c     maximum number of different sites in system
      
      mxsite=max(mxsite,1)
      
c     maximum number of chemical bonds per node
      
      mxbond=max(1,(mxbond+minnode-1)/minnode)
      
c     maximum number of bond angles per node
      
      mxangl=max(1,(mxangl+minnode-1)/minnode)
      
c     maximum number of torsion angles per node
      
      mxdihd=max(1,(mxdihd+minnode-1)/minnode)
      
c     maximum number of inversion potentials per node
      
      mxinv=max(1,(mxinv+minnode-1)/minnode)
      
c     maximum number of constraints per node
      
      mxcons=max(1,2*((mxcons+minnode-1)/minnode))
      
c     maximum number of tethered atoms per node
      
      msteth=max(1,(msteth+minnode-1)/minnode)
      
c     maximum size for working arrays for bonds, angles, dihedrals
c     inversion potentials, tethers and core-shell units
      
      msbad=max(mxbond,mxangl,mxdihd,mxinv,msteth,mxshl)
      
c     maximum number of grid points in potentials arrays
      
      if(mxgrid.eq.0)then
        
        mxgrid=max(1000,int(rvdw/0.01d0+0.5d0)+4)
        
      endif
      
      mxegrd=0
      if(lewald.or.lspme.or.lhke.or.lcshft)mxegrd=mxgrid
      
c     maximum dimension of rdf arrays
      
      mxrdf=max(128,int(rcut/drdf))
      
c     maximum dimension of zdensity arrays
      
      mxzdn=max(128,int(zlen/dzdn))
      
c     maximum number of rigid groups in system
      
      mxgrp=max(mxgrp,1)
      
c     maximum number of rigid groups per node
      
      msgrp=max(1,(mxgrp+minnode-1)/minnode)
      
c     maximum number of sites per rigid unit
      
      mxngp=max(mxngp,3)
      
c     maximum number of sites in rigid units
      
      mxgatm=max(1,mxgatm)
      
c     maximum number of timesteps in stack arrays
      
      mxstak=max(100,mxstak)
      
c     maximum number of variables in stack arrays
      
      mxnstk=45+mxatyp
      
c     dimension of shake shared atoms array
      
      mxlshp=max(mxcons*2,1)
      
c     set dimension of working arrays in ewald sum
      
      mxewld=1
      mxebuf=1
      if(lewald)then
        
        mxftab=1
        mxewld=msatms
        mxebuf=(2*kmaxa+1)*(2*kmaxb+1)*(2*kmaxc+1)-1
        if(lfree.or.lghost)mxebuf=3*mxebuf
        if(mxnode.le.16.and.mxebuf.le.5000)mxebuf=1
        
      endif
      
c     set dimension of working arrays in spme
      
      mxspme=1
      if(lspme)then
        
        mxspme=mxatms
        mxftab=2*(kmaxd+kmaxe+kmaxf)
        
      endif
      
c     set dimension of working arrays for HK ewald
      
      mxhko=1
      mxhke=1
      if(lhke)then
        
        mxhko=2
        mxewld=msatms
        mxhke=msatms
        if(nhko.gt.0)mxhko=max(2,nhko)
        mxebuf=(2*kmaxa+1)*(2*kmaxb+1)-1
        if(mxnode.le.16.and.mxebuf.le.5000)mxebuf=1
        
      endif
      
      if(lsolva)mxebuf_sol=mxebuf
      if(lfree.or.lghost)then
        
        mxebuf_fre=mxebuf
        mxewld_fre=mxewld
        
      endif

c     maximum dimension of principal transfer buffer
      
      mxbuff=max(6*mxatms,8*(mxcons+1),8*(mxgrp+1),mxnstk*mxstak,
     x  mxebuf,mxgrid,2*kmaxa*kmaxb*kmaxc,2*kmaxd*kmaxe*kmaxf,
     x  10000)
      
c     maximum size of verlet neighbour/link cell list for each atom
c     decide if link-cells in use or not
      
      cut=rcut+delr
      dens=dble(mxatms)/volm
      ratio=1.5d0*dens*(4.d0*pi/3.d0)*cut**3
      mxlist=min(nint(ratio),(mxatms+1)/2)
      if(imcon.eq.0) then
        
        cell(1)=max(xhi+2.d0*cut,3.d0*cut)
        cell(5)=max(yhi+2.d0*cut,3.d0*cut)
        cell(9)=max(zhi+2.d0*cut,3.d0*cut)
        
      endif
      if(imcon.eq.6)then
        
        cell(9)=max(zhi+2.d0*cut,3.d0*cut,cell(9))
        
      endif
      
      if(nolink)then
        
        loglnk=.false.
        
      else
        
        loglnk=.true.
        call dcell(cell,celprp)
        ilx=int(celprp(7)/cut)
        ily=int(celprp(8)/cut)
        ilz=int(celprp(9)/cut)
        if(ilx.lt.3.or.ily.lt.3.or.ilz.lt.3)loglnk=.false.
        ncells=ilx*ily*ilz
        if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7)loglnk=.false.
        if(mxneut.gt.0.and.ncells.le.36) loglnk=.false.
        
      endif
      
      mxcell=1
      if(loglnk)then
        
        mxlist=14*nint(1.5d0*dens*celprp(10)/dble(ncells))
        mxcell=(ilx+2)*(ily+2)*(ilz+2)
        
      endif
      
      if(mxneut.gt.0)mxlist=(mxneut+1)/2
      mxlist=2*mxlist
      if(mxtbp.gt.0.or.mxfbp.gt.0.or.mxter.gt.0)then
        
        if(mxtbp.gt.0)cut=min(cut,rctbp)
        if(mxfbp.gt.0)cut=min(cut,rcfbp)
        if(mxter.gt.0)cut=min(cut,rcter)
        ilx=max(3,int(celprp(7)/cut))
        ily=max(3,int(celprp(8)/cut))
        ilz=max(3,int(celprp(9)/cut))
        mxcell=max(mxcell,(ilx+2)*(ily+2)*(ilz+2))
        
      endif
      mxcell=int(dble(mxcell)*densvar/100.d0)
      mxlist=int(dble(mxlist)*densvar/100.d0)
      mxlist=max(500,mxlist)
      
c     maximum size for coordinate difference arrays
      
      mxxdf=max(mxlist,mxatms,mxcons,mxn1*mxn1*(mxneut+1)/2)
      
c     maximum number of core-shell unit types
      
      mxtshl=max(mxtshl,1)
      
c     potential of mean force array parameter
      
      mxpmf=max(mxpmf,1)
      
c     number of pmf constraints on a processor
      
      mspmf=max(1,(mxpmf+minnode-1)/minnode)
      
c     maximum number of sites to define pmf units
      
      mxspmf=max(mxspmf,1)
      
c     maximum iterations in quaternion integration
      
      mxquat=100
      
c     maximum number of shake cycles
      
      mxshak=100
      
c     maximum b-spline interpolation order
      
      mxspl=12
      
c     increment mxneut
      
      if(mxneut.gt.0)mxneut=mxneut+1
      
      return
      
      end subroutine parset
      
      subroutine fldscan(idnode,mxn1,rctbp,rcfbp,rcter)
      
c***********************************************************************
c     
c     dl_poly routine for scanning the field file to determine the
c     required parameters
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c     
c***********************************************************************
      
      integer, parameter :: mmk=1000
      
      character*8 name,keyword,chr(mmk)
      logical check,ltable,lmetab,safe,lneut,loop1,loop2
      real(8) rctbp,rcter,rcfbp,rct,ppp
      integer mxn1,nxn1,idnode,nold
      integer itmols,ksite,numsit,isite,nrept,ifrz,i,j
      integer ishls,ibonds,numcon,numang,icon,iang,idih,numdih
      integer numinv,iinv,numgrp,kgrp,numgsit,numteth,iteth
      integer ipmf,jpmf,npmf,itpvdw,itptbp,itpfbp
      integer itpter,k,nfld,nummols,idum,numshl,nneu
      integer numbonds,itpmet,iii,ngrid
      
      mxtmls=0
      mxatms=0
      mxgrp=0
      mxtcon=0
      mxtbnd=0
      mxtang=0
      mxtdih=0
      mxtinv=0
      mxpmf=0
      mxspmf=0
      mxungp=0
      mxngp=0
      mxneut=0
      mxmet=0
      mxatyp=0
      mxn1=0
      nxn1=0
      nold=-1
      mxgatm=0
      mxteth=0
      msteth=0
      mxvdw=0
      mxtbp=0
      mxter=0
      mxexcl=0
      mxsite=0
      mxbond=0
      mxcons=0
      mxangl=0
      mxdihd=0
      mxinv=0
      mxshl=0
      mxtshl=0
      mxfbp=0
      mxgrid=0 
      rctbp=0.d0
      rcter=0.d0
      rcfbp=0.d0
      safe=.true.
      loop1=.true.
      loop2=.true.
      lneut=.false.
      ltable=.false.
      lmetab=.false.
      
c     open force field data file
      
      if(idnode.eq.0)open (nfield,file='FIELD')
      
      call getrec(safe,idnode,nfield)
      if(.not.safe)call abortscan(52,idnode)
      
c     read and process directives from field file
      
      do while(loop1)
        
        call getrec(safe,idnode,nfield)
        if(.not.safe)call abortscan(52,idnode)
        call lowcase(record,lenrec)
        
        if(findstring('neut',record,idum))then
          
          lneut=.true.
          
        elseif(findstring('molecu',record,idum))then
          
          mxtmls=intstr(record,lenrec,idum)
          
          do itmols=1,mxtmls
            
            loop2=.true.
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            
            do while(loop2)
              
              call getrec(safe,idnode,nfield)
              if(.not.safe)call abortscan(52,idnode)
              call lowcase(record,lenrec)
              
              ksite=0
              
              if(findstring('nummol',record,idum))then
                
                nummols=intstr(record,lenrec,idum)
                
              elseif(findstring('atoms',record,idum))then
                
                numsit=intstr(record,lenrec,idum)
                mxatms=mxatms+numsit*nummols
                mxsite=mxsite+numsit
                ksite=0
                do isite=1,numsit
                  
                  if(ksite.lt.numsit)then
                    
                    call getrec(safe,idnode,nfield)
                    if(.not.safe)call abortscan(52,idnode)
                    
                    call getword(name,record,8,lenrec)
                    ppp=dblstr(record,lenrec,idum)
                    ppp=dblstr(record,lenrec,idum)
                    nrept=intstr(record,lenrec,idum)
                    ifrz=intstr(record,lenrec,idum)
                    nneu=intstr(record,lenrec,idum)
                    if(nrept.eq.0)nrept=1
                    if(lneut)then
                      if(nneu.ne.nold) nxn1=0
                      nxn1=nxn1+nrept
                      mxn1=max(mxn1,nxn1)
                      nold=nneu
                    endif
                    
                    if(mxatyp.eq.0)then
                      
                      mxatyp=1
                      chr(1)=name
                      
                    else
                      
                      check=.true.
                      do j=1,mxatyp
                        
                        if(name.eq.chr(j))check=.false.
                        
                      enddo
                      if(check)then
                        
                        mxatyp=mxatyp+1
                        if(mxatyp.le.mmk)chr(mxatyp)=name
                        
                      endif
                      
                    endif
                    if(nrept.eq.0)nrept=1
                    ksite=ksite+nrept
                    
                  endif
                  
                enddo
                
                if(mmk.lt.mxatyp)call abortscan(34,idnode)
                
                if(lneut)mxneut=mxneut+nneu*nummols
                
              elseif(findstring('shell',record,idum))then
                
                numshl=intstr(record,40,idum)
                mxtshl=mxtshl+numshl
                mxshl=mxshl+nummols*numshl
                
                do ishls=1,numshl
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('bonds',record,idum))then
                
                numbonds=intstr(record,lenrec,idum)
                mxtbnd=mxtbnd+numbonds
                mxbond=mxbond+nummols*numbonds
                
                do ibonds=1,numbonds
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('constr',record,idum))then
                
                numcon=intstr(record,lenrec,idum)
                mxtcon=mxtcon+numcon
                mxcons=mxcons+nummols*numcon
                
                do icon=1,numcon
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('angles',record,idum))then
                
                numang=intstr(record,lenrec,idum)
                mxtang=mxtang+numang
                mxangl=mxangl+nummols*numang
                
                do iang=1,numang
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('dihedr',record,idum))then
                
                numdih=intstr(record,lenrec,idum)
                mxtdih=mxtdih+numdih
                mxdihd=mxdihd+nummols*numdih
                
                do idih=1,numdih
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('invers',record,idum))then
                
                numinv=intstr(record,lenrec,idum)
                mxtinv=mxtinv+numinv
                mxinv=mxinv+nummols*numinv
                
                do iinv=1,numinv
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('rigid',record,idum))then
                
                numgrp=intstr(record,lenrec,idum)
                mxungp=mxungp+numgrp
                mxgrp=mxgrp+numgrp*nummols
                
                do kgrp=1,numgrp
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                  numgsit=intstr(record,lenrec,idum)
                  mxgatm=mxgatm+numgsit*nummols
                  mxngp=max(mxngp,numgsit)
                  do j=1,numgsit
                    
                    iii=intstr(record,lenrec,idum)
                    if(iii.eq.0)then
                      call getrec(safe,idnode,nfield)
                      if(.not.safe)call abortscan(52,idnode)
                      iii=intstr(record,lenrec,idum)
                    endif
                    
                  enddo
                  
                enddo
                
              elseif(findstring('teth',record,idum))then
                
                numteth=intstr(record,lenrec,idum)
                mxteth=mxteth+numteth
                msteth=msteth+numteth*nummols
                
                do iteth=1,numteth
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  
                enddo
                
              elseif(findstring('pmf',record,idum))then
                
                do ipmf=1,2
                  
                  call getrec(safe,idnode,nfield)
                  if(.not.safe)call abortscan(52,idnode)
                  call lowcase(record,lenrec)
                  npmf=intstr(record,lenrec,idum)
                  mxspmf=mxspmf+npmf
                  
                  do jpmf=1,npmf
                    
                    call getrec(safe,idnode,nfield)
                    if(.not.safe)call abortscan(52,idnode)
                    
                  enddo
                  
                enddo
                
                mxpmf=mxpmf+nummols
                
              elseif(findstring('finish',record,idum))then
                
                loop2=.false.
                
              endif
              
            enddo
            
          enddo
          
        elseif(findstring('vdw',record,idum))then
          
          if(findstring('tab',record,idum))ltable=.true.
          mxvdw=intstr(record,lenrec,idum)
          do itpvdw=1,mxvdw
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            call lowcase(record,lenrec)
            if(findstring('tab',record,idum))ltable=.true.
            
          enddo
          mxvdw=max(mxvdw,(mxatyp*(mxatyp+1))/2)
          
          if(ltable)then
            
            if(idnode.eq.0)open(ntable,file='TABLE')
            
            call getrec(safe,idnode,ntable)
            if(.not.safe)call abortscan(24,idnode)
            call getrec(safe,idnode,ntable)
            if(.not.safe)call abortscan(24,idnode)
            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            mxgrid=max(mxgrid,intstr(record,lenrec,idum))
            
            close (ntable)
            
          endif
          
        elseif(findstring('metal',record,idum))then
          
          if(findstring('eam',record,idum))lmetab=.true.
          mxmet=intstr(record,lenrec,idum)
          do itpmet=1,mxmet
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            call lowcase(record,lenrec)
            if(findstring('eam',record,idum))lmetab=.true.
            
          enddo
          mxmet=max(mxmet,(mxatyp*(mxatyp+1))/2)
          
          if(lmetab)then
            
            if(idnode.eq.0)open(ntable,file='TABEAM')
            
            call getrec(safe,idnode,ntable)
            if(.not.safe)call abortscan(24,idnode)
            call getrec(safe,idnode,ntable)
            if(.not.safe)call abortscan(24,idnode)
            do i=1,intstr(record,lenrec,idum)
              
              call getrec(safe,idnode,ntable)
              if(.not.safe)call abortscan(24,idnode)
              ngrid=intstr(record,lenrec,idum)
              mxgrid=max(mxgrid,ngrid+4)
              do j=1,(ngrid+3)/4
                
                call getrec(safe,idnode,ntable)
                if(.not.safe)call abortscan(24,idnode)
                
              enddo
              
            enddo
            
            close (ntable)
            
          endif
          
        elseif(findstring('tbp',record,idum))then
          
          mxtbp=intstr(record,lenrec,idum)
          
          do itptbp=1,mxtbp
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            call getword(name,record,8,lenrec)
            call getword(name,record,8,lenrec)
            call getword(name,record,8,lenrec)
            call getword(keyword,record,4,lenrec)
            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            rct=dblstr(record,lenrec,idum)
            rctbp=max(rctbp,rct)
            
          enddo
          
        elseif(findstring('fbp',record,idum))then
          
          mxfbp=intstr(record,lenrec,idum)
          do itpfbp=1,mxfbp
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            call getword(name,record,8,lenrec)
            call getword(name,record,8,lenrec)
            call getword(name,record,8,lenrec)
            call getword(name,record,8,lenrec)
            call getword(keyword,record,4,lenrec)
            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            rct=dblstr(record,lenrec,idum)
            rcfbp=max(rcfbp,rct)
            
          enddo
          
        elseif(findstring('tersof',record,idum))then
          
          mxter=intstr(record,lenrec,idum)
          
          do itpter=1,mxter
            
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            call getrec(safe,idnode,nfield)
            if(.not.safe)call abortscan(52,idnode)
            rct=dblstr(record,lenrec,idum)
            rcter=max(rcter,rct)
            
          enddo
          
        elseif(findstring('extern',record,idum))then
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)call abortscan(52,idnode)
          nfld=intstr(record,lenrec,idum)
          if(nfld.eq.0)nfld=5
          call getrec(safe,idnode,nfield)
          if(.not.safe)call abortscan(52,idnode)
          
          do k=1,nfld
            
            ppp=dblstr(record,lenrec,idum)
            if(idum.gt.lenrec.and.k.lt.nfld)then
              call getrec(safe,idnode,nfield)
              if(.not.safe)call abortscan(52,idnode)
            endif
            
          enddo
          
        elseif(findstring('close',record,idum))then
          
          loop1=.false.
          
        endif
        
      enddo
      
      if(idnode.eq.0)close (nfield)
      
      if(mxpmf.gt.0)mxpmf=mxatms
      if(mxtcon.gt.0)mxexcl=max(mxexcl,6)
      if(mxtbnd.gt.0)mxexcl=max(mxexcl,6)
      if(mxtang.gt.0)mxexcl=max(mxexcl,16)
      if(mxtdih.gt.0)mxexcl=max(mxexcl,50)
      if(mxtinv.gt.0)mxexcl=max(mxexcl,50)
      if(mxneut.gt.0)mxexcl=max(mxexcl,10*mxn1*mxn1)
      if(mxgrp.gt.0)mxexcl=max(mxexcl,mxngp)
      
      return
      
      end subroutine fldscan
      
      subroutine cfgscan
     x  (idnode,nconf,imcon,volm,xhi,yhi,zhi,cell,buffer)
      
c***********************************************************************
c     
c     dl_poly subroutine for scanning the initial configuration
c     file to determine the number of atoms present
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith  june       1997
c     
c     note: volm is volume containing all particles, not system volume
c     
c***********************************************************************
      
      character*80 header
      character*8 name
      logical lvolm
      real(8) cell,celprp,buffer,extra,volm,xhi,yhi,zhi
      real(8) xxx,yyy,zzz,uuu,vvv,www,coz
      integer idnode,nconf,imcon,i,levcfg
      dimension cell(9),celprp(10),buffer(10),extra(5)
      
      imcon=0
      xhi=0.d0
      yhi=0.d0
      zhi=0.d0
      volm=0.d0
      do i=1,9
        
        cell(i)=0.d0
        
      enddo
      if(idnode.eq.0)then
        
        open (nconf,file='CONFIG')
        
c     read the CONFIG file header
        
        read(nconf,'(a80)',end=100)header
        read(nconf,'(2i10)',end=100)levcfg,imcon
        lvolm=(imcon.eq.0.or.imcon.eq.6)
        
c     specify molecular dynamics simulation cell
        
        if(imcon.gt.0)then
          
          read(nconf,'(3f20.0)',end=100)cell(1),cell(2),cell(3)
          read(nconf,'(3f20.0)',end=100)cell(4),cell(5),cell(6)
          read(nconf,'(3f20.0)',end=100)cell(7),cell(8),cell(9)
          call dcell(cell,celprp)
          
        endif
        
        if(.not.lvolm)then
          
          volm=celprp(10)
          
          if(imcon.eq.4)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.5)then
            
            volm=0.5d0*celprp(10)
            
          elseif(imcon.eq.7)then
            
            volm=0.5d0*celprp(10)
            
          endif
          
        endif
        
        i=0
        do while(.true.)
          
          i=i+1
          if(levcfg.eq.0)then
            
            read(nconf,'(a8)',end=100) name
            read(nconf,'(3f20.0)')xxx,yyy,zzz
            
          else if(levcfg.eq.1)then
            
            read(nconf,'(a8)',end=100) name
            read(nconf,'(3f20.0)')xxx,yyy,zzz
            read(nconf,'(3f20.0)')uuu,vvv,www
            
          else
            
            read(nconf,'(a8)',end=100) name
            read(nconf,'(3f20.0)')xxx,yyy,zzz
            read(nconf,'(3f20.0)')uuu,vvv,www
            read(nconf,'(3f20.0)')uuu,vvv,www
            
          endif
          
          if(lvolm)then
            
            if(i.eq.1)then
              
              xhi=abs(xxx)
              yhi=abs(yyy)
              zhi=abs(zzz)
              
            else
              
              xhi=max(xhi,abs(xxx))
              yhi=max(yhi,abs(yyy))
              zhi=max(zhi,abs(zzz))
              
            endif
            
          endif
          
        enddo
        
  100   continue
        
        if(imcon.eq.0)then
          
          volm=8.d0*xhi*yhi*zhi
          
        else if(imcon.eq.6)then
          
          coz=(cell(1)*cell(4)+cell(2)*cell(5)+cell(3)*cell(6))/
     x      (celprp(1)*celprp(2))
          volm=2.d0*zhi*celprp(1)*celprp(2)*sqrt(1.d0-coz**2)
          
        endif
        
        close (nconf)
        
      endif
      
      extra(1)=dble(imcon)
      extra(2)=xhi
      extra(3)=yhi
      extra(4)=zhi
      extra(5)=volm
      call gdsum(extra,5,buffer)
      call gdsum(cell,9,buffer)
      imcon=nint(extra(1))
      xhi=extra(2)
      yhi=extra(3)
      zhi=extra(4)
      volm=extra(5)
      
      return
      
      end subroutine cfgscan
      
      subroutine conscan
     x  (redirect,lewald,lspme,lhke,nolink,lcshft,lsolva,lfree,lfrmas,
     x  lghost,idnode,imcon,nhko,rcut,rvdw,delr,densvar,drdf,dzdn,
     x  zlen,cell)
      
c***********************************************************************
c     
c     dl_poly subroutine for scanning the contents of the control file
c     
c     copyright - daresbury laboratory 1997
c     author    - w. smith  june       1997
c     
c***********************************************************************
      
      logical safe,lewald,lspme,lhke,peek,nolink,lcshft,lmetad
      logical lsolva,lfree,lfrmas,lghost,redirect
      real(8) cell,celprp,rcut,rvdw,delr,eps,alpha,fac,tol,tol1
      real(8) densvar,drdf,dzdn,zlen
      integer nhko,idnode,imcon,idum,jmp
      integer nlatt,kmax1,kmax2,kmax3,kmaxpow2
      dimension celprp(10),cell(9)
      
      nhko=0
      mxstak=0
      kmaxa=0
      kmaxb=1
      kmaxc=1
      kmaxd=1
      kmaxe=1
      kmaxf=1
      rcut=0.d0
      rvdw=0.d0
      delr=0.d0
      drdf=0.05d0
      dzdn=0.05d0
      zlen=0.d0
      densvar=1.d2
      peek=.true.
      lhke=.false.
      lspme=.false.
      lewald=.false.
      lcshft=.false.
      nolink=.false.
      lghost=.false.
      lfree=.false.
      lfrmas=.false.
      lsolva=.false.
      lmetad=.false.
      redirect=.false.
      
c     open the simulation input file
      
      if(idnode.eq.0)open (nread,file='CONTROL')
      
      call getrec(safe,idnode,nread)
      if(.not.safe)call abortscan(17,idnode)
      
      do while(peek)
        
        call getrec(safe,idnode,nread)
        if(.not.safe)call abortscan(17,idnode)
        call lowcase(record,lenrec)
        if(record(1).ne.'#')then
          
          if(findstring('stack',record,idum))then
            
            mxstak=intstr(record,lenrec,idum)
            
          elseif(findstring('no link',record,idum))then
            
            nolink=.true.
            
          elseif(findstring('metafreeze',record,idum))then
            
            lmetad=.true.
            do while(lmetad)
              call getrec(safe,idnode,nread)
              if(.not.safe)call abortscan(17,idnode)
              call lowcase(record,lenrec)
              lmetad=.not.findstring('endmet',record,idum)
            enddo
            
          elseif(findstring('redirect',record,idum))then
            
            redirect=.true.
            
          elseif(findstring('densvar',record,idum))then
            
            densvar=dblstr(record,lenrec,idum)
            
          elseif(findstring('shift',record,idum).or.
     x        findstring('reaction',record,idum))then
            
            lcshft=.true.
            
          elseif(findstring('ewald',record,idum).or.
     x        findstring('spme',record,idum).or.
     x        findstring('hke',record,idum))then
            
c     read Ewald or HK-Ewald or SPM-Ewald sum parameters
            
            lhke=findstring('hke',record,idum)
            lspme=findstring('spme',record,idum)
            lewald=findstring('ewald',record,idum)
            
            if(findstring('precision',record,idum))then
              
              eps=dblstr(record,lenrec,idum)
              if(lhke) then
                
                nhko=intstr(record,lenrec,idum)
                nlatt=intstr(record,lenrec,idum)
                nlatt=min(nlatt,2)
                
              endif
              
              if(rcut.lt.1.d-6)rcut=10.d0
              
c     compute alpha and the kmax
              
              if(lewald.or.lspme)then
                
                call dcell(cell,celprp)
                eps=min(abs(eps),0.5d0)
                tol=sqrt(abs(log(eps*rcut)))
                alpha=sqrt(abs(log(eps*rcut*tol)))/rcut
                tol1=sqrt(-log(eps*rcut*(2.d0*tol*alpha)**2))
                fac=1.d0
                if(imcon.eq.4.or.imcon.eq.5.or.imcon.eq.7) 
     x            fac=2.d0**(1.d0/3.d0)
                kmax1=nint(0.25d0+fac*celprp(1)*alpha*tol1/pi)
                kmax2=nint(0.25d0+fac*celprp(2)*alpha*tol1/pi)
                kmax3=nint(0.25d0+fac*celprp(3)*alpha*tol1/pi)
                
              elseif(lhke)then
                
                if(nhko.eq.0)then
                  if(eps.le.1.d-6)then
                    alpha=3.46d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=3.14d0/rcut
                  else
                    alpha=2.76d0/rcut
                  endif
                elseif(nhko.eq.1)then
                  if(eps.le.1.d-6)then
                    alpha=4.37d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=4.08d0/rcut
                  else
                    alpha=3.75d0/rcut
                  endif                
                elseif(nhko.eq.2)then
                  if(eps.le.1.d-6)then
                    alpha=5.01d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=4.74d0/rcut
                  else
                    alpha=4.44d0/rcut
                  endif
                elseif(nhko.eq.3)then
                  if(eps.le.1.d-6)then
                    alpha=5.55d0/rcut
                  elseif(eps.le.1.d-5)then
                    alpha=5.28d0/rcut
                  else
                    alpha=5.00d0/rcut
                  endif
                endif
                alpha=alpha/dble(2*nlatt+1)
                if(abs(cell(9)).lt.1.d-8)cell(9)=1.d0
                call dcell(cell,celprp)
                tol=2.d0*alpha*sqrt(abs(log(eps*alpha)))
                tol1=2.d0*alpha*sqrt(abs(log(eps*alpha*tol)))
                kmax1=nint(0.25d0+0.5d0*celprp(1)*tol1/pi)
                kmax2=nint(0.25d0+0.5d0*celprp(2)*tol1/pi)
                kmax3=1
                
              endif
              
            else
              
              alpha=dblstr(record,lenrec,idum)
              kmax1=intstr(record,lenrec,idum)
              kmax2=intstr(record,lenrec,idum)
              
              if(lhke)then
                
                kmax3=1
                nhko=intstr(record,lenrec,idum)
                
              else
                
                kmax3=intstr(record,lenrec,idum)
                
              endif
              
            endif
            
c     for spme double kmax and set to next power of 2, with current
c     upper limit of 512
            
            if(lspme)then
              
              kmaxpow2=1
              do while (kmax1.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2=kmaxpow2 * 2
              end do
              kmaxd=2 * kmaxpow2
              
              kmaxpow2=1
              do while (kmax2.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2=kmaxpow2 * 2
              end do
              kmaxe=2 * kmaxpow2
              
              kmaxpow2=1
              do while (kmax3.gt.kmaxpow2.and.kmaxpow2.lt.256)
                kmaxpow2=kmaxpow2 * 2
              end do
              kmaxf=2 * kmaxpow2
              
            elseif(lhke) then
              
              kmaxa=kmax1
              kmaxb=kmax2
              kmaxc=1
              
            else
              
              kmaxa=kmax1
              kmaxb=kmax2
              kmaxc=kmax3
              
            endif
            
          elseif(findstring('cut',record,idum))then
            
            rcut=dblstr(record,lenrec,idum)
            
          elseif(findstring('rvdw',record,idum))then
            
            rvdw=dblstr(record,lenrec,idum)
            
          elseif(findstring('delr',record,idum))then
            
            delr=dblstr(record,100,idum)
            
          else if(findstring('rdf',record,idum))then
            
            if(.not.findstring('print',record,idum))then
              
              jmp=intstr(record,lenrec,idum)
              drdf=dblstr(record,lenrec,idum)
              
            endif
            
          else if(findstring('zden',record,idum))then
            
            jmp=intstr(record,lenrec,idum)
            dzdn=dblstr(record,lenrec,idum)
            zlen=dblstr(record,lenrec,idum)
            if(dzdn.lt.1.d-8)then
              
              dzdn=0.1d0
              zlen=0.1d0*dble(128)
              
            elseif(zlen.lt.1.d-8)then
              
              zlen=dzdn*dble(128)
              
            endif
            
          elseif(findstring('solva',record,idum))then
            
            lsolva=.true.
            
          elseif(findstring('decomp',record,idum))then
            
            lsolva=.true.
            
          elseif(findstring('free',record,idum))then
            
            lfree=.true.
            
          elseif(findstring('excit',record,idum))then
            
            lghost=.true.
            lsolva=.true.
            
          elseif(findstring('reset_mass',record,idum))then
            
            lfrmas=.true.
            
          elseif(findstring('switch',record,idum))then
            
            lghost=.true.
            lsolva=.true.
            
          elseif(findstring('finish',record,idum))then
            
            peek=.false.
            
          endif
          
        endif
        
      enddo
      
      if(idnode.eq.0)close (nread)
      if(abs(rvdw).le.1.d-10)rvdw=rcut
      if(drdf.lt.1.d-8)drdf=0.05d0
      if(dzdn.lt.1.d-8)dzdn=0.05d0
      
      return
      
      end subroutine conscan
      
      subroutine abortscan(key,idnode)
      
c*********************************************************************
c     
c     dl_poly subroutine for controlled exit of file scan
c     
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c*********************************************************************
      
      integer key,idnode
      
      write(nrite,'(/,/,1x,a,i5)') 
     x  'DL_POLY terminated due to error ', key
      
      if(key.eq.17)then
        
        write(nrite,'(/,/,1x,a)')
     x    'error - strange exit from CONTROL file processing'
        
      else if(key.eq.52)then
        
        write(nrite,'(/,/,1x,a)')
     x    'error - end of FIELD file encountered'
        
      else if(key.eq.24)then
        
        write(nrite,'(/,/,1x,a)')
     x    'error - end of file encountered in TABLE file'
        
      else if(key.eq.34)then
        
        write(nrite,'(/,/,1x,a)')
     x    'error - character array memory allocation failure'
        
      endif
      
      if(idnode.eq.0) then
        close (nrite)
        close (nhist)
        close (nread)
        close (nconf)
        close (nstats)
        close (nrest)
        close (nfield)
        close (ntable)
      endif
      
      call gsync()
      call exitcomms()
      
      return
      end subroutine abortscan
      
      subroutine dcell(aaa,bbb)
      
c***********************************************************************
c     
c     dl_poly subroutine to calculate the dimensional properties of
c     a simulation cell specified by the input matrix aaa.
c     the results are returned in the array bbb, with :
c     
c     bbb(1 to 3) - lengths of cell vectors
c     bbb(4 to 6) - cosines of cell angles
c     bbb(7 to 9) - perpendicular cell widths
c     bbb(10)     - cell volume
c     
c     copyright daresbury laboratory 1992
c     author - w. smith         july 1992
c     
c***********************************************************************
      
      real(8) aaa,bbb,axb1,axb2,axb3,bxc1,bxc2,bxc3,cxa1,cxa2,cxa3
      
      dimension aaa(9),bbb(10)
      
c     calculate lengths of cell vectors
      
      bbb(1)=sqrt(aaa(1)*aaa(1)+aaa(2)*aaa(2)+aaa(3)*aaa(3))
      bbb(2)=sqrt(aaa(4)*aaa(4)+aaa(5)*aaa(5)+aaa(6)*aaa(6))
      bbb(3)=sqrt(aaa(7)*aaa(7)+aaa(8)*aaa(8)+aaa(9)*aaa(9))
      
c     calculate cosines of cell angles
      
      bbb(4)=(aaa(1)*aaa(4)+aaa(2)*aaa(5)+aaa(3)*aaa(6))/(bbb(1)*bbb(2))
      bbb(5)=(aaa(1)*aaa(7)+aaa(2)*aaa(8)+aaa(3)*aaa(9))/(bbb(1)*bbb(3))
      bbb(6)=(aaa(4)*aaa(7)+aaa(5)*aaa(8)+aaa(6)*aaa(9))/(bbb(2)*bbb(3))
      
c     calculate vector products of cell vectors
      
      axb1=aaa(2)*aaa(6)-aaa(3)*aaa(5)
      axb2=aaa(3)*aaa(4)-aaa(1)*aaa(6)
      axb3=aaa(1)*aaa(5)-aaa(2)*aaa(4)
      bxc1=aaa(5)*aaa(9)-aaa(6)*aaa(8)
      bxc2=aaa(6)*aaa(7)-aaa(4)*aaa(9)
      bxc3=aaa(4)*aaa(8)-aaa(5)*aaa(7)
      cxa1=aaa(8)*aaa(3)-aaa(2)*aaa(9)
      cxa2=aaa(1)*aaa(9)-aaa(3)*aaa(7)
      cxa3=aaa(2)*aaa(7)-aaa(1)*aaa(8)
      
c     calculate volume of cell
      
      bbb(10)=abs(aaa(1)*bxc1+aaa(2)*bxc2+aaa(3)*bxc3)
      
c     calculate cell perpendicular widths
      
      bbb(7)=bbb(10)/sqrt(bxc1*bxc1+bxc2*bxc2+bxc3*bxc3)
      bbb(8)=bbb(10)/sqrt(cxa1*cxa1+cxa2*cxa2+cxa3*cxa3)
      bbb(9)=bbb(10)/sqrt(axb1*axb1+axb2*axb2+axb3*axb3)
      
      return
      end subroutine dcell
      
      end module setup_module

      module dlp_io_module
c***********************************************************************
c     
c     DL_POLY routines needed for minimalistic I/O for
c     CONFIG, FIELD & HISTORY files (called by wrappers to VOTCA calls) 
c     
c     copyright - daresbury laboratory  2013
c     author    - Andrey Brukhno   July 2013
c
c     some routines are based on the corresponding code of
c     DL_POLY Classic by W. Smith and T.R. Forester
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

!      use error_module
!      use parse_module
      use setup_module

      implicit none

c AB: maximum length of records/string is defined by lenrec (parse_module.f)
c AB: set some useful length parameters for CHARACTER types

      INTEGER, PARAMETER :: half = 8
      INTEGER, PARAMETER :: word = 16
      INTEGER, PARAMETER :: srec = 80
      INTEGER, PARAMETER :: lrec = 255
      INTEGER, PARAMETER :: lstr = 800

c AB: set the precision parameters for INTEGER and REAL types

!      INTEGER, PARAMETER :: int = SELECTED_REAL_KIND(5)
!      INTEGER, PARAMETER :: lng = SELECTED_REAL_KIND(10)

!      INTEGER, PARAMETER :: flt = SELECTED_REAL_KIND(6,38)
!      INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(15,308)

      character(len=srec), allocatable, target :: mtmolname(:)
      character(len=half), allocatable, target :: mtatmname(:)
      character(len=half), allocatable, target :: mtatmtype(:)
      character(len=half), allocatable, target :: fratmname(:)
      character(len=half), allocatable, target :: fratmtype(:)

      real(8), allocatable, target :: xyz(:,:)
      real(8), allocatable, target :: vel(:,:)
      real(8), allocatable, target :: frc(:,:)

      real(8), allocatable, target :: mtatmmass(:)
      real(8), allocatable, target :: mtatmchrg(:)
      real(8), allocatable, target :: fratmmass(:)
      real(8), allocatable, target :: fratmchrg(:)

      integer, allocatable, target :: mtnmolecs(:)
      integer, allocatable, target :: mtnmolatm(:)
      integer, allocatable, target :: mtnmolgrp(:)
      integer, allocatable, target :: mtatmnrep(:)
      integer, allocatable, target :: mtatmngrp(:)

      integer :: imemall,imemfld,imemfrm

      
      TYPE, bind(C) :: FullCellT

      real(kind=C_DOUBLE) :: uc(9)
      real(kind=C_DOUBLE) :: rc(9)
      real(kind=C_DOUBLE) :: ac(9)
      real(kind=C_DOUBLE) :: pc(10)
   
      END TYPE FullCellT
      

      TYPE, bind(C) :: FrameSiteT

      character(len=1,kind=c_char) :: name(9)
      character(len=1,kind=c_char) :: type(9)

      real(kind=C_DOUBLE) :: r(3)
      real(kind=C_DOUBLE) :: v(3)
      real(kind=C_DOUBLE) :: f(3)

      real(kind=C_DOUBLE) :: m,q
      integer(kind=C_INT) :: id,im,ig

      END TYPE FrameSiteT
      

      TYPE, bind(C) :: FrameSpecsT

      character(len=1,kind=c_char) :: title(81)
      real(kind=C_DOUBLE)    :: cell(9)
      real(kind=C_DOUBLE)    :: tstep,energy
      integer(kind=C_INT)    :: nstep,nsites,imcon,keytrj

      END TYPE FrameSpecsT


      TYPE, bind(C) :: FieldSiteT

      character(len=1,kind=c_char) :: name(9)
      character(len=1,kind=c_char) :: type(9)

      real(kind=C_DOUBLE) :: m,q
      integer(kind=C_INT) :: idmol,idgrp,ifrzn,nrept

      END TYPE FieldSiteT


      TYPE, bind(C) :: MolecSpecsT

        character(len=1,kind=c_char) :: name(81)
        integer(kind=C_INT)    :: id,nrept,nsites,ngroups

      END TYPE MolecSpecsT


      TYPE, bind(C) :: FieldSpecsT

        character(len=1,kind=c_char) :: title(81)
        character(len=1,kind=c_char) :: units(9)

        integer(kind=C_INT) :: ineut
        integer(kind=C_INT) :: nmols,natms
        integer(kind=C_INT) :: nbonds,nangles,ndhdrs,nexcls

      END TYPE FieldSpecsT


      contains


      subroutine IniFullCell(Cell)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FullCellT), pointer :: Cell

      Cell%uc(1) = 0.0_C_DOUBLE
      Cell%uc(2) = 0.0_C_DOUBLE
      Cell%uc(3) = 0.0_C_DOUBLE
      Cell%uc(4) = 0.0_C_DOUBLE
      Cell%uc(5) = 0.0_C_DOUBLE
      Cell%uc(6) = 0.0_C_DOUBLE
      Cell%uc(7) = 0.0_C_DOUBLE
      Cell%uc(8) = 0.0_C_DOUBLE
      Cell%uc(9) = 0.0_C_DOUBLE

      Cell%rc(1) = 0.0_C_DOUBLE
      Cell%rc(2) = 0.0_C_DOUBLE
      Cell%rc(3) = 0.0_C_DOUBLE
      Cell%rc(4) = 0.0_C_DOUBLE
      Cell%rc(5) = 0.0_C_DOUBLE
      Cell%rc(6) = 0.0_C_DOUBLE
      Cell%rc(7) = 0.0_C_DOUBLE
      Cell%rc(8) = 0.0_C_DOUBLE
      Cell%rc(9) = 0.0_C_DOUBLE

      Cell%ac(1) = 0.0_C_DOUBLE
      Cell%ac(2) = 0.0_C_DOUBLE
      Cell%ac(3) = 0.0_C_DOUBLE
      Cell%ac(4) = 0.0_C_DOUBLE
      Cell%ac(5) = 0.0_C_DOUBLE
      Cell%ac(6) = 0.0_C_DOUBLE
      Cell%ac(7) = 0.0_C_DOUBLE
      Cell%ac(8) = 0.0_C_DOUBLE
      Cell%ac(9) = 0.0_C_DOUBLE

      Cell%pc(1) = 0.0_C_DOUBLE
      Cell%pc(2) = 0.0_C_DOUBLE
      Cell%pc(3) = 0.0_C_DOUBLE
      Cell%pc(4) = 0.0_C_DOUBLE
      Cell%pc(5) = 0.0_C_DOUBLE
      Cell%pc(6) = 0.0_C_DOUBLE
      Cell%pc(7) = 0.0_C_DOUBLE
      Cell%pc(8) = 0.0_C_DOUBLE
      Cell%pc(9) = 0.0_C_DOUBLE
      Cell%pc(10)= 0.0_C_DOUBLE

      return
      end subroutine IniFullCell


      subroutine IniCell(Cell)

      use iso_c_binding!, only : C_DOUBLE

      implicit none

      real(C_DOUBLE) :: Cell(9)

      Cell(1) = 0.0_C_DOUBLE
      Cell(2) = 0.0_C_DOUBLE
      Cell(3) = 0.0_C_DOUBLE
      Cell(4) = 0.0_C_DOUBLE
      Cell(5) = 0.0_C_DOUBLE
      Cell(6) = 0.0_C_DOUBLE
      Cell(7) = 0.0_C_DOUBLE
      Cell(8) = 0.0_C_DOUBLE
      Cell(9) = 0.0_C_DOUBLE

      return
      end subroutine IniCell

      subroutine SetCell(Cell,Cell0)

      use iso_c_binding!, only : C_DOUBLE

      implicit none

      real(C_DOUBLE) :: Cell(9)
      real(8) :: Cell0(9)

      Cell(1) = real(Cell0(1),kind=C_DOUBLE)
      Cell(2) = real(Cell0(2),kind=C_DOUBLE)
      Cell(3) = real(Cell0(3),kind=C_DOUBLE)
      Cell(4) = real(Cell0(4),kind=C_DOUBLE)
      Cell(5) = real(Cell0(5),kind=C_DOUBLE)
      Cell(6) = real(Cell0(6),kind=C_DOUBLE)
      Cell(7) = real(Cell0(7),kind=C_DOUBLE)
      Cell(8) = real(Cell0(8),kind=C_DOUBLE)
      Cell(9) = real(Cell0(9),kind=C_DOUBLE)

      return
      end subroutine SetCell


      subroutine IniFrmSite(FrmSite)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FrameSiteT) :: FrmSite

      call inichar(FrmSite%name,9)
      call inichar(FrmSite%type,9)

      FrmSite%name(1) = C_NULL_CHAR
      FrmSite%type(1) = C_NULL_CHAR

      FrmSite%r(1) = 0.0_C_DOUBLE
      FrmSite%r(2) = 0.0_C_DOUBLE
      FrmSite%r(3) = 0.0_C_DOUBLE

      FrmSite%v(1) = 0.0_C_DOUBLE
      FrmSite%v(2) = 0.0_C_DOUBLE
      FrmSite%v(3) = 0.0_C_DOUBLE

      FrmSite%f(1) = 0.0_C_DOUBLE
      FrmSite%f(2) = 0.0_C_DOUBLE
      FrmSite%f(3) = 0.0_C_DOUBLE

      FrmSite%m = 0.0_C_DOUBLE
      FrmSite%q = 0.0_C_DOUBLE

      FrmSite%id = 0_c_int
      FrmSite%im = 0_c_int
      FrmSite%ig = 0_c_int

      return
      end subroutine IniFrmSite


      subroutine IniFrmSites(FrmSite,ns)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FrameSiteT), pointer :: FrmSite(:)

      integer ns,is

      do is=1,ns

      call inichar(FrmSite(is)%name,9)
      call inichar(FrmSite(is)%type,9)

      FrmSite(is)%name(1) = C_NULL_CHAR
      FrmSite(is)%type(1) = C_NULL_CHAR

      FrmSite(is)%r(1) = 0.0_C_DOUBLE
      FrmSite(is)%r(2) = 0.0_C_DOUBLE
      FrmSite(is)%r(3) = 0.0_C_DOUBLE

      FrmSite(is)%v(1) = 0.0_C_DOUBLE
      FrmSite(is)%v(2) = 0.0_C_DOUBLE
      FrmSite(is)%v(3) = 0.0_C_DOUBLE

      FrmSite(is)%f(1) = 0.0_C_DOUBLE
      FrmSite(is)%f(2) = 0.0_C_DOUBLE
      FrmSite(is)%f(3) = 0.0_C_DOUBLE

      FrmSite(is)%m = 0.0_C_DOUBLE
      FrmSite(is)%q = 0.0_C_DOUBLE

      FrmSite(is)%id = 0_c_int
      FrmSite(is)%im = 0_c_int
      FrmSite(is)%ig = 0_c_int

      end do

      return
      end subroutine IniFrmSites


      subroutine IniFrmSpecs(FrmBase)

      TYPE(FrameSpecsT),pointer :: FrmBase

      call inichar(FrmBase%title,81)
      FrmBase%title(1) = C_NULL_CHAR

      call IniCell(FrmBase%cell)

      FrmBase%tstep  = 0.0_c_double
      FrmBase%energy = 0.0_c_double
      FrmBase%nstep  = 0_c_int
      FrmBase%nsites = 0_c_int
      FrmBase%imcon  = 0_c_int
      FrmBase%keytrj = 0_c_int

      return
      end subroutine IniFrmSpecs


      subroutine IniFldSite(FldSite)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FieldSiteT) :: FldSite

      call inichar(FldSite%name,9)
      call inichar(FldSite%type,9)

      FldSite%name(1) = C_NULL_CHAR
      FldSite%type(1) = C_NULL_CHAR

      FldSite%m = 0.0_C_DOUBLE
      FldSite%q = 0.0_C_DOUBLE

      FldSite%idmol   = 0_c_int
      FldSite%idgrp   = 0_c_int
      FldSite%ifrzn   = 0_c_int
      FldSite%nrept   = 0_c_int

      return
      end subroutine IniFldSite


      subroutine IniFldSites(FldSite,ns)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FieldSiteT), pointer :: FldSite(:)

      integer ns,is

      do is=1,ns

      call inichar(FldSite(is)%name,9)
      call inichar(FldSite(is)%type,9)

      FldSite(is)%name(1) = C_NULL_CHAR
      FldSite(is)%type(1) = C_NULL_CHAR

      FldSite(is)%m = 0.0_C_DOUBLE
      FldSite(is)%q = 0.0_C_DOUBLE

      FldSite(is)%idmol   = 0_c_int
      FldSite(is)%idgrp   = 0_c_int
      FldSite(is)%ifrzn   = 0_c_int
      FldSite(is)%nrept   = 0_c_int

      end do

      return
      end subroutine IniFldSites


      subroutine IniMolSpec(MolBase)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(MolecSpecsT) :: MolBase

      call inichar(MolBase%name,81)
      MolBase%name(1) = C_NULL_CHAR

      MolBase%id      = 0_c_int
      MolBase%nrept   = 0_c_int
      MolBase%nsites  = 0_c_int
      MolBase%ngroups = 0_c_int

      return
      end subroutine IniMolSpec


      subroutine IniMolSpecs(MolBase,ns)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(MolecSpecsT), pointer :: MolBase(:)

      integer ns,is

      do is=1,ns

      call inichar(MolBase(is)%name,81)
      MolBase(is)%name(1) = C_NULL_CHAR

      MolBase(is)%id      = 0_c_int
      MolBase(is)%nrept   = 0_c_int
      MolBase(is)%nsites  = 0_c_int
      MolBase(is)%ngroups = 0_c_int

      end do

      return
      end subroutine IniMolSpecs


      subroutine IniFldSpecs(FldBase)

      use iso_c_binding!, only : C_NULL_CHAR

      implicit none

      TYPE(FieldSpecsT), pointer :: FldBase

      call inichar(FldBase%title,81)
      call inichar(FldBase%units,9)

      FldBase%title(1) = C_NULL_CHAR
      FldBase%units(1) = C_NULL_CHAR

      FldBase%ineut   = 0_c_int
      FldBase%nmols   = 0_c_int
      FldBase%natms   = 0_c_int
      FldBase%nbonds  = 0_c_int
      FldBase%nangles = 0_c_int
      FldBase%ndhdrs  = 0_c_int
      FldBase%nexcls  = 0_c_int

      return
      end subroutine IniFldSpecs


c***********************************************************************
c     
c     DL_POLY routines to write CONFIG, FIELD & HISTORY files (for CG-mapping) 
c     NOTE: maximum length of strings is defined by lenrec (parse_module.f)
c     
c     copyright daresbury laboratory   2013
c     author - Andrey Brukhno May-June 2013
c     
c***********************************************************************
      
      subroutine fldread(idnode,fname,
     x header,units,lneut,nmolt,nmols,natms,istate)

c***********************************************************************
c     
c     dl_poly routine for scanning the field file to determine the
c     required parameters (parse_module.f)
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith  november   1994
c
c     amended  by Andrey Brukhno   May 2013
c     
c***********************************************************************
      use iso_c_binding
!      use dlp_io_module

      implicit none

      integer, parameter :: mmk=1000
      
c AB: here the pointers may be associated with native Fortran (derived) types
c AB: but these are mostly for interoperability with structures allocated in C(++)

!      TYPE(FieldSpecsT), pointer :: FldBase
!      TYPE(MolecSpecsT), pointer :: MolBase(:)
!      TYPE(FieldSiteT),  pointer :: FldSite(:)

!      TYPE(FullCellT),   target  :: SysCellF
!      TYPE(FrameSpecsT), target  :: FrmBaseF
!      TYPE(FrameSiteT),  allocatable, target :: FrmSiteF(:)

!      TYPE(FullCellT),   pointer :: SysCell
!      TYPE(FrameSpecsT), pointer :: FrmBase
!      TYPE(FrameSiteT),  pointer :: FrmSite(:)

c AB: these are now global within the module scope

!      character*(80) mtmolname(nmolt)
!      character*(8)  mtatmname(natms)
!      real(8)  mtatmmass(natms),mtatmcharge(natms)
!      integer mtatmnrep(natms),mtatmngrp(natms)
!      integer mtnmolecs(nmolt)
!      integer mtnmolatm(nmolt)
!      integer mtnmolgrp(nmolt)

      character(len=*) :: fname

      character(len=8)  :: chr(mmk)

      character(len=80) :: header,lname1,lname2,recstr
      character(len=8)  :: name,type,keyword,units

      real(8) :: atmmass,atmcharge
      real(8) :: rctbp,rcter,rcfbp,rct,ppp

      integer :: idnode,nmolt,nmols,natms,istate,ios

      integer :: mxn1,nxn1,nold,lunits,lchr,i,j,ich
      integer :: itmols,jtmols,ksite,numsit,isite,jsite,nrept,ifrz
      integer :: ishls,ibonds,numcon,numang,icon,iang,idih,numdih
      integer :: numinv,iinv,numgrp,kgrp,numgsit,numteth,iteth
      integer :: ipmf,jpmf,npmf,itpvdw,itptbp,itpfbp
      integer :: itpter,k,nfld,nummols,numatms,idum,numshl
      integer :: numbonds,itpmet,iii,ngrid,nneu,ngrp

      logical :: lneut,safe,loop1,loop2,check,ltable,lmetab
      
c AB: allocate all the needed big arrays (see above; maybe done externally!)

!      if( .not.memall(.true.) ) then
!        write(*,*)'---'
!        write(*,*)'fldread(): Failed allocating memory - FULL STOP!'
!        write(*,*)'==='
!        stop
!        return
!      endif

c AB: initialise all the needed pointers/arrays (see above; maybe done externally!)

!      call IniFldSpecs(FldBase)
!      call IniMolSpecs(MolBase,mxtmls)
!      call IniFldSite(FldSite,mxsite)

!      do im = 1,nmolt
!        call IniMolSpec(MolBase(im))
!      end do

!      do ia = 1,natms
!        call IniFldSite(FldSite(ia))
!      end do

c AB: below is just a test for pointer-based initialisation of HISTORY-related types

!      nullify(SysCell,FrmBase,FrmSite)
!      Allocate(FrmSiteF(mxatms))

!      SysCell => SysCellF
!      FrmBase => FrmBaseF
!      FrmSite => FrmSiteF

!      call IniFullCell(SysCell)
!      call IniFrameSpecs(FrmBase)
!      call IniFrmSites(FrmSite,mxatms)

!      do isite = 1,mxatms
!        call IniFrmSite(FrmSite(isite))
!      end do

      istate = 0
      ios    = 0

      lunits = 0
      nmols  = 0
      lchr   = 0

      numatms=0

!      mxtmls=0
!      mxatms=0
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
!      mxsite=0
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
      
      if(idnode.eq.0) open (nfield,file=trim(fname),iostat=ios)
      
      if( ios.ne.0 ) then
        write(*,'(/,3a)')
     x  'fldread(): ERROR: FIELD (topology) file "',trim(fname),
     x  '" could not be open, escaping!'

        istate = -3

        return
      endif

      call inichar(record,lenrec)

c AB: read in the very first record from FIELD file
      call getrec(safe,idnode,nfield)
!      if(.not.safe) call abortscan(52,idnode)
      if(.not.safe) then 
        istate = -5
        return
      endif

      call strip(record,lenrec)

c AB: it must be the title of the simulation task
      lchr = lenchar(record,80)
      call cpchr2str(record,header,lchr)

ccc AB: debugging only (comment out)
!      write(*,*)
!      write(*,*)'FIELD scan: first line reads',lchr
!      write(*,*)
!      write(*,*)'"',(record(ich),ich=1,lchr),'"'
ccc
c     read and process directives from field file
      
      do while(loop1)

        call inichar(record,lenrec)

        call getrec(safe,idnode,nfield)
!        if(.not.safe)call abortscan(52,idnode)
        if(.not.safe) then 
          istate = -5
          return
        endif

        call lowcase(record,lenrec)

c AB: looking for the units directive
        if(findstring('unit',record,idum))then
          
ccc AB: next is for testing (comment out)
!          write(*,'(a,80a)')' FIELD scan: found "units" directive - ',
!     *    (record(ich),ich=1,80)
ccc
c AB: checking the units definition against the known choices
          if(findstring('ev',record,idum))then
             units='eV'
             lunits=1
          elseif(findstring('kcal',record,idum))then
             units='kcal'
             lunits=2
          elseif(findstring('kj',record,idum))then
             units='kJ/mol'
             lunits=3
          elseif(findstring('k',record,idum))then
             units='K'
             lunits=4
          elseif(findstring('int',record,idum))then
             units='internal'
             lunits=0
          else
!             call abortscan(999,idnode)
!        if(.not.safe) then 
             istate = -4
             return
!        endif
          endif

          call stripstr(units)

!          lchr = min(8,len(trim(units)))
!          call cpstr2chr(units,FldBase%units,lchr)
!          FldBase%units(lchr+1) = C_NULL_CHAR

ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found "units" directive - "',
!     *         trim(units)
!!     *        ,'" energy units =?= "',(FldBase%units(ich),ich=1,lchr),'"'
ccc       
        elseif(findstring('neut',record,idum))then
                    
ccc AB: next is for testing (comment out)
!        write(*,'(a,80a)')' FIELD scan: found "neut" directive - ',
!     *  (record(ich),ich=1,80)
ccc
          lneut=.true.
                    
!          FldBase%ineut = 1

        elseif(findstring('molecul',record,idum))then
          
!          mxtmls=intstr(record,lenrec,idum)
!
ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found "molecules" directive - ',mxtmls
!
!          if(mxtmls.lt.1 .or. mxtmls.gt.maxtmol) then
!             call abortscan(998,idnode)
!          endif
!
!          nmolt = mxtmls
ccc
          nmolt = intstr(record,lenrec,idum)

          if(nmolt.lt.1 .or. nmolt.gt.mxtmls) then
!            call abortscan(998,idnode)
!        if(.not.safe) then 
            istate = -4
            return
!        endif
          endif
          
!          FldBase%nmols = int(nmolt,kind=C_INT)

!          do itmols=1,mxtmls
          do itmols=1,nmolt
            
            loop2=.true.
            call getrec(safe,idnode,nfield)
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
              istate = -5
              return
            endif

c AB: store the molecule name in mtmolname(...) array (max 80 chars)
            call strip(record,lenrec)

            lchr = lenchar(record,80)
            call cpchr2str(record,mtmolname(itmols),lchr)
!            call cpchr2chr(record,MolBase(itmols)%name,lchr)
!            MolBase(itmols)%name(lchr+1) = C_NULL_CHAR
!            MolBase(itmols)%id = int(itmols,kind=C_INT)

ccc AB: next is for testing (comment out)
!            write(*,*)'FIELD scan: found molecule name - "'//
!     *                trim(mtmolname(itmols))
!!     *                //'" =?= "',(MolBase(itmols)%name(ich),ich=1,lchr),'"'
ccc
c AB: next is for testing (comment out)
!            write(*,*)
!            write(*,*)'"',(record(ich),ich=1,lchr),'"'
!            write(*,*)' =?= '
!            write(*,*)'"',(MolBase(itmols)%name(ich),ich=1,lchr),'"'
!            write(*,*)
ccc
            do while(loop2)
              
              call getrec(safe,idnode,nfield)
!              if(.not.safe)call abortscan(52,idnode)
              if(.not.safe) then 
                istate = -5
                return
              endif

              call lowcase(record,lenrec)
              
              ksite=0
              
              if(findstring('nummol',record,idum))then
                
                nummols=intstr(record,lenrec,idum)

                mtnmolecs(itmols) = nummols

!                MolBase(itmols)%nrept = int(nummols,kind=C_INT)

                nmols = nmols+nummols

ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found "nummols" directive - ',nummols
!     x               ,MolBase(itmols)%nrept,nmols
ccc
              elseif(findstring('atoms',record,idum))then
                
                numsit=intstr(record,lenrec,idum)
             
!                if(numsit.lt.1 .or. numsit.gt.mxatms) then
                if(numsit.lt.1 .or. numsit.gt.mxsite) then
!                  call abortscan(998,idnode)
!        if(.not.safe) then 
                  istate = -4
                  return
!        endif
                endif
   
                mtnmolatm(itmols)=numsit

!                MolBase(itmols)%nsites = int(numsit,kind=C_INT)

ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found "atoms" directive - ',numsit
!     *    ,MolBase(itmols)%nsites
!     *    ,MolBase(itmols)%nsites*MolBase(itmols)%nrept
ccc
!                mxatms=mxatms+numsit*nummols
!                mxsite=mxsite+numsit

                jsite=0
                if(itmols.gt.1) then
                  do jtmols=1,itmols-1
                    jsite=jsite+mtnmolatm(jtmols)
                  enddo
                endif

                ngrp=1
                nneu=1

                ksite=0
                do isite=1,numsit

                  if(ksite.lt.numsit)then
                    
                    jsite=jsite+1
                    
                    call getrec(safe,idnode,nfield)
!                    if(.not.safe)call abortscan(52,idnode)
                    if(.not.safe) then 
                       istate = -5
                       return
                    endif

                    call getword(name,record,8,lenrec)

                    if(len(trim(name)).lt.1) then
!                    if(trim(name).eq.'') then
!                       call abortscan(999,idnode)
!                       if(.not.safe) then 
                       istate = -5
                       return
!                       endif
                    endif

                    atmmass = dblstr(record,lenrec,idum)
                    
ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found atom mass = ',atmmass
ccc
                    if(atmmass.lt.0.d0) then
!                       call abortscan(997,idnode)
!                       if(.not.safe) then 
                       istate = -4
                       return
!                       endif
                    endif

                    atmcharge = dblstr(record,lenrec,idum)
                
ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found atom charge = ',atmcharge
ccc
                    nrept=intstr(record,lenrec,idum)

                    if(nrept.lt.1) nrept=1

                    numatms = numatms+nrept

ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: found atom nrept = ',nrept
ccc
                    if(numatms.gt.mxatms) then
!                       call abortscan(998,idnode)
!                       if(.not.safe) then 
                       istate = -4
                       return
!                       endif
                    endif

                    ifrz=intstr(record,lenrec,idum)
                    nneu=intstr(record,lenrec,idum)

                    if(lneut)then
                      if(nneu.ne.nold) nxn1=0
                      nxn1=nxn1+nrept
                      mxn1=max(mxn1,nxn1)
                      nold=nneu
c AB: number of groups within molecule
                      ngrp=max(ngrp,nneu)
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

                    lchr = lenchar(record,lenrec)

                    if( lchr.gt.0 ) then

                      ich =intstr(record,lenrec,idum)

                      call getword(type,record,8,lenrec)
                      
                      if(len(trim(type)).gt.0) then

                        call getword(type,record,8,lenrec)

                        call stripstr(type)
                        
                        if(len(trim(type)).lt.1) then
                          type = trim(name)
                        endif

                      else
                        type = trim(name)
                      endif

                    else
                      type = trim(name)
                    endif

c AB: we only go over blocks of molecule & atom types here, therefore
c AB: we need only one entry per repetitive site in each molecule type!

!                 do jsite=ksite+1,ksite+nrept
                    mtatmname(jsite) = trim(name)
                    mtatmtype(jsite) = trim(type)
                    mtatmmass(jsite) = atmmass
                    mtatmchrg(jsite) = atmcharge
                    mtatmnrep(jsite) = nrept
                    mtatmngrp(jsite) = nneu

!                    lchr = min(8,len(trim(name)))
!                    call cpstr2chr(name,FldSite(jsite)%name,lchr)
!                    FldSite(jsite)%name(lchr+1) = C_NULL_CHAR
!
!                    lchr = min(8,len(trim(type)))
!                    call cpstr2chr(type,FldSite(jsite)%type,lchr)
!                    FldSite(jsite)%type(lchr+1) = C_NULL_CHAR
!
!                    FldSite(jsite)%m = real(atmmass,kind=C_DOUBLE)
!                    FldSite(jsite)%q = real(atmcharge,kind=C_DOUBLE)
!                    FldSite(jsite)%idmol = int(isite,kind=C_INT)
!                    FldSite(jsite)%idgrp = int(nneu,kind=C_INT)
!                    FldSite(jsite)%ifrzn = int(ifrz,kind=C_INT)
!                    FldSite(jsite)%nrept = int(nrept,kind=C_INT)

!                 enddo

ccc AB: next is for testing (comment out)
!          write(*,'(3a,i4,a,i4,3a)')
!!     *    'FIELD scan: found atom named "',trim(lname1),'" #',
!     *    'FIELD scan: found atom named "',trim(mtatmname(jsite)),'" #',
!!     *    ksite+1,' in molecule type ',itmols,' "',
!     *    jsite,' in molecule type ',itmols
!     *    ,' "',trim(mtmolname(itmols)),'"'
ccc
ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: filled in atom "',
!     *          FldSite(jsite)%name,'" / "',FldSite(jsite)%type,
!     *          '" # ',isite,' - props for ',nrept,!atmmass,atmcharge,
!     *          ' identical atoms: ',
!     *          FldSite(jsite)%m,FldSite(jsite)%q
ccc
                    ksite=ksite+nrept

                  endif
                  
                enddo

                mtnmolgrp(itmols) = ngrp

!                MolBase(itmols)%ngroups = int(ngrp,kind=C_INT)

ccc AB: next is for testing (comment out)
!          write(*,*)'FIELD scan: filled in atom props for ',jsite,
!     *              ' unique atoms ...'
ccc
!                if(mmk.lt.mxatyp)   call abortscan(34,idnode)
                if(mmk.lt.mxatyp) then 
                   istate = -4
                   return
                endif

                if(lneut)mxneut=mxneut+nneu*nummols
                
              elseif(findstring('shell',record,idum))then
                
                numshl=intstr(record,40,idum)
                mxtshl=mxtshl+numshl
                mxshl=mxshl+nummols*numshl
                
                do ishls=1,numshl
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
   
                enddo
                
              elseif(findstring('bonds',record,idum))then
                
                numbonds=intstr(record,lenrec,idum)
                mxtbnd=mxtbnd+numbonds
                mxbond=mxbond+nummols*numbonds
                
                do ibonds=1,numbonds
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif

                enddo
                
              elseif(findstring('constr',record,idum))then
                
                numcon=intstr(record,lenrec,idum)
                mxtcon=mxtcon+numcon
                mxcons=mxcons+nummols*numcon
                
                do icon=1,numcon
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif

                enddo
                
              elseif(findstring('angles',record,idum))then
                
                numang=intstr(record,lenrec,idum)
                mxtang=mxtang+numang
                mxangl=mxangl+nummols*numang
                
                do iang=1,numang
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
                  
                enddo
                
              elseif(findstring('dihedr',record,idum))then
                
                numdih=intstr(record,lenrec,idum)
                mxtdih=mxtdih+numdih
                mxdihd=mxdihd+nummols*numdih
                
                do idih=1,numdih
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
                  
                enddo
                
              elseif(findstring('invers',record,idum))then
                
                numinv=intstr(record,lenrec,idum)
                mxtinv=mxtinv+numinv
                mxinv=mxinv+nummols*numinv
                
                do iinv=1,numinv
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
                  
                enddo
                
              elseif(findstring('rigid',record,idum))then
                
                numgrp=intstr(record,lenrec,idum)
                mxungp=mxungp+numgrp
                mxgrp=mxgrp+numgrp*nummols
                
                do kgrp=1,numgrp
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
                  
                  numgsit=intstr(record,lenrec,idum)
                  mxgatm=mxgatm+numgsit*nummols
                  mxngp=max(mxngp,numgsit)
                  do j=1,numgsit
                    
                    iii=intstr(record,lenrec,idum)
                    if(iii.eq.0)then
                      call getrec(safe,idnode,nfield)
!                      if(.not.safe)call abortscan(52,idnode)
                      if(.not.safe) then 
                         istate = -5
                         return
                      endif

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
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif
                  
                enddo
                
              elseif(findstring('pmf',record,idum))then
                
                do ipmf=1,2
                  
                  call getrec(safe,idnode,nfield)
!                  if(.not.safe)call abortscan(52,idnode)
                  if(.not.safe) then 
                     istate = -5
                     return
                  endif

                  call lowcase(record,lenrec)
                  npmf=intstr(record,lenrec,idum)
                  mxspmf=mxspmf+npmf
                  
                  do jpmf=1,npmf
                    
                    call getrec(safe,idnode,nfield)
!                    if(.not.safe)call abortscan(52,idnode)
                    if(.not.safe) then 
                       istate = -5
                       return
                    endif
                   
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
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            call lowcase(record,lenrec)
            if(findstring('tab',record,idum))ltable=.true.
            
          enddo
          mxvdw=max(mxvdw,(mxatyp*(mxatyp+1))/2)
          
          if(ltable)then
            goto 131

            if(idnode.eq.0)open(ntable,file='TABLE')
            
            call getrec(safe,idnode,ntable)
!            if(.not.safe)call abortscan(24,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            call getrec(safe,idnode,ntable)
!            if(.not.safe)call abortscan(24,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            ppp=dblstr(record,lenrec,idum)
            ppp=dblstr(record,lenrec,idum)
            mxgrid=max(mxgrid,intstr(record,lenrec,idum))
            
            close (ntable)

 131        continue
          endif
          
        elseif(findstring('metal',record,idum))then
          
          if(findstring('eam',record,idum))lmetab=.true.
          mxmet=intstr(record,lenrec,idum)
          do itpmet=1,mxmet
            
            call getrec(safe,idnode,nfield)
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            call lowcase(record,lenrec)
            if(findstring('eam',record,idum))lmetab=.true.
            
          enddo
          mxmet=max(mxmet,(mxatyp*(mxatyp+1))/2)
          
          if(lmetab)then
            
            if(idnode.eq.0)open(ntable,file='TABEAM')
            
            call getrec(safe,idnode,ntable)
!            if(.not.safe)call abortscan(24,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            call getrec(safe,idnode,ntable)
!            if(.not.safe)call abortscan(24,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            do i=1,intstr(record,lenrec,idum)
              
              call getrec(safe,idnode,ntable)
!              if(.not.safe)call abortscan(24,idnode)
              if(.not.safe) then 
                 istate = -5
                 return
              endif

              ngrid=intstr(record,lenrec,idum)
              mxgrid=max(mxgrid,ngrid+4)
              do j=1,(ngrid+3)/4
                
                call getrec(safe,idnode,ntable)
!                if(.not.safe)call abortscan(24,idnode)
                if(.not.safe) then 
                   istate = -5
                   return
                endif
                
              enddo
              
            enddo
            
            close (ntable)
            
          endif
          
        elseif(findstring('tbp',record,idum))then
          
          mxtbp=intstr(record,lenrec,idum)
          
          do itptbp=1,mxtbp
            
            call getrec(safe,idnode,nfield)
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

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
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

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
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            call getrec(safe,idnode,nfield)
!            if(.not.safe)call abortscan(52,idnode)
            if(.not.safe) then 
               istate = -5
               return
            endif

            rct=dblstr(record,lenrec,idum)
            rcter=max(rcter,rct)
            
          enddo
          
        elseif(findstring('extern',record,idum))then
          
          call getrec(safe,idnode,nfield)
!          if(.not.safe)call abortscan(52,idnode)
          if(.not.safe) then 
             istate = -5
             return
          endif

          nfld=intstr(record,lenrec,idum)
          if(nfld.eq.0)nfld=5
          call getrec(safe,idnode,nfield)
!          if(.not.safe)call abortscan(52,idnode)
          if(.not.safe) then 
             istate = -5
             return
          endif
          
          do k=1,nfld
            
            ppp=dblstr(record,lenrec,idum)
            if(idum.gt.lenrec.and.k.lt.nfld)then
              call getrec(safe,idnode,nfield)
!              if(.not.safe)call abortscan(52,idnode)
              if(.not.safe) then 
                 istate = -5
                 return
              endif
            endif
            
          enddo
          
        elseif(findstring('close',record,idum))then
          
          loop1=.false.
          
        endif
        
      enddo
      
      if(idnode.eq.0) close (nfield)
      
      if(mxpmf.gt.0)mxpmf=mxatms
      if(mxtcon.gt.0)mxexcl=max(mxexcl,6)
      if(mxtbnd.gt.0)mxexcl=max(mxexcl,6)
      if(mxtang.gt.0)mxexcl=max(mxexcl,16)
      if(mxtdih.gt.0)mxexcl=max(mxexcl,50)
      if(mxtinv.gt.0)mxexcl=max(mxexcl,50)
      if(mxneut.gt.0)mxexcl=max(mxexcl,10*mxn1*mxn1)
      if(mxgrp.gt.0)mxexcl=max(mxexcl,mxngp)
      
      return
      end subroutine fldread

      subroutine fldwrite
     x (idnode,fname,header,units,lneut,nmoltin,
     x  mtmolnamein,mtnmolecsin,mtnmolatmin,mtatmnamein,mtatmtypein,
     x  mtatmmassin,mtatmchrgin,mtatmnrepin,mxatmcg,istate)

c AB: NOTE - CG bead names & types should have been switched places (in the call)
c AB: so that beads of same type have same names in FIELD_CG
c AB: otherwise beads with different names but same type 
c AB: still add to the number of pair potentials (vdw)

      implicit none

      character(len=*) :: mtmolnamein(*)
      character(len=*) :: mtatmnamein(*),mtatmtypein(*)
      character(len=*) :: fname,header

      character(len=255) :: output
      character(len=80)  :: buffer,string
      character(len=8)   :: units,atname1,atname2

      character(len=8), allocatable :: atnamediff(:)

      real(8) :: mtatmmassin(*),mtatmchrgin(*)
      integer :: mtnmolecsin(*),mtnmolatmin(*),mtatmnrepin(*)

      integer :: idnode,nmoltin,mxatmcg,ntatmcg,npatmcg
      integer :: istate,ios,mst,lstr,im,ia,ib,ie,ii,is,it,ip,j
      logical :: lneut,isflt,in

      atname1 = ' '
      atname2 = ' '

      istate = 0
      ios    = 0

ccc AB: debugging only (comment out)
!      write(*,'(/,3a)')
!     x 'fldwrite(): DL_POLY Field (topology) file creation routine'
ccc
      open(nfield,file=trim(fname),status='new',iostat=ios)

      if(ios.ne.0) then

         write(*,'(/,3a)')
     x   'fldwrite(): WARNING: FIELD (topology) file "',trim(fname),
     x   '" is being replaced ...'
         
         open(nfield,file=trim(fname),status='replace',iostat=ios)

         if( ios.ne.0 ) then
           write(*,'(/,3a)')
     x     'fldwrite(): ERROR: FIELD (topology) file "',trim(fname),
     x     '" could not be open, escaping!'

           istate = -3

           return
         endif

      endif

ccc AB: debugging only (comment out)
!      write(*,'(/,3a)')'fldwrite(): starting Field file "',
!     x trim(fname),'"'
!
!      write(*,'(/,5a)')'fldwrite(): system title is "',trim(header),
!     x '" and units are "',trim(units),'"'
ccc

      lstr = len(trim(header))
      if( lstr.gt.80 ) then
        lstr = 80
        header = header(:lstr)
        write(*,'(/,4a)')'fldwrite(): WARNING - too long header,',
     x  ' truncated to "',trim(header),'" (80 symbols)'
      endif

      write(nfield,'(a)')trim(header)
      write(nfield,'(2a)')'UNITS   ',trim(units)

      if( lneut ) write(nfield,'(a)')'neutral groups'

      if( nmoltin.gt.0 ) then
        write(nfield,'(a,i4)')'Molecular types ',nmoltin
      else
        write(*,'(a,i3,a,/)')
     x   'fldwrite(): ERROR - no molecular types in the input (',
     x   nmoltin,' )!'
!     x   nmoltin,' ) - FULL STOP!'

        close(nfield)
        istate = -4

!        stop
        return
      endif

      if( ALLOCATED(atnamediff) ) DEALLOCATE (atnamediff)

      ALLOCATE (atnamediff(1:mxatms), STAT=mst )
      if( mst.ne.0 ) then

        output='fldwrite(): ERROR - not enough memory for char(8) '//
     x         'atnamediff'

        write(buffer,*)1
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)mxatms
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'==='

        istate = -6

!        stop
        return
      endif

      it = 0
      is = 0
      ip = 0
      DO im=1,nmoltin

        write(nfield,'(a)')mtmolnamein(im)

c AB: make the output look nice: numerics take only the needed space
        write(buffer,*)mtnmolecsin(im)
        ii = indexnum(buffer,isflt)

        write(nfield,'(2a)')'nummolecs  ',trim(buffer(ii:))
!        write(nfield,'(a,i8)')'nummolecs ',mtnmolecsin(im)

        write(buffer,*)mtnmolatmin(im)
        ii = indexnum(buffer,isflt)

        write(nfield,'(2a)')'atoms      ',trim(buffer(ii:))
!        write(nfield,'(a,i8)')'atoms ',mtnmolatmin(im)

        ip = ip+mtnmolatmin(im)*mtnmolecsin(im)

        Do ia=1,mtnmolatmin(im)

          atname2 = mtatmnamein(is+1)

          in = .false.
          do j=1,is
            in = mtatmnamein(j).eq.atname2
            if( in ) exit
          enddo

          if( .not.in ) then
            it = it+1
            atnamediff(it) = atname2
          endif

          is = is+1

          write(nfield,'(a8,2f12.4,4i5,2(a,a8),i8)')mtatmnamein(is),
     x    mtatmmassin(is),mtatmchrgin(is),mtatmnrepin(is),
     x    0,1,ia,' ',mtatmnamein(is),' ',mtatmtypein(is),is

        Enddo

!        write(nfield,'(a)')'shell # #'

!        write(nfield,'(a)')'rigid bodies #'

!        write(nfield,'(a)')'bonds #'

!        write(nfield,'(a)')'constraints #'

!        write(nfield,'(a)')'pmf Bond/A'

!        write(nfield,'(a)')'angles #'

!        write(nfield,'(a)')'dihedrals #'

!        write(nfield,'(a)')'inversions #'

!        write(nfield,'(a)')'teth #'

        write(nfield,'(a)')'finish'

      ENDDO

      ntatmcg = it
      npatmcg = it+it*(it-1)/2

      if(ip.ne.mxatmcg) then

        write(*,*)'---'
        write(*,*)
     x  'fldwrite(): ERROR - mismatch in number of sites : ',
     x  ip,' =?= ',mxatmcg,' expected!'
!     x  ,'- FULL STOP!'
        write(*,*)'==='

        istate = -4

!        stop
        return
      endif

      if(npatmcg.gt.0) then

c AB: make the output look nice: numerics take only the needed space
        write(buffer,*)npatmcg
        ii = indexnum(buffer,isflt)

        write(nfield,'(2a)')'vdw  ',trim(buffer(ii:))
!        write(nfield,'(a,i8)')'vdw ',npatmcg

        do ia=1,ntatmcg
          do ib=ia,ntatmcg

            write(nfield,'(2(a8,a),a8,2f12.4)')
     x      atnamediff(ia),' ',atnamediff(ib),' ','     tab',1.0,1.0
!     x      mtatmnamein(ia),' ',mtatmnamein(ib),' ','     tab',1.0,1.0

          enddo
        enddo

      endif

      write(nfield,'(a)')'close'

      close(nfield)

      if( ALLOCATED(atnamediff) ) DEALLOCATE (atnamediff)

      return
      end subroutine fldwrite


      subroutine cfgread
     x  (idnode,fname,header,imcon,levcfg,natoms,atname,
     x   coord,veloc,force,volm,xhi,yhi,zhi,cell,energy)
!     x   coord,veloc,force,volm,xhi,yhi,zhi,cell,buffer,energy)
      
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
c     the subroutine has been amended to read in the config/frame
c     by Andrey Brukhno June 2013
c     
c***********************************************************************
      
      character(len=*) :: fname,header

!      character(len=80) :: header
      character(len=*) :: atname(*)
      character(len=8) :: name

      logical :: lvolm

      real(8) :: coord(3,*)
      real(8) :: veloc(3,*)
      real(8) :: force(3,*)

      real(8) :: cell,celprp,buffer,extra,volm,xhi,yhi,zhi
      real(8) :: xxx,yyy,zzz,uuu,vvv,www,coz,energy

      integer :: idnode,ficfg,imcon,i,levcfg,natoms,line,line1
      dimension cell(9),celprp(10),buffer(10),extra(5)
      
      ficfg = 25

      xhi=0.d0
      yhi=0.d0
      zhi=0.d0
      volm=0.d0

      imcon=0
      lvolm=(imcon.eq.0.or.imcon.eq.6)

      do i=1,9
        cell(i)=0.d0
      enddo

      line  = 0
      line1 = 0

      if(idnode.eq.0) then
        
        open (ficfg,file=trim(fname),err=200)

        write(*,*)
        write(*,*)'cfgread(): reading CONFIG file "',trim(fname),'"'

c     read the CONFIG file header

        line = line+1
        read(ficfg,'(a)',end=200,err=200)header

        line = line+1
        read(ficfg,*,end=200,err=200)levcfg,imcon,natoms,energy

        if( levcfg.lt.0 .or. levcfg.gt.2 ) then
          write(*,*)'cfgread(): unrecognized CONFIG format: levcfg =',
     *    levcfg,' not in [0 1 2] - FULL STOP!'

          close(ficfg)
          stop
          return
!        elseif( natoms.lt.0 .or. natoms.gt.mxatms ) then
        elseif( abs(natoms).gt.mxatms ) then
          write(*,*)'cfgread(): too many atoms in CONFIG: natoms =',
     *    abs(natoms),' >',mxatms,
     *    ' (recompile for larger arrays) - FULL STOP!'

          close(ficfg)
          stop
          return
        endif

c     specify molecular dynamics simulation cell
        
        if(imcon.eq.-1)then

          line = line+1
          read(ficfg,*,end=200,err=200)cell(1)

          cell(5)=cell(1)
          cell(9)=cell(1)

          call dcell(cell,celprp)

          imcon=1
        
        elseif(imcon.eq.-2)then

          line = line+1
          read(ficfg,*,end=200,err=200)cell(1),cell(5),cell(9)

          call dcell(cell,celprp)

          imcon=2

        elseif(imcon.gt.0 .and. imcon.lt.8)then
          
          line = line+1
          read(ficfg,*,end=200,err=200)cell(1),cell(2),cell(3)

          line = line+1
          read(ficfg,*,end=200,err=200)cell(4),cell(5),cell(6)

          line = line+1
          read(ficfg,*,end=200,err=200)cell(7),cell(8),cell(9)

          call dcell(cell,celprp)
          
        else
          write(*,*)'cfgread(): unrecognized CONFIG format: imcon =',
     *    imcon,' not in [-2 -1 0 ... 7] - FULL STOP!'

          close(ficfg)
          stop
          return
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
        
!        natoms = min(natoms,mxatms)

        line1 = line

        i=0
        do 
          i=i+1
!        do i=1,natoms

c AB: for extracting coords only for |natoms| - the file may contain more!
          if(natoms.lt.0 .and. i.gt.-natoms) exit

c AB: for reading in coords for all atoms in the file unless arrays overflow
          if(i.gt.mxatms) exit

!          name = trim(atname(i))

          if(levcfg.eq.0) then
            
            line = line+1
            read(ficfg,'(a8)',end=100,err=200) atname(i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      coord(1,i),coord(2,i),coord(3,i)
            
          else if(levcfg.eq.1) then
            
            line = line+1
            read(ficfg,'(a8)',end=100,err=200) atname(i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      coord(1,i),coord(2,i),coord(3,i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      veloc(1,i),veloc(2,i),veloc(3,i)
         
          else
            
            line = line+1
            read(ficfg,'(a8)',end=100,err=200) atname(i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      coord(1,i),coord(2,i),coord(3,i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      veloc(1,i),veloc(2,i),veloc(3,i)
            
            line = line+1
            read(ficfg,*,end=100,err=200)
     x      force(1,i),force(2,i),force(3,i)

          endif
          
          if(lvolm)then
            
            if(i.eq.1)then
              
              xhi=abs(coord(1,i))
              yhi=abs(coord(2,i))
              zhi=abs(coord(3,i))
              
            else
              
              xhi=max(xhi,abs(coord(1,i)))
              yhi=max(yhi,abs(coord(2,i)))
              zhi=max(zhi,abs(coord(3,i)))
              
            endif

            if( abs(imcon).lt.3 .and. 
     *        ( xhi.gt.cell(1) .or. yhi.gt.cell(5) .or. yhi.gt.cell(9) )
     *        ) then

              write(*,'(2(a,3f10.5),a)')
     x        'cfgread(): check max dims =',xhi,yhi,zhi,
     x        ' vs cell dims =',cell(1),cell(5),cell(9),' - using PBC?'
              
            endif

          endif
          
        enddo
        
  100   continue
                
        close (ficfg)

        write(*,'(2(a,i6),a,/)')
     *  ' cfgread(): extracted coords for ',i-1,
     *  ' atoms (leftover: ',abs(natoms)-i+1,' atoms)'
!        write(*,'(a,i6,3a,i6,a,/)')
!     *  ' atoms from CONFIG file "',trim(fname),
!     *  '" (leftover: ',abs(natoms)-i+1,' atoms)'
        natoms = i-1

        if(imcon.eq.0)then
          
          volm=8.d0*xhi*yhi*zhi
          
        else if(imcon.eq.6)then
          
          coz=(cell(1)*cell(4)+cell(2)*cell(5)+cell(3)*cell(6))/
     x      (celprp(1)*celprp(2))
          volm=2.d0*zhi*celprp(1)*celprp(2)*sqrt(1.d0-coz**2)
          
        endif

      endif
      
      extra(1)=dble(imcon)
      extra(2)=xhi
      extra(3)=yhi
      extra(4)=zhi
      extra(5)=volm
!      call gdsum(extra,5,buffer)
!      call gdsum(cell,9,buffer)
      imcon=nint(extra(1))
      xhi=extra(2)
      yhi=extra(3)
      zhi=extra(4)
      volm=extra(5)
      
      return

 200  continue
      if( line.eq.0 ) then
        write(*,*)'cfgread(): could not open CONFIG file "',trim(fname),
     *  '" - FULL STOP!'
      else
        write(*,*)'cfgread(): could not read CONFIG file "',trim(fname),
     *  '" at line ',line,' - FULL STOP!'
      endif

      close(ficfg)

      stop
      return

      end subroutine cfgread

      subroutine cfgwrite
     x  (idnode,fname,header,imcon,levcfg,natoms,atname,
     x   coord,veloc,force,volm,xhi,yhi,zhi,cell,buffer)
      
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
c     the subroutine has been amended to writing out the config/frame
c     by Andrey Brukhno June 2013
c     
c***********************************************************************
      
      character(len=*) :: fname,header

!      character(len=80) :: header
      character(len=*) :: atname(*)
      character(len=8) :: name

      logical :: lvolm

      real(8) ::  coord(3,*)
      real(8) ::  veloc(3,*)
      real(8) ::  force(3,*)

      real(8) :: cell,celprp,buffer,extra,volm,xhi,yhi,zhi
      real(8) :: xxx,yyy,zzz,uuu,vvv,www,coz

      integer :: idnode,ficfg,imcon,i,levcfg,natoms
      dimension cell(9),celprp(10),buffer(10),extra(5)
      
      ficfg = 25

      xhi=0.d0
      yhi=0.d0
      zhi=0.d0
      volm=0.d0

!      imcon=0
      lvolm=(imcon.eq.0.or.imcon.eq.6)

!      do i=1,9
!        cell(i)=0.d0
!      enddo

      if(idnode.eq.0) then
        
        open (ficfg,file=trim(fname))
        
c     read the CONFIG file header
        
        write(ficfg,'(a)')trim(header)
        write(ficfg,'(3i10,1p,g20.12)')levcfg,imcon,natoms,0.0
        
c     specify molecular dynamics simulation cell
        
        if(imcon.gt.0)then
          
          write(ficfg,'(3f20.12)')cell(1),cell(2),cell(3)
          write(ficfg,'(3f20.12)')cell(4),cell(5),cell(6)
          write(ficfg,'(3f20.12)')cell(7),cell(8),cell(9)
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
        
        natoms = min(natoms,mxatms)

!        i=0
!        do while(.true.)
!          i=i+1
        do i=1,natoms

          name = trim(atname(i))

          if(levcfg.eq.0) then
            
            write(ficfg,'(a8,i10)') name,i
            write(ficfg,'(3g20.10)')
     x      coord(1,i),coord(2,i),coord(3,i)
            
          else if(levcfg.eq.1) then
            
            write(ficfg,'(a8,i10)') name,i
            write(ficfg,'(3g20.10)')
     x      coord(1,i),coord(2,i),coord(3,i)
            write(ficfg,'(3g20.12)')
     x      veloc(1,i),veloc(2,i),veloc(3,i)
         
          else
            
            write(ficfg,'(a8,i10)') name,i
            write(ficfg,'(3g20.10)')
     x      coord(1,i),coord(2,i),coord(3,i)
            write(ficfg,'(3g20.12)')
     x      veloc(1,i),veloc(2,i),veloc(3,i)
            write(ficfg,'(3g20.12)')
     x      force(1,i),force(2,i),force(3,i)

          endif
          
          if(lvolm)then
            
            if(i.eq.1)then
              
              xhi=abs(coord(1,i))
              yhi=abs(coord(2,i))
              zhi=abs(coord(3,i))
              
            else
              
              xhi=max(xhi,abs(coord(1,i)))
              yhi=max(yhi,abs(coord(2,i)))
              zhi=max(zhi,abs(coord(3,i)))
              
            endif

            write(*,'(a,3f10.5)')
     x      ' writecfg(): max dims =',xhi,yhi,zhi
            
          endif
          
        enddo
        
!  100   continue
        
        if(imcon.eq.0)then
          
          volm=8.d0*xhi*yhi*zhi
          
        else if(imcon.eq.6)then
          
          coz=(cell(1)*cell(4)+cell(2)*cell(5)+cell(3)*cell(6))/
     x      (celprp(1)*celprp(2))
          volm=2.d0*zhi*celprp(1)*celprp(2)*sqrt(1.d0-coz**2)
          
        endif
        
        close (ficfg)
        
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
      
      end subroutine cfgwrite

      subroutine config_write(fname,header,levcfg,imcon,natms,engcfg,
     x  cell,atname,coord,veloc,force,istate)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing CONFIG files
c     
c     copyright - daresbury laboratory 
c     author    - w. smith aug 2007
c
c     the subroutine has been amended to write out the config/frame
c     by Andrey Brukhno June 2013
c     
c***********************************************************************
      
      implicit none
      
      integer,parameter :: nconf=99

      character(len=*) :: fname,header
      character(len=*) :: atname(*)

      integer :: i,natms,levcfg,imcon,nstep,istate

      real(8) :: engcfg
      real(8) :: cell(9)

      real(8) :: coord(3,*)
      real(8) :: veloc(3,*)
      real(8) :: force(3,*)

      open(nconf,file=fname,form='formatted')
      
      write(nconf,'(a)')trim(header)
      write(nconf,'(3i10,1p,g20.12)') levcfg,imcon,natms,engcfg
      if(imcon.gt.0) write(nconf,'(3f20.12)') cell
      
      do i=1,natms
        
        write(nconf,'(a8,i10)') atname(i),i

        write(nconf,'(3g20.10)') 
     x    coord(1,i),coord(2,i),coord(3,i)

        if(levcfg.gt.0) write(nconf,'(3g20.12)')
     x    veloc(1,i),veloc(2,i),veloc(3,i)

        if(levcfg.gt.1) write(nconf,'(3g20.12)') 
     x    force(1,i),force(2,i),force(3,i)

      enddo
      
      close (nconf)
      
      return
      end subroutine config_write
      

      subroutine frmread
     x  (idnode,istrj,iflg,fname,title,imcon,keytrj,natms,nstep,tstep,
     x   energy,cell,name,weight,chge,xyz,vel,frc,istate)
!     x  (isnew,iflg,fname,title,name,imcon,keytrj,natms,
!     x   nstep,tstep,cell,chge,weight,xyz,vel,frc,istate)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading the formatted history file 
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c     amended  by Andrey Brukhno   May 2013
c
c     single processor version
c     
c***********************************************************************
      implicit none

      integer, parameter :: nhst=55,ncfg=77

      character(len=*) :: fname,title
      character(len=8) :: name(*),step

      real(8) :: cell(*),weight(*),chge(*)
      real(8) :: xyz(3,*),vel(3,*),frc(3,*)

      real(8) :: tstep,energy,vx,vy,vz,fx,fy,fz
      
      integer :: idnode,iflg,imcon,keytrj,natms,nstep,istate
      integer :: ipbc,matms,i,j

      integer,save :: fifrm,line,lfrm,ktrj,iopen=-1

      logical :: istrj
      logical :: isnew=.true.

c AB: isnew  = .true.  => the first instance of looking at HISTORY  (open it!)
c AB: iflg > 0  => reading HISTORY header only, to get the params (do NOT close it!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)

c AB: isnew  = .false. => further instances of looking at HISTORY (it IS open!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)
c AB: iflg = -1 => reading HISTORY done after next frame read (close it!)
c AB: iflg < -1 => reading HISTORY done, no frame to be read  (close it!)

c AB: istate = 0 => all right, everything goes according to plan
c AB: istate = 1 => all right, levtrj/keytrj reset to 0      (as found in HISTORY)
c AB: istate = 2 => all right, levtrj/keytrj reset to 0 or 1 (as found in HISTORY)
c AB: istate = 3 => all right, imcon reset to PBC/cell spec found in HISTORY
c AB: istate = 4 => all right, number of site entries reset to that found in HISTORY (<=mxatms)
c AB: istate =-1 => all right, end of HISTORY file reached with last frame read in full
c AB: istate =-2 => error: end of HISTORY file reached before a frame record finished
c AB: istate =-3 => error: HISTORY file not found
c AB: istate =-4 => error: number of site entries in HISTORY is greater than expected (>mxatms)

      istate = 0

      if( idnode.ne.0 ) return

      energy = 0.d0

c     open HISTORY/CONFIG file if new job

      line = 0
      lfrm = 0

      if( iflg.gt.0 .and. .not.isnew ) then

        if(istrj) then 
          close(nhst)
        else 
          close(ncfg)
        endif

        isnew = .true.

      endif

      if( isnew ) then

c AB: number of lines per frame header (adjustment)
        lfrm = lfrm+2
        
        if( istrj ) then
c AB: reading the header of trajectory file HISTORY

          fifrm  = nhst

          open(fifrm,file=trim(fname),status='old',err=100)

          read(fifrm,'(a80)',end=200,err=200) title
          line = line+1

          write(*,'(/,5a)')'frmread(): history file "',trim(fname),
     x    '" with header "',trim(title),'"'

!          read(fifrm,'(2i10)',end=200,err=200) ktrj,imcon
          read(fifrm,*,end=200,err=200) ktrj,ipbc,matms
          line = line+1

        else
c AB: reading the header of configuration file CONFIG

          fifrm  = ncfg

          open(fifrm,file=trim(fname),status='old',err=100)

          read(fifrm,'(a80)',end=200,err=200) title
          line = line+1

          write(*,'(/,5a)')'frmread(): config file "',trim(fname),
     x    '" with header "',trim(title),'"'

!          read(fifrm,'(2i10)',end=200,err=200) ktrj,imcon
          read(fifrm,*,end=200,err=200) ktrj,ipbc,matms,energy
          line = line+1

        endif

        isnew = .false.

        if(iflg.gt.0) then
c AB: iflg > 0  => reading the header only, to get the params (do NOT close it!)
c AB: keytrj,imcon,mxatms are taken from header
c AB: => memory (re)allocation should precede the actual reading of frames

          write(*,'(/,3a)')
     x       'frmread(): file "',trim(fname),
     x       '" has been read for parameters (only header was read)'

          keytrj = ktrj
          imcon  = ipbc
          natms  = matms

          iflg  = 0
!          isnew = .false.

          return
        endif

        if( keytrj.gt.ktrj ) then

          if( keytrj.eq.1 ) then

            write(*,*)'frmread(): WARNING - no velocities in file!'
!            write(*,'(a)')'frmread(): ERROR - no velocities in file'
            istate = 1
!            close (fifrm)
!            isnew=.true.
!            return

            write(*,*)'frmread(): resetting keytrj ',
     x      keytrj,' -> ',ktrj

          else

            write(*,*)'frmread(): WARNING - no vels./forces in file!'
!            write(*,'(a)')'frmread(): ERROR - no forces in file'
            istate = 2
!            close (fifrm)
!            isnew = .true.
!            return

            write(*,*)'frmread(): resetting keytrj ',
     x      keytrj,' -> ',ktrj

          endif

          keytrj = ktrj

        endif

        if( imcon.ne.ipbc ) then

          write(*,*)'frmread(): WARNING - cell PBC is not as expected!'
!          write(*,'(a)')'frmread(): ERROR - no forces in file'
          istate = 3
!          close (fifrm)
!          isnew = .true.
!          return

          write(*,*)'frmread(): resetting imcon ',
     x    imcon,' -> ',ipbc

        endif

        imcon = ipbc

      elseif(iflg.lt.-1) then
c AB: iflg < -1 => reading HISTORY/CONFIG done, frame is not read (close it!)

        close(fifrm)
        isnew = .true.
        iflg  = -3

        write(*,'(/,3a,/)')
     x       'frmread(): file "',trim(fname),
     x       '" has been closed (all reading done)'

        return
      endif
      
      if( istrj ) then
c AB: reading frame from HISTORY file (the parameters for the next frame)

        lfrm = lfrm+1

!        read(fifrm,'(a8,4i10,f12.6)',end=200,err=200)
        read(fifrm,*,end=200,err=200)
     x       step,nstep,matms,ktrj,imcon,tstep
        line = line+1

      else
c AB: reading frame from CONFIG file (no additional parameters after the header)
c AB: close the file after reading in a frame (only one frame is present)
        iflg  = -1
        matms = natms
      endif

      if( matms.gt.mxatms ) then
c AB: mxatms is not taken from here and is too small (= 0 ?)

        write(*,'(a)')'frmread(): ERROR - too many atoms to handle'
        write(*,'(3a,i8,a,i8,/)')'frmread(): file ',
     x   trim(fname),' contains',matms,' atoms > max =',mxatms

        close (fifrm)
        isnew  = .true.
        iflg   = -3

        istate = -4

        natms  = matms

        return
      elseif( matms.ne.natms ) then
c AB: the expected number of atoms is greater than the actually found
c AB: we need more or less automatic reading process, so correct whatever we can!

!        write(*,'(a)')'Error - too many atoms in MD cell'
       write(*,'(/,a)')'frmread(): WARNING - number of atoms was reset'
       write(*,'(3a,2(i8,a),/)')'frmread(): file "',
     x   trim(fname),'" contains',matms,
     x   ' atoms (=?= ',natms,' expected)'

        istate = 4
!        close (fifrm)
!        isnew=.true.
!        return

        natms = matms

      endif

      if(imcon.gt.0) then

        lfrm = lfrm+3

        read(fifrm,*,end=200,err=200) cell(1),cell(2),cell(3)
        line = line+1
        read(fifrm,*,end=200,err=200) cell(4),cell(5),cell(6)
        line = line+1
        read(fifrm,*,end=200,err=200) cell(7),cell(8),cell(9)
        line = line+1

      endif

      lfrm = lfrm+matms*(2+ktrj)

!      do i = 1,natms
      do i = 1,matms
      
       if( istrj ) then
c AB: reading frame from HISTORY file (atom name, index, mass and charge)

!        read(fifrm,'(a8,i10,2f12.6)',end=200,err=200)
        read(fifrm,*,end=200,err=200)
     x    name(i),j,weight(i),chge(i)
        line = line+1

       else
c AB: reading frame from CONFIG file (atom name and index only)

        read(fifrm,*,end=200,err=200)
     x    name(i),j
        line = line+1

       endif

        read(fifrm,*,end=200,err=200) xyz(1,i),xyz(2,i),xyz(3,i)
        line = line+1

        if(keytrj.gt.0)then
          read(fifrm,*,end=200,err=200) vel(1,i),vel(2,i),vel(3,i)
        line = line+1
        else if(ktrj.gt.0)then
          read(fifrm,*,end=200,err=200) vx,vy,vz
        line = line+1
        endif

        if(keytrj.gt.1)then
          read(fifrm,*,end=200,err=200) frc(1,i),frc(2,i),frc(3,i)
        line = line+1
        else if(ktrj.gt.1)then
          read(fifrm,*,end=200,err=200) fx,fy,fz
        line = line+1
        endif

      enddo

      if( iflg.lt.0 ) then
c AB: iflg = -1 => reading HISTORY/CONFIG done after the frame is read (close it!)

        close (fifrm)
        isnew = .true.
        iflg  = -3
!        iflg = 0

      endif

      return

  100 continue

      write(*,'(/,3a,/)')
     x 'frmread(): ERROR - file "',trim(fname),
     x '" not found, escaping...'
!     x '" not found - FULL STOP!'

      istate = -3
        iflg = -3
!      stop
      return

  200 continue

!      if( mod((line-2),(natms*(2+ktrj)+lfrm)).ne.0. ) then
      if( line.ne.0 .and. line.ne.lfrm ) then
        write(*,'(/,3a,2(i8,a),/)')
     x   'frmread(): ERROR - file "',trim(fname),
!     x   '" ended abnormally - FULL STOP! '
     x   '" ended abnormally after frame line '
     x   ,line,' (=?= ',lfrm,' expected), escaping...'

        close (fifrm)
        isnew = .true.
        iflg  = -3
!        iflg = 0

        istate = -2
!        stop
        return
      endif

      write(*,'(/,3a,2(i8,a))')
     x 'frmread(): WARNING - reached the end of file "',
     x trim(fname),'" at step ',nstep,', frame line ',line

      close (fifrm)
      isnew = .true.
      iflg  = -3
!      iflg = 0

      istate = -1

      return
      end subroutine frmread

      subroutine frmwrite
     x (idnode,istrj,iflg,fname,title,imcon,keytrj,natms,nstep,tstep,
     x  energy,cell,atname,charge,weight,coord,veloc,force,istate)

c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c
c***********************************************************************
      implicit none
      
      integer, parameter :: nhst=66, ncfg=88
      integer, parameter :: ndump=100

      character(len=*) :: fname,title
      character(len=8) :: atname(*)

      real(8) :: cell(9),weight(*),charge(*)
      real(8) :: coord(3,*),veloc(3,*),force(3,*)

      real(8) :: tstep,energy

      integer :: idnode,iflg,imcon,keytrj,natms,nstep,istate,i
      
      integer,save :: fofrm,iopen=-1,ios=0,idump=0
!      data    iopen/-1/
!      data    ios/0/

      logical :: istrj
      logical,save :: isnew=.true.
!      data    isnew/.true./

      istate = 0

      if( idnode.ne.0 ) return

!      if(istrj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if( isnew ) then
          
          idump = 0

          if( istrj ) then
c AB: writing trajectory file HISTORY
            fofrm  = nhst
          else
c AB: writing configuration file CONFIG
            fofrm  = ncfg
          endif

          if( iopen.lt.0 ) then

            open(fofrm,file=trim(fname),status='new',iostat=ios)

            if(ios.ne.0) then

              write(*,'(/,3a)')
     x        'frmwrite(): Warning - file "',trim(fname),
     x        '" is being replaced ...'

              open(fofrm,file=trim(fname),status='replace')

            endif
          
            if( istrj ) then
c AB: writing the header of trajectory file HISTORY

              write(fofrm,'(a)') title
              write(fofrm,'(3i10)') keytrj,imcon,natms

            else

c AB: writing the header of configuration file CONFIG/REVCON
              write(fofrm,'(a)') title
              write(fofrm,'(3i10,1p,g20.12)') keytrj,imcon,natms,energy

            endif

            if( iflg.gt.1 ) then
c AB: iflg > 0  => writing the header only (do NOT close it!)

              write(*,'(/,3a)')
     x         'frmwrite(): file "',trim(fname),
     x         '" has been started with parameters (only header)'

              isnew = .false.
              iopen = 0
              iflg  = 0

              return
            endif

          elseif( iopen.gt.0 ) then

            open(fofrm,file=trim(fname),status='old',position='append')

          endif

          isnew = .false.
          iopen = 0

        elseif(iflg.lt.-1) then

          close(fofrm)
          isnew = .true.
          iopen = -1
          iflg  = -3

          write(*,'(/,3a,/)')
     x     'frmwrite(): file "',trim(fname),
     x     '" has being closed (all writing done)'

          return
        endif
        

        if( istrj ) then
c AB: writing trajectory file HISTORY

          write(fofrm,'(a8,4i10,f12.6)') 'timestep',
     x     nstep,natms,keytrj,imcon,tstep

          if(imcon.gt.0) write(fofrm,'(3g12.4)') cell

          do i = 1,natms

            write(fofrm,'(a8,i10,2f12.6)')
     x        atname(i),i,weight(i),charge(i)

            write(fofrm,'(1p,3e12.4)')
     x        coord(1,i),coord(2,i),coord(3,i)

            if(keytrj.gt.0)then
              write(fofrm,'(1p,3e12.4)')
     x        veloc(1,i),veloc(2,i),veloc(3,i)
            endif

            if(keytrj.gt.1)then
              write(fofrm,'(1p,3e12.4)')
     x          force(1,i),force(2,i),force(3,i)
            endif
          
          enddo

        else
c AB: writing frame into CONFIG file (no additional parameters after the header)
c AB: close the file after writing down a frame (only one frame should be present)
        
          iflg = -1

          if(imcon.gt.0)then
          
            write(fofrm,'(3f20.12)')cell(1),cell(2),cell(3)
            write(fofrm,'(3f20.12)')cell(4),cell(5),cell(6)
            write(fofrm,'(3f20.12)')cell(7),cell(8),cell(9)
!          call dcell(cell,celprp)
          
          endif

          do i=1,natms

            write(fofrm,'(a8,i10)') atname(i),i

            write(fofrm,'(3g20.10)')
     x        coord(1,i),coord(2,i),coord(3,i)
            
            if(keytrj.gt.0) then
            
              write(fofrm,'(3g20.12)')
     x          veloc(1,i),veloc(2,i),veloc(3,i)
         
            endif

            if(keytrj.gt.1) then
            
              write(fofrm,'(3g20.12)')
     x          force(1,i),force(2,i),force(3,i)

            endif
          
          enddo

        endif
        
        idump = idump+1
        
!        if( .not.isnew .and. mod(nstep,ndump).eq.0 ) then
        if( .not.isnew ) then

          if(iflg.lt.0) then
c AB: iflg = -1 => writing HISTORY/CONFIG done (close it!)

            close (fofrm)
            isnew = .true.
            iopen = -1
            iflg  = -3
!            iflg = 0

          elseif( mod(idump,ndump).eq.0 ) then
          
c     close history file at regular intervals

            close (fofrm)
            isnew = .true.
            iopen = 1
!            iflg  = -3

          endif

        endif
      
      return
      end subroutine frmwrite


      subroutine trjread
     x  (new,iflg,fname,title,imcon,keytrj,natms,nstep,tstep,
     x   cell,name,weight,chge,xyz,vel,frc,istate)
!     x  (new,iflg,fname,title,name,imcon,keytrj,natms,
!     x   nstep,tstep,cell,chge,weight,xyz,vel,frc,istate)
      
c***********************************************************************
c     
c     dl_poly subroutine for reading the formatted history file 
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c     amended  by Andrey Brukhno   May 2013
c
c     single processor version
c     
c***********************************************************************
      implicit none

      integer, parameter :: nhist=77

      character(len=*) :: fname,title
      character(len=8) :: name(*),step

      real(8) :: cell(*),chge(*),weight(*)
      real(8) :: xyz(3,*),vel(3,*),frc(3,*)

      real(8) :: tstep,vx,vy,vz,fx,fy,fz
      
      integer :: iflg,imcon,keytrj,natms,nstep,istate
      integer :: ipbc,matms,i,j

      integer,save :: line,lfrm,ktrj

      logical :: new

c AB: new  = .true.  => the first instance of looking at HISTORY  (open it!)
c AB: iflg > 0  => reading HISTORY header only, to get the params (do NOT close it!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)

c AB: new  = .false. => further instances of looking at HISTORY (it IS open!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)
c AB: iflg = -1 => reading HISTORY done after next frame read (close it!)
c AB: iflg < -1 => reading HISTORY done, no frame to be read  (close it!)

c AB: istate = 0 => all right, everything goes according to plan
c AB: istate = 1 => all right, levtrj/keytrj reset to 0      (as found in HISTORY)
c AB: istate = 2 => all right, levtrj/keytrj reset to 0 or 1 (as found in HISTORY)
c AB: istate = 3 => all right, imcon reset to PBC/cell spec found in HISTORY
c AB: istate = 4 => all right, number of site entries reset to that found in HISTORY (<=mxatms)
c AB: istate =-1 => all right, end of HISTORY file reached with last frame read in full
c AB: istate =-2 => error: end of HISTORY file reached before a frame record finished
c AB: istate =-3 => error: HISTORY file not found
c AB: istate =-4 => error: number of site entries in HISTORY is greater than expected (>mxatms)

      istate = 0

c     open HISTORY file if new job

      if( new ) then

c AB: number of lines per frame header (adjustment)
        lfrm = 1
        line = 0
        
        open(nhist,file=trim(fname),status='old',err=100)
        
        read(nhist,'(a80)',end=200,err=200) title
        line = line+1

        write(*,'(/,5a)')'trjread(): History file "',trim(fname),
     x   '" with header "',trim(title),'"'

!        read(nhist,'(2i10)',end=200,err=200) ktrj,imcon
        read(nhist,*,end=200,err=200) ktrj,ipbc,matms
        line = line+1

        new  = .false.

        if(ipbc.gt.0) lfrm = 4

        if(iflg.gt.0) then
c AB: iflg > 0  => reading HISTORY header only, to get the params (do NOT close it!)
c AB: keytrj,imcon,mxatms are taken from header
c AB: => memory (re)allocation should precede the actual reading of frames

          write(*,'(/,3a)')
     x       'trjread(): trajectory file "',trim(fname),
     x       '" has been read for parameters (only header was read)'

          keytrj = ktrj
          imcon  = ipbc
          natms  = matms

          iflg = 0

          return
        endif

        if( keytrj.gt.ktrj ) then

          if( keytrj.eq.1 ) then

            write(*,*)'trjread(): WARNING - no velocities in file!'
!            write(*,'(a)')'trjread(): ERROR - no velocities in file'
            istate = 1
!            close (nhist)
!            new=.true.
!            return

            write(*,*)'trjread(): resetting keytrj ',
     x      keytrj,' -> ',ktrj

          else

            write(*,*)'trjread(): WARNING - no vels./forces in file!'
!            write(*,'(a)')'trjread(): ERROR - no forces in file'
            istate = 2
!            close (nhist)
!            new=.true.
!            return

            write(*,*)'trjread(): resetting keytrj ',
     x      keytrj,' -> ',ktrj

          endif

          keytrj = ktrj

        endif

        if( imcon.ne.ipbc ) then

          write(*,*)'trjread(): WARNING - unit cell is not as expected!'
!          write(*,'(a)')'trjread(): ERROR - no forces in file'
          istate = 3
!          close (nhist)
!          new=.true.
!          return

          write(*,*)'trjread(): resetting imcon ',
     x    imcon,' -> ',ipbc

        endif

        imcon = ipbc

      elseif(iflg.lt.-1) then
c AB: iflg < -1 => reading HISTORY done, no frame to be read (close it!)

        close(nhist)
        new  = .true.
        iflg = 0

        write(*,'(/,3a,/)')
     x       'trjread(): trajectory file "',trim(fname),
     x       '" has been closed (all reading done correctly)'

        return
      endif
      
!      read(nhist,'(a8,4i10,f12.6)',end=200,err=200)
      read(nhist,*,end=200,err=200)
     x     step,nstep,matms,ktrj,imcon,tstep
      line = line+1
      
      if( matms.gt.mxatms ) then
c AB: mxatms is not taken from here and is too small (= 0 ?)

        write(*,'(a)')'trjread(): ERROR - too many atoms to handle'
        write(*,'(3a,i8,a,i8,/)')'trjread(): History file ',
     x   trim(fname),' contains',matms,' atoms > max =',mxatms

        istate = -4
        close (nhist)
        new=.true.
        
        natms = matms

        return
      elseif( matms.ne.natms ) then
c AB: the expected number of atoms is greater than the actually found
c AB: we need more or less automatic reading process, so correct whatever we can!

!        write(*,'(a)')'Error - too many atoms in MD cell'
       write(*,'(/,a)')'trjread(): WARNING - number of atoms was reset'
       write(*,'(3a,2(i8,a),/)')'trjread(): History file "',
     x   trim(fname),'" contains',matms,
     x   ' atoms (=?= ',natms,' expected)'

        istate = 4
!        close (nhist)
!        new=.true.
!        return

        natms = matms

      endif

      if(imcon.gt.0) then

        read(nhist,*,end=200,err=200) cell(1),cell(2),cell(3)
        line = line+1
        read(nhist,*,end=200,err=200) cell(4),cell(5),cell(6)
        line = line+1
        read(nhist,*,end=200,err=200) cell(7),cell(8),cell(9)
        line = line+1

      endif

!      do i = 1,natms
      do i = 1,matms

!        read(nhist,'(a8,i10,2f12.6)',end=200,err=200)
        read(nhist,*,end=200,err=200)
     x    name(i),j,weight(i),chge(i)
        line = line+1

        read(nhist,*,end=200,err=200) xyz(1,i),xyz(2,i),xyz(3,i)
        line = line+1

        if(keytrj.ge.1)then
          read(nhist,*,end=200,err=200) vel(1,i),vel(2,i),vel(3,i)
        line = line+1
        else if(ktrj.ge.1)then
          read(nhist,*,end=200,err=200) vx,vy,vz
        line = line+1
        endif

        if(keytrj.ge.2)then
          read(nhist,*,end=200,err=200) frc(1,i),frc(2,i),frc(3,i)
        line = line+1
        else if(ktrj.ge.2)then
          read(nhist,*,end=200,err=200) fx,fy,fz
        line = line+1
        endif

      enddo

      if( iflg.lt.0 ) then
c AB: iflg = -1 => reading HISTORY done after current frame read (close it!)

        close (nhist)
        new  = .true.
        iflg = 0

      endif

      return

  100 continue

      write(*,'(/,3a,/)')
     x 'trjread(): ERROR - History file "',trim(fname),
     x '" not found - FULL STOP!'

      istate = -3
      stop
      return

  200 continue

      if( mod((line-2),(natms*(2+ktrj)+lfrm)).ne.0. ) then
        write(*,'(/,3a,3i7,/)')
     x   'trjread(): ERROR - History file "',trim(fname),
     x   '" ended abnormally - FULL STOP! '
     x   ,line,lfrm,(natms*(2+ktrj)+lfrm)

        close (nhist)
        new  = .true.
        iflg = 0

        istate = -2
        stop
        return
      endif

      write(*,'(/,3a,i8)')
     x 'trjread(): WARNING - reached the end of History file "',
     x trim(fname),'" at step ',nstep

      close (nhist)
      new  = .true.
      iflg = 0

      istate = -1

      return
      end subroutine trjread

      subroutine trjwrite
     x (idnode,ltraj,fname,title,imcon,istraj,keytrj,natms,cell,
     x  atname,charge,weight,coord,veloc,force,ndump,nstraj,nstep,tstep)

c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c
c***********************************************************************

!      use setup_module
!      use config_module

      implicit none
      
      character(len=*) :: fname,title
      character(len=8) :: atname(*)

      logical :: ltraj
      logical,save :: newjob
!=.true.
      data    newjob/.true./

      integer,save :: iopen,ios
!=-1,ios=0
      data    iopen/-1/
      data    ios/0/

      integer :: idnode,imcon,istraj,keytrj,natms,nstraj,nstep,ndump,i

      real(8) :: tstep
      real(8) :: cell(9)

      real(8) :: weight(*)
      real(8) :: charge(*)

      real(8) :: coord(3,*)
      real(8) :: veloc(3,*)
      real(8) :: force(3,*)
      
      if(ltraj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if(newjob)then
          
          newjob = .false.

          if( iopen.lt.0 ) then

            open(nhist,file=trim(fname),status='new',iostat=ios)

            if(ios.ne.0) then

              write(*,'(/,3a)')
     x        'trjwrite(): Warning - History file "',trim(fname),
     x        '" is being replaced ...'

              open(nhist,file=trim(fname),status='replace')

            endif
          
            write(nhist,'(a)') title
            write(nhist,'(3i10)') keytrj,imcon,natms

          elseif( iopen.gt.0 ) then

            open(nhist,file=trim(fname),status='old',position='append')

          endif

          iopen = 0

        elseif(iopen.gt.1) then

          close(nhist)
          newjob = .true.

          write(*,'(/,3a)')
     x         'trjwrite(): WARNING - History file "',trim(fname),
     x         '" is being closed, escaping ...'
          return
        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then
          
          write(nhist,'(a8,4i10,f12.6)') 'timestep',
     x         nstep,natms,keytrj,imcon,tstep

          if(imcon.gt.0) write(nhist,'(3g12.4)') cell

          do i = 1,natms

            write(nhist,'(a8,i10,2f12.6)')
     x        atname(i),i,weight(i),charge(i)

            write(nhist,'(1p,3e12.4)')
     x        coord(1,i),coord(2,i),coord(3,i)

            if(keytrj.ge.1)then
              write(nhist,'(1p,3e12.4)')
     x          veloc(1,i),veloc(2,i),veloc(3,i)
            endif

            if(keytrj.ge.2)then
              write(nhist,'(1p,3e12.4)')
     x          force(1,i),force(2,i),force(3,i)
            endif

          enddo

        endif

c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0) then
          
          close (nhist)
          newjob = .true.
          iopen  = 1

        endif
        
      endif
      
      return
      end subroutine trjwrite


      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c***********************************************************************

      implicit none

      real(8) :: a,b,d,r
      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert


      subroutine pbc_xyz(dxyz,cell,rcell)

      implicit none
      real(8) :: dxyz(3),cell(9),rcell(9)
      real(8) :: dxr,dyr,dzr
      
      dxr = dxyz(1)*rcell(1) + dxyz(2)*rcell(4) + dxyz(3)*rcell(7)
      dyr = dxyz(1)*rcell(2) + dxyz(2)*rcell(5) + dxyz(3)*rcell(8)
      dzr = dxyz(1)*rcell(3) + dxyz(2)*rcell(6) + dxyz(3)*rcell(9)
      
      dxr = dxr-nint(dxr)
      dyr = dyr-nint(dyr)
      dzr = dzr-nint(dzr)
      
      dxyz(1) = dxr*cell(1) + dyr*cell(4) + dzr*cell(7)
      dxyz(2) = dxr*cell(2) + dyr*cell(5) + dzr*cell(8)
      dxyz(3) = dxr*cell(3) + dyr*cell(6) + dzr*cell(9)
      
      return
      end subroutine pbc_xyz


c***********************************************************************
c     
c     dl_poly routines for allocating and deallocating arrays
c     used by the CG mapping routines above
c
c     copyright - daresbury laboratory 2013
c     author Andrey Brukhno       June 2013
c     
c***********************************************************************

      logical function memall(reset)

      implicit none
      logical :: reset

      imemall = 0
      imemfld = 0
      imemfrm = 0

      memall = .true.

      if( .not.memstr('mtmolname',mtmolname,1,mxtmls,reset) )then
         memall = .false.
      elseif(.not.memstr('mtatmname',mtatmname,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memstr('mtatmname',mtatmtype,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memstr('fratmname',fratmname,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memstr('fratmtype',fratmtype,1,mxatms,reset))then
         memall = .false.
!      elseif(.not.memstr('',,1,,reset))then
!         memall = .false.
      elseif(.not.memdbl2('xyz',xyz,1,3,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl2('vel',vel,1,3,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl2('frc',frc,1,3,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl('mtatmmass',mtatmmass,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl('mtatmchrg',mtatmchrg,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl('fratmmass',fratmmass,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memdbl('fratmchrg',fratmchrg,1,mxatms,reset))then
         memall = .false.
!      elseif(.not.memdbl('',,1,,reset))then
!         memall = .false.
      elseif(.not.memint('mtnmolecs',mtnmolecs,1,mxtmls,reset))then
         memall = .false.
      elseif(.not.memint('mtnmolatm',mtnmolatm,1,mxtmls,reset))then
         memall = .false.
      elseif(.not.memint('mtnmolatm',mtnmolgrp,1,mxtmls,reset))then
         memall = .false.
      elseif(.not.memint('mtatmnrep',mtatmnrep,1,mxatms,reset))then
         memall = .false.
      elseif(.not.memint('mtatmnrep',mtatmngrp,1,mxatms,reset))then
         memall = .false.
!      elseif(.not.memint('',,1,,reset))then
!         memall = .false.
!      elseif()then
      endif

      if( reset ) then
        imemall = 2
        imemfld = 1
        imemfrm = 1
      endif

      return
      end function memall

      logical function memfld(reset)

      implicit none
      logical :: reset

      imemfld = 0

      memfld = .true.

      if( .not.memstr('mtmolname',mtmolname,1,mxtmls,reset) )then
         memfld = .false.
      elseif(.not.memstr('mtatmname',mtatmname,1,mxatms,reset))then
         memfld = .false.
      elseif(.not.memstr('mtatmname',mtatmtype,1,mxatms,reset))then
         memfld = .false.
!      elseif(.not.memstr('',,1,,reset))then
!         memfld = .false.
      elseif(.not.memdbl('mtatmmass',mtatmmass,1,mxatms,reset))then
         memfld = .false.
      elseif(.not.memdbl('mtatmchrg',mtatmchrg,1,mxatms,reset))then
         memfld = .false.
!      elseif(.not.memdbl('',,1,,reset))then
!         memfld = .false.
      elseif(.not.memint('mtnmolecs',mtnmolecs,1,mxtmls,reset))then
         memfld = .false.
      elseif(.not.memint('mtnmolatm',mtnmolatm,1,mxtmls,reset))then
         memfld = .false.
      elseif(.not.memint('mtnmolatm',mtnmolgrp,1,mxtmls,reset))then
         memfld = .false.
      elseif(.not.memint('mtatmnrep',mtatmnrep,1,mxatms,reset))then
         memfld = .false.
      elseif(.not.memint('mtatmnrep',mtatmngrp,1,mxatms,reset))then
         memfld = .false.
!      elseif(.not.memint('',,1,,reset))then
!         memfld = .false.
!      elseif()then
      endif

      if( reset ) imemfld = 1

      return
      end function memfld

      logical function memfrm(reset)

      implicit none
      logical :: reset

      imemfrm = 0

      memfrm = .true.

      if(.not.memstr('fratmname',fratmname,1,mxatms,reset))then
         memfrm = .false.
      elseif(.not.memstr('fratmtype',fratmtype,1,mxatms,reset))then
         memfrm = .false.
!      elseif(.not.memstr('',,1,,reset))then
!         memfrm = .false.
      elseif(.not.memdbl2('xyz',xyz,1,3,1,mxatms,reset))then
         memfrm = .false.
      elseif(.not.memdbl2('vel',vel,1,3,1,mxatms,reset))then
         memfrm = .false.
      elseif(.not.memdbl2('frc',frc,1,3,1,mxatms,reset))then
         memfrm = .false.
      elseif(.not.memdbl('fratmmass',fratmmass,1,mxatms,reset))then
         memfrm = .false.
      elseif(.not.memdbl('fratmchrg',fratmchrg,1,mxatms,reset))then
         memfrm = .false.
!      elseif(.not.memdbl('',,1,,reset))then
!         memfrm = .false.
!      elseif(.not.memint('',,1,,reset))then
!         memfrm = .false.
!      elseif()then
      endif

      if( reset ) imemfrm = 1

      return
      end function memfrm


      logical function incmxtmls(mxtmol0,mxtmol1)

      implicit none
      integer :: mxtmol0,mxtmol1

      incmxtmls = .true.

      if( .not.incstr('mtmolname',mtmolname,1,mxtmol0,mxtmol1) )then
         incmxtmls = .false.
!      elseif(.not.incstr('',,1,mxtmol0,mxtmol1))then
!         incmxtmls = .false.
      elseif(.not.incint('mtnmolecs',mtnmolecs,1,mxtmol0,mxtmol1))then
         incmxtmls = .false.
      elseif(.not.incint('mtnmolatm',mtnmolatm,1,mxtmol0,mxtmol1))then
         incmxtmls = .false.
      elseif(.not.incint('mtnmolatm',mtnmolgrp,1,mxtmol0,mxtmol1))then
         incmxtmls = .false.
!      elseif(.not.incint('',,1,mxtmol0,mxtmol1))then
!         incmxtmls = .false.
      endif

      return
      end function incmxtmls


      logical function incmxatms(mxatm0,mxatm1)

      implicit none
      integer :: mxatm0,mxatm1

      incmxatms = .true.

      if( .not.incstr('mtatmname',mtatmname,1,mxatm0,mxatm1) )then
         incmxatms = .false.
      elseif(.not.incstr('mtatmtype',mtatmtype,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incstr('fratmname',fratmname,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incstr('fratmtype',fratmtype,1,mxatm0,mxatm1))then
         incmxatms = .false.
!      elseif(.not.incstr('',,1,,reset))then
!         incmxatms = .false.
      elseif(.not.incdbl2('xyz',xyz,1,3,3,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl2('vel',vel,1,3,3,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl2('frc',frc,1,3,3,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl('mtatmmass',mtatmmass,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl('mtatmchrg',mtatmchrg,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl('fratmmass',fratmmass,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incdbl('fratmchrg',fratmchrg,1,mxatm0,mxatm1))then
         incmxatms = .false.
!      elseif(.not.incdbl2('',,1,3,3,1,mxatm0,mxatm1))then
!         incmxatms = .false.
      elseif(.not.incint('mtatmnrep',mtatmnrep,1,mxatm0,mxatm1))then
         incmxatms = .false.
      elseif(.not.incint('mtatmnrep',mtatmngrp,1,mxatm0,mxatm1))then
         incmxatms = .false.
!      elseif(.not.incint('',,1,,reset))then
!         incmxatms = .false.
      endif

      return
      end function incmxatms


      logical function memstr(str_name,str_array,min0,max0,reset)

      implicit none

c AB: name of the array
      character(len=*)   :: str_name
      character(len=255) :: output
      character(len=80)  :: buffer

      integer :: min0,max0,memstate,ib
      logical :: reset,isflt

c AB: the actual array to allocate
      character(len=*), allocatable :: str_array(:)

      memstr   = .true.
      memstate = 0

      if( ALLOCATED(str_array) ) DEALLOCATE(str_array)
      
      if( .not.reset ) return

      ALLOCATE( str_array(min0:max0), STAT=memstate )

      if( memstate.ne.0 ) then

        output='memstr(): Error - not enough memory for char*(*) '//
     x         trim(str_name)

        write(buffer,*)min0
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)max0
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'---'

        memstr  = .false.

      endif

      return
      end function memstr
      
      
      logical function incstr(str_name,str_array,min0,max0,max1)

      implicit none

c AB: name of the array
      character(len=*)    :: str_name
      character(len=srec) :: str_nameN

      integer :: min0,max0,max1,maxt,i

c AB: the actual array to allocate
      character(len=*),    allocatable :: str_array(:)
      character(len=lrec), allocatable :: str_arrayN(:)

      incstr   = .true.

      str_nameN = trim(str_name)//'[inc]'

      if( memstr(str_nameN,str_arrayN,min0,max1,.true.) ) then
c AB: success allocating memory for temporary array
c AB: store the data from the original array before re-allocating it

        maxt = min(max0,max1)
        do i=min0,maxt
          str_arrayN(i)=trim(str_array(i))
        enddo

        if( memstr(str_name,str_array,min0,max1,.true.) ) then
c AB: success allocating memory for increased array
c AB: restore the data from the temporary array into the original one

          do i=min0,max1
            if(i.le.max0) then
              str_array(i)=trim(str_arrayN(i))
            else
              str_array(i)=''
            endif
          enddo

          if( ALLOCATED(str_arrayN) ) DEALLOCATE(str_arrayN)

        else
c AB: failed to allocate memory for increased array
         incstr = .false.
        endif

      else
c AB: failed to allocate memory for temporary larger array
         incstr = .false.
      endif

      return
      end function incstr


      logical function memint(str_name,int_array,min0,max0,reset)

      implicit none

c AB: name of the array
      character(len=*)   :: str_name
      character(len=255) :: output
      character(len=80)  :: buffer

      integer :: min0,max0,memstate,ib
      logical :: reset,isflt

c AB: the actual array to allocate
      integer, allocatable :: int_array(:)

      memint   = .true.
      memstate = 0

      if( ALLOCATED(int_array) ) DEALLOCATE(int_array)
      
      if( .not.reset ) return

      ALLOCATE( int_array(min0:max0), STAT=memstate )

      if( memstate.ne.0 ) then

        output='memstr(): Error - not enough memory for integer '//
     x         trim(str_name)

        write(buffer,*)min0
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)max0
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'---'

        memint  = .false.

      endif

      return
      end function memint
      

      logical function incint(str_name,int_array,min0,max0,max1)

      implicit none

c AB: name of the array
      character(len=*)     :: str_name
      character(len=srec)  :: str_nameN

      integer :: min0,max0,max1,maxt,i

c AB: the actual array to allocate
      integer, allocatable :: int_array(:)
      integer, allocatable :: int_arrayN(:)

      incint   = .true.

      str_nameN = trim(str_name)//'[inc]'

      if( memint(str_nameN,int_arrayN,min0,max1,.true.) ) then
c AB: success allocating memory for temporary array
c AB: store the data from the original array before re-allocating it

        maxt = min(max0,max1)
        do i=min0,maxt
          int_arrayN(i)=int_array(i)
        enddo

        if( memint(str_name,int_array,min0,max1,.true.) ) then
c AB: success allocating memory for increased array
c AB: restore the data from the temporary array into the original one

          do i=min0,max1
            if(i.le.max0) then
              int_array(i)=int_arrayN(i)
            else
              int_array(i)=0
            endif
          enddo

          if( ALLOCATED(int_arrayN) ) DEALLOCATE(int_arrayN)

        else
c AB: failed to allocate memory for increased array
         incint = .false.
        endif

      else
c AB: failed to allocate memory for temporary larger array
         incint = .false.
      endif

      return
      end function incint


      logical function memflt(str_name,flt_array,min0,max0,reset)

      implicit none

c AB: name of the array
      character(len=*)   :: str_name
      character(len=255) :: output
      character(len=80)  :: buffer

      integer :: min0,max0,memstate,ib
      logical :: reset,isflt

c AB: the actual array to allocate
      real(4), allocatable :: flt_array(:)

      memflt   = .true.
      memstate = 0

      if( ALLOCATED(flt_array) ) DEALLOCATE(flt_array)
      
      if( .not.reset ) return

      ALLOCATE( flt_array(min0:max0), STAT=memstate )

      if( memstate.ne.0 ) then

        output='memstr(): Error - not enough memory for real(4) '//
     x         trim(str_name)

        write(buffer,*)min0
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)max0
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'---'

        memflt  = .false.

      endif

      return
      end function memflt
      

      logical function incflt(str_name,flt_array,min0,max0,max1)

      implicit none

c AB: name of the array
      character(len=*)    :: str_name
      character(len=srec) :: str_nameN

      integer :: min0,max0,max1,maxt,ib,i

c AB: the actual array to allocate
      real(4), allocatable :: flt_array(:)
      real(4), allocatable :: flt_arrayN(:)

      incflt   = .true.

      str_nameN = trim(str_name)//'[inc]'

      if( memflt(str_nameN,flt_arrayN,min0,max1,.true.) ) then
c AB: success allocating memory for temporary array
c AB: store the data from the original array before re-allocating it

        maxt = min(max0,max1)
        do i=min0,maxt
          flt_arrayN(i)=flt_array(i)
        enddo

        if( memflt(str_name,flt_array,min0,max1,.true.) ) then
c AB: success allocating memory for increased array
c AB: restore the data from the temporary array into the original one

          do i=min0,max1
            if(i.gt.max0) then
              flt_array(i) = 0.0
            else
              flt_array(i) = flt_arrayN(i)
            endif
          enddo

          if( ALLOCATED(flt_arrayN) ) DEALLOCATE(flt_arrayN)

        else
c AB: failed to allocate memory for increased array
         incflt = .false.
        endif

      else
c AB: failed to allocate memory for temporary larger array
         incflt = .false.
      endif

      return
      end function incflt


      logical function memdbl(str_name,dbl_array,min0,max0,reset)

      implicit none

c AB: name of the array
      character(len=*)   :: str_name
      character(len=255) :: output
      character(len=80)  :: buffer

      integer :: min0,max0,memstate,ib
      logical :: reset,isflt

c AB: the actual array to allocate
      real(8), allocatable :: dbl_array(:)

      memdbl   = .true.
      memstate = 0

      if( ALLOCATED(dbl_array) ) DEALLOCATE(dbl_array)
      
      if( .not.reset ) return

      ALLOCATE( dbl_array(min0:max0), STAT=memstate )

      if( memstate.ne.0 ) then

        output='memstr(): Error - not enough memory for real(8) '//
     x         trim(str_name)

        write(buffer,*)min0
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)max0
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'---'

        memdbl  = .false.

      endif

      return
      end function memdbl
      

      logical function incdbl(str_name,dbl_array,min0,max0,max1)

      implicit none

c AB: name of the array
      character(len=*)    :: str_name
      character(len=srec) :: str_nameN

      integer :: min0,max0,max1,maxt,ib,i

c AB: the actual array to allocate
      real(8), allocatable :: dbl_array(:)
      real(8), allocatable :: dbl_arrayN(:)

      incdbl = .true.

      str_nameN = trim(str_name)//'[inc]'

      if( memdbl(str_nameN,dbl_arrayN,min0,max1,.true.) ) then
c AB: success allocating memory for temporary array
c AB: store the data from the original array before re-allocating it

        maxt = min(max0,max1)
        do i=min0,maxt
          dbl_arrayN(i)=dbl_array(i)
        enddo

        if( memdbl(str_name,dbl_array,min0,max1,.true.) ) then
c AB: success allocating memory for increased array
c AB: restore the data from the temporary array into the original one

          do i=min0,max1
            if(i.gt.max0) then
              dbl_array(i) = 0.d0
            else
              dbl_array(i) = dbl_arrayN(i)
            endif
          enddo

          if( ALLOCATED(dbl_arrayN) ) DEALLOCATE(dbl_arrayN)

        else
c AB: failed to allocate memory for increased array
         incdbl = .false.
        endif

      else
c AB: failed to allocate memory for temporary larger array
         incdbl = .false.
      endif

      return
      end function incdbl


      logical function memdbl2(str_name,dbl_array2,
     x                         min1,max1,min2,max2,reset)

      implicit none

c AB: name of the array
      character(len=*)   :: str_name

      character(len=255) :: output
      character(len=80)  :: buffer

      integer :: min1,max1,min2,max2,memstate,ib
      logical :: reset,isflt

c AB: the actual array to allocate
      real(8), allocatable :: dbl_array2(:,:)

      memdbl2  = .true.
      memstate = 0

      if( ALLOCATED(dbl_array2) ) DEALLOCATE(dbl_array2)
      
      if( .not.reset ) return

      ALLOCATE( dbl_array2(min1:max1,min2:max2), STAT=memstate )

      if( memstate.ne.0 ) then

        output='memstr(): Error - not enough memory for real(8) '//
     x         trim(str_name)

        write(buffer,*)min1
        ib = indexnum(buffer,isflt)

        output = trim(output)//'('//trim(buffer(ib:))

        write(buffer,*)max1
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))

        write(buffer,*)min2
        ib = indexnum(buffer,isflt)

        output = trim(output)//','//trim(buffer(ib:))

        write(buffer,*)max2
        ib = indexnum(buffer,isflt)

        output = trim(output)//':'//trim(buffer(ib:))//')'

        write(*,*)'---'
        write(*,*) trim(output)
        write(*,*)'---'

        memdbl2  = .false.

      endif

      return
      end function memdbl2
      

      logical function incdbl2(str_name,dbl_array2,
     x                         min1,max10,max11,min2,max20,max21)

      implicit none

c AB: name of the array
      character(len=*)     :: str_name
      character(len=srec)  :: str_nameN

      integer :: min1,max10,max11,min2,max20,max21,maxt1,maxt2,i,j

c AB: the actual array to allocate
      real(8), allocatable :: dbl_array2(:,:)
      real(8), allocatable :: dbl_array2N(:,:)

      incdbl2 = .true.

      str_nameN = trim(str_name)//'[inc]'

      if( memdbl2(str_nameN,dbl_array2N,min1,max11,min2,max21,.true.) ) 
     x  then
c AB: success allocating memory for temporary array
c AB: store the data from the original array before re-allocating it

        maxt1 = min(max10,max11)
        maxt2 = min(max20,max21)

        do j=min2,maxt2
          do i=min1,maxt1
            dbl_array2N(i,j) = dbl_array2(i,j)
          enddo
        enddo

        if( memdbl2(str_name,dbl_array2,min1,max11,min2,max21,.true.) ) 
     x    then
c AB: success allocating memory for increased array
c AB: restore the data from the temporary array into the original one

          do j=min2,max21
            if(j.gt.max20) then
              do i=min1,max11
                dbl_array2(i,j) = 0.d0
              enddo
            else
             do i=min1,max11
              if(i.gt.max10) then
                dbl_array2(i,j) = 0.d0
              else
                dbl_array2(i,j) = dbl_array2N(i,j)
              endif
             enddo
            endif
          enddo

          if( ALLOCATED(dbl_array2N) ) DEALLOCATE(dbl_array2N)

        else
c AB: failed to allocate memory for increased array
         incdbl2 = .false.
        endif

      else
c AB: failed to allocate memory for temporary larger array
         incdbl2 = .false.
      endif

      return
      end function incdbl2


      subroutine inichar(chars,imax)
      implicit none
      character(len=1) :: chars(*)
      integer :: i,imax
      
      if(imax.gt.lenrec) imax=lenrec

      do i=1,imax
         chars(i) = ' '
      enddo
      
      return
      end subroutine inichar

      integer function lenchar(chars,imax)

c***********************************************************************
c     
c     DL_POLY routine to find the length of a char*1 array up to blanks
c
c     copyright daresbury laboratory 2013
c     author Andrey Brukhno      May 2013
c     
c***********************************************************************
      implicit none

      character(len=1) :: chars(*)
      integer :: i,imax
      
      if(imax.gt.lenrec) imax=lenrec

      lenchar=0

      do i=imax,1,-1
        if( chars(i).ne.' ' ) then 
          lenchar=i
          exit
!        elseif( chars(i).ne.' ' ) then 
!          lenchar=i
        endif
      enddo
     
      return
      end function lenchar

      integer function lenchar2(chars,chsep,imax)

c***********************************************************************
c     
c     DL_POLY routine to find the length of a char*1 array up to blanks
c
c     copyright daresbury laboratory 2013
c     author Andrey Brukhno      May 2013
c     
c***********************************************************************
      implicit none

      character(len=1) :: chars(*),chsep
      integer :: i,imax
      
      if(imax.gt.lenrec) imax=lenrec

      lenchar2=imax

      do i=1,imax
        if( chars(i).eq.chsep ) then 
          lenchar2=i-1
          exit
        endif
      enddo
      
      return
      end function lenchar2

      subroutine stripstr(string)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a char*(n) array string
c     maximum length is lenrec characters (check above)
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c
c     amended for character strings of arbitrary length
c     by Andrey Brukhno          May 2013
c     
c***********************************************************************
      implicit none

      character(len=*) :: string
      integer :: i,imax,j

      imax = len(trim(string))

      do i=1,imax

        if( string(1:1).eq.' ' .and. len(trim(string)).gt.1 ) then 
          string = string(2:)
        else
          exit
        endif

!        if(string(1:1).eq.' ')then
!          do j=1,imax-1
!            string(j:j)=string(j+1:j+1)
!          enddo
!          string(imax:imax)=' '
!        endif

      enddo

      return
      end subroutine stripstr

      subroutine stripQ(string)

c***********************************************************************
c     
c     DL_POLY routine to strip outmost quotation signs around a char*(*) string
c     works for character strings of arbitrary length
c     
c     copyright daresbury laboratory 2013
c     author - Andrey Brukhno    May 2013
c     
c***********************************************************************

      implicit none

      character(len=*) :: string
      character(len=1) :: char1
      integer :: imax

c AB: strip the leading blanks, while skipping the trailing ones
      call stripstr(string)

      char1 = string(1:1)

      if( char1.eq.'"' .or. char1.eq."'" ) then

c AB: find the rightmost closing Q-sign
        imax = index(trim(string),char1,.true.)

        if( imax.gt.1 ) then
c AB: remove only a Q-pair

          string = string(2:imax-1)

c AB: strip the leading blanks if any remain
          call stripstr(string)

        endif

      endif

      return
      end subroutine stripQ

      subroutine stripQleft(string)

c***********************************************************************
c     
c     DL_POLY routine to strip leftmost quotation signs and take string in between
c     works for character strings of arbitrary length
c     
c     copyright daresbury laboratory 2013
c     author - Andrey Brukhno    May 2013
c     
c***********************************************************************

      implicit none

      character(len=*) :: string
      integer :: ii,jj,ib,ie

      ii = index(string,'"')
      jj = index(string,"'")

      if( ii.eq.0 ) then
        ib = jj
      elseif( ii.lt.jj .or. jj.eq.0 ) then
        ib = ii
      else
        ib = jj
      endif

      if( ib.gt.0 ) then
c AB: remove quotation signs only if a pair found

        ie = index(string(ib+1:),string(ib:ib))

        if( ie.gt.0 ) string = string(ib+1:ib+ie-1)

      endif

      return
      end subroutine stripQleft

      subroutine cpchr2chr(chars1,chars2,imax)
      implicit none
      character(len=1) :: chars1(*),chars2(*)
      integer :: imax,l1,l2,i,j
      
      if(imax.gt.lenrec) imax=lenrec

      l1 = lenchar(chars1,imax)
      l2 = lenchar(chars2,imax)

!      if(l1.gt.l2) l1=l2
      do i=1,l1
         chars2(i) = chars1(i)
      enddo

      do j=i,l2
         chars2(j) = ' '
      enddo
      
      return
      end subroutine cpchr2chr

      subroutine cpstr2chr(chars1,chars2,imax)
      implicit none
      character(len=*) :: chars1
      character(len=1) :: chars2(*)
      integer :: imax,l1,l2,i,j
      
      if(imax.gt.lenrec) imax=lenrec

      l1 = len(trim(chars1))
      l2 = lenchar(chars2,imax)

      if(l1.gt.imax) l1=imax
      do i=1,l1
         chars2(i) = chars1(i:i)
      enddo

      do j=i,l2
         chars2(j) = ' '
      enddo
      
      return
      end subroutine cpstr2chr

      subroutine cpchr2str(chars1,chars2,imax)
      implicit none
      character(len=1) :: chars1(*)
      character(len=*) :: chars2
      integer :: imax,l1,l2,i,j
      
      if(imax.gt.lenrec) imax=lenrec

      l1 = lenchar(chars1,imax)
      l2 = len(trim(chars2))

      do i=1,l1
         chars2(i:i) = chars1(i)
      enddo

!      if(l2.gt.imax) l2=imax
      do j=i,l2
         chars2(j:j) = ' '
      enddo
      
      return
      end subroutine cpchr2str

      character(len=lenrec) function chr2str(chars,imax)

c***********************************************************************
c     
c     DL_POLY routine to convert a char*1 array into a character string
c
c     copyright daresbury laboratory 2013
c     author Andrey Brukhno      May 2013
c     
c***********************************************************************
      implicit none

      character(len=1) :: chars(*)
      integer :: i,j,imax,lchr,lstr
      
      if(imax.gt.lenrec) imax=lenrec

      lchr = lenchar(chars,imax)
!      lstr = min(imax,lchr)
!      lstr = min(lstr,lenrec)

!      do i=1,lstr
      do i=1,lchr
        chr2str(i:i)=chars(i)
      enddo
  
      do j=i,lenrec
        chr2str(j:j)=' '
      enddo
    
      return
      end function chr2str

      subroutine lowcstr(string)

c***********************************************************************
c     
c     DL_POLY routine to lowercase a char*1 array string of up to 255 characters.
c     Transportable to non-ASCII machines (?)
c     
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c
c     amended for character strings of arbitrary length
c     by Andrey Brukhno          May 2013
c     
c***********************************************************************

      implicit none

      character(len=*) :: string
      character(len=1) :: letter
      integer :: i,length

      length = len(trim(string))

!      do i=1,min(255,length)
      do i=1,length

        letter=string(i:i)

        if(letter.eq.'A')then
          letter='a'
        else if(letter.eq.'B')then
          letter='b'
        else if(letter.eq.'C')then
          letter='c'
        else if(letter.eq.'D')then
          letter='d'
        else if(letter.eq.'E')then
          letter='e'
        else if(letter.eq.'F')then
          letter='f'
        else if(letter.eq.'G')then
          letter='g'
        else if(letter.eq.'H')then
          letter='h'
        else if(letter.eq.'I')then
          letter='i'
        else if(letter.eq.'J')then
          letter='j'
        else if(letter.eq.'K')then
          letter='k'
        else if(letter.eq.'L')then
          letter='l'
        else if(letter.eq.'M')then
          letter='m'
        else if(letter.eq.'N')then
          letter='n'
        else if(letter.eq.'O')then
          letter='o'
        else if(letter.eq.'P')then
          letter='p'
        else if(letter.eq.'Q')then
          letter='q'
        else if(letter.eq.'R')then
          letter='r'
        else if(letter.eq.'S')then
          letter='s'
        else if(letter.eq.'T')then
          letter='t'
        else if(letter.eq.'U')then
          letter='u'
        else if(letter.eq.'V')then
          letter='v'
        else if(letter.eq.'W')then
          letter='w'
        else if(letter.eq.'X')then
          letter='x'
        else if(letter.eq.'Y')then
          letter='y'
        else if(letter.eq.'Z')then
          letter='z'
        endif

        string(i:i)=letter

      enddo

      return
      end subroutine lowcstr

      integer function indexnum(string,isflt)

c***********************************************************************
c     
c     DL_POLY routine to find the start of numerical data within string
c
c     Find the first digit in the string, if any
c     If found, analyse the characters around it to figure out 
c     if the figure complies with the floating point format.
c     Include the minus sign, if any
c
c     copyright daresbury laboratory 2013
c     author Andrey Brukhno      May 2013
c     
c***********************************************************************

      implicit none

      character(len=*) :: string
      character(len=1) :: char1
      integer :: i,lstr,ib,ie
      logical :: isflt

      isflt = .false.

      indexnum=0

      lstr=len(trim(string))
      do i=1,lstr

        char1 = string(i:i)

        if( char1.eq.'0' ) then
          indexnum = i
          exit
        elseif( char1.eq.'1' ) then
          indexnum = i
          exit
        elseif( char1.eq.'2' ) then
          indexnum = i
          exit
        elseif( char1.eq.'3' ) then
          indexnum = i
          exit
        elseif( char1.eq.'4' ) then
          indexnum = i
          exit
        elseif( char1.eq.'5' ) then
          indexnum = i
          exit
        elseif( char1.eq.'6' ) then
          indexnum = i
          exit
        elseif( char1.eq.'7' ) then
          indexnum = i
          exit
        elseif( char1.eq.'8' ) then
          indexnum = i
          exit
        elseif( char1.eq.'9' ) then
          indexnum = i
          exit
        endif

      enddo

      if( indexnum.gt.0 ) then

        ib = indexnum+1

        i = indexnum-1
        if(i.gt.0) then

          char1 = string(i:i)

c AB: check for the floating point before the figure
          if( char1.eq.'.' ) then 

            isflt = .true.
            indexnum = i
            i = i-1
            if(i.lt.1) return
            
            char1 = string(i:i)
            
          endif

c AB: check for the 'minus' sign before the figure ('plus' doesn't matter)
          if( char1.eq.'-' ) then 

            indexnum = i
            if( isflt ) return

          endif

        endif

c AB: the beginning has been found, but floating point remains unclear

        lstr = index(string(ib:),' ')-1
        if( lstr.lt.1 ) return
c AB: blank right after the first digit found, so it's done

        lstr = ib+lstr-1
        ie   = lstr-1

c AB: check for the floating point in the figure (before the next blank space)
        if( index(string(ib:lstr),'.').gt.0 ) then
          isflt = .true.
          return
        elseif(ie.gt.ib-1) then

          if( index(string(ib:ie),'e').gt.0 ) then
c AB: check for the exponent in the figure (third before the last non-blank/space)
            isflt = .true.
            return
          elseif( index(string(ib:ie),'E').gt.0 ) then
            isflt = .true.
            return
          elseif( index(string(ib:ie),'d').gt.0 ) then
            isflt = .true.
            return
          elseif( index(string(ib:ie),'D').gt.0 ) then
            isflt = .true.
            return
          endif

        endif

      endif

      return
      end function indexnum


      end module dlp_io_module

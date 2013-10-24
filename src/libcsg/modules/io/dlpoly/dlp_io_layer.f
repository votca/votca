      
c***********************************************************************
c     
c     DL_POLY wrapper-routines needed for minimalistic I/O for 
c     CONFIG, FIELD & HISTORY files (wrappers for VOTCA calls) 
c     NOTE: maximum length of strings is defined by lenrec=255 (parse_module.f)
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************

      subroutine field_scan(iflgC,mxatmsC,mtatmsC,mxtmlsC)

!      use iso_c_binding
      use, intrinsic :: iso_c_binding

!      use setup_module
      use dlp_io_module

      implicit none

      integer(kind=C_INT),intent(inout) :: iflgC,mxatmsC,mtatmsC,mxtmlsC

      integer :: idnode,mxn1
      real(8) :: rctbp,rcfbp,rcter

      idnode = 0
      mxn1   = 0
      rctbp  = 0.0
      rcfbp  = 0.0
      rcter  = 0.0

c AB: scan FIELD file for the max numbers, before allocating the memory

      call fldscan(idnode,mxn1,rctbp,rcfbp,rcter)

      mxatmsC = int(mxatms,kind=C_INT)
      mtatmsC = int(mxsite,kind=C_INT)
      mxtmlsC = int(mxtmls,kind=C_INT)

      imemall = 0
      imemfld = 0
      imemfrm = 0

      if( iflgC.eq.2_C_INT .and. .not.memall(.true.) ) then
        write(*,*)'---'
        write(*,*)
     x       'field_scan(2): Failed allocating memory!'
        write(*,*)'==='

        iflgC = -6_C_INT
!        stop
!        return
      elseif( iflgC.eq.1_C_INT .and. .not.memfld(.true.) ) then
        write(*,*)'---'
        write(*,*)
     x       'field_scan(1): Failed allocating memory!'
        write(*,*)'==='

        iflgC = -6_C_INT
!        stop
!        return
      else
        write(*,*)
        write(*,*)
     x  'field_scan(1/2): (re)allocating memory - SUCCESS!'

        iflgC = 0_C_INT
      endif

      imemall = imemfld+imemfrm
!      iflgC   =-int(imemall,kind=C_INT)
!      write(*, *)
!      write(*, *)'field_scan(): ISO C bound'
!      write(*, *)'field_scan(): max nums = ',mxatms,mxsite,mxtmls
!!      write(*, *)'field_scan(): max nums = ',mxatmsC,mtatmsC,mxtmlsC

      return
      end subroutine field_scan


      subroutine field_read(istatC,FieldBaseC,MolecBaseC,FieldSiteC)
     x bind (C, name="field_read_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for reading FIELD file into C-struct
c     equivalent derived types and passing to VOTCA DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: istatC,iflgC

      TYPE(C_PTR),target, intent(inout) :: FieldBaseC
      TYPE(C_PTR),target, intent(inout) :: MolecBaseC
      TYPE(C_PTR),target, intent(inout) :: FieldSiteC

      TYPE(FieldSpecsT), pointer :: FldBase
      TYPE(MolecSpecsT), pointer :: MolBase(:)
      TYPE(FieldSiteT),  pointer :: FldSite(:)

!      character(kind=C_CHAR),intent(in) :: fnameC(81)

      character(len=srec) :: fname,title
      character(len=half) :: units

      integer :: natms,nmols,nmolt,matms,mtatms
      integer :: idnode,iflg,istate,lchr,ich,im,ia,ima

      logical, save :: is_mem_new=.false.,lneut=.false.
!      logical, save :: is_mem_new,lneut
!      data is_mem_new/.false./,lneut/.false./

      fname  = 'FIELD'
      title  = ' '
      units  = ' '

      idnode = 0
      istate = 0
      lchr   = 0

      iflgC = istatC
      iflg  = int(iflgC)

c AB: scan FIELD file for the max numbers, before allocating the memory
!      call fldscan(idnode,mxn1,rctbp,rcfbp,rcter)

      CALL C_F_POINTER(C_LOC(FieldBaseC), FldBase)

c AB: allocate all the needed big arrays

      is_mem_new = (iflg.gt.0)

      if( is_mem_new ) then
        nmolt  = int(FldBase%nmols)
        natms  = int(FldBase%natms)
        mxtmls = nmolt
        mxsite = natms
      else
!        nmoltC = int(mxtmls,kind=C_INT)
!        natmsC = int(mxsite,kind=C_INT)
        nmolt  = mxtmls
        natms  = mxsite
      endif

      is_mem_new = is_mem_new .or. (imemfld.lt.1)

      if( is_mem_new ) then
        if( .not.memfld(.true.) ) then
          write(*,*)'---'
          write(*,*)
     x       'field_read(1): Failed allocating memory!'
          write(*,*)'==='

          istatC = -6_C_INT
!          stop
          return
        else
          write(*,*)
          write(*,*)
     x       'field_read(1): (re)allocating memory - SUCCESS!'
        endif
      endif

      is_mem_new = .false.
      imemall = imemfld+imemfrm

      CALL C_F_POINTER(C_LOC(MolecBaseC), MolBase, shape=[nmolt])
      CALL C_F_POINTER(C_LOC(FieldSiteC), FldSite, shape=[natms])

      matms = mxatms

!      lchr = lenchar2(fnameC,C_NULL_CHAR,81)
!      lchr = min(lchr,lenchar(fnameC,80))
!      call cpchr2str(fnameC,fname,lchr)

ccc AB: next is for testing (comment out)
!      write(*,*)
!      write(*,*)'field_read(.): ISO C bound'
!      write(*,*)'field_read(.): max nums = ',matms,natms,nmolt
!!      write(*,*)'field_read(): ',matms,natms,nmolt,' "',trim(fname),'"'
!
!      lchr = lenchar2(FldBase%title,C_NULL_CHAR,81)
!      lchr = min(lchr,lenchar(FldBase%title,80))
!      call cpchr2str(FldBase%title,title,lchr)
!      write(*,*)
!      write(*,*)'field_read(.): check C-header "',trim(title),'"',lchr
!!     x     (FldBase%title(ich),ich=1,lchr),'"',lchr
!      write(*,*)
ccc
      call IniFldSpecs(FldBase)
      call IniMolSpecs(MolBase,nmolt)
      call IniFldSites(FldSite,natms)

c AB: read the necessary bits from FIELD file

      call fldread(idnode,fname,title,units,lneut,
     x             nmolt,nmols,natms,istate)

      istatC = int(istate,kind=C_INT)

      lchr = min(80,len_trim(title))
      call cpstr2chr(title,FldBase%title,lchr)
      FldBase%title(lchr+1) = C_NULL_CHAR

      lchr = min(8,len_trim(units))
      call cpstr2chr(units,FldBase%units,lchr)
      FldBase%units(lchr+1) = C_NULL_CHAR

      FldBase%ineut = 0_C_INT
      if( lneut ) FldBase%ineut = 1_C_INT
      
      FldBase%nmols = int(nmolt,kind=C_INT)
      FldBase%natms = int(natms,kind=C_INT)

      matms = 0
      ia = 0

      do im=1,mxtmls

         lchr = len_trim(mtmolname(im))
         call cpstr2chr(trim(mtmolname(im)),MolBase(im)%name,lchr)
         MolBase(im)%name(lchr+1) = C_NULL_CHAR

         MolBase(im)%id      = int(im,kind=C_INT)
         MolBase(im)%nrept   = int(mtnmolecs(im),kind=C_INT)
         MolBase(im)%nsites  = int(mtnmolatm(im),kind=C_INT)
         MolBase(im)%ngroups = int(mtnmolgrp(im),kind=C_INT)

ccc AB: debugging only (comment out)
!         write(*,*)'field_read(): check molecule type "',
!     x     (MolBase(im)%name(ich),ich=1,lchr),'" =?= "',
!     x     trim(mtmolname(im)),'", ',MolBase(im)%id,
!     x     MolBase(im)%nrept,MolBase(im)%nsites,MolBase(im)%ngroups
!!     x     mtnmolecs(im),mtnmolatm(im),mtnmolgrp(im)
ccc
         mtatms = MolBase(im)%nsites

         do ima=1,mtatms
            ia = ia+1

            lchr = len_trim(mtatmname(ia))
            call cpstr2chr(trim(mtatmname(ia)),FldSite(ia)%name,lchr)
            FldSite(ia)%name(lchr+1) = C_NULL_CHAR

            lchr = len_trim(mtatmtype(ia))
            call cpstr2chr(trim(mtatmtype(ia)),FldSite(ia)%type,lchr)
            FldSite(ia)%type(lchr+1) = C_NULL_CHAR

            FldSite(ia)%m = real(mtatmmass(ia),kind=C_DOUBLE)
            FldSite(ia)%q = real(mtatmchrg(ia),kind=C_DOUBLE)

            FldSite(ia)%nrept = int(mtatmnrep(ia),kind=C_INT)
            FldSite(ia)%idgrp = int(mtatmngrp(ia),kind=C_INT)

            matms = matms+mtnmolecs(im)*mtatmnrep(ia)

         end do

      end do

      if( ia.ne.natms ) then
        write(*,*)
        write(*,*)'field_read(-1): check atoms # in moltypes - ',
     x            ia,' =?= ',natms

        istatC = -1_C_INT
        iflg = -1
      endif

      if( matms.ne.mxatms ) then
        write(*,*)
        write(*,*)'field_read(-2): check atoms # in total - ',
     x            matms,' =?= ',mxatms

        istatC = -2_C_INT
        iflg = -2
      endif

c AB: deallocate all the big arrays (see above)
      if( iflg.lt.0 .or. istate.lt.0 ) then

        if( iflg.eq.-1 .and. .not.memfld(.false.) ) then

          write(*,*)'---'
          write(*,*)
     x    'field_read(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!          stop
!          return
        elseif( iflg.eq.-2 .and. .not.memall(.false.) ) then

          write(*,*)'---'
          write(*,*)
     x    'field_read(-2): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -8_C_INT
!          stop
!          return
        else
          write(*,*)
          write(*,*)
     x    'field_read(-1/-2): deallocating memory - SUCCESS!'
        endif

        imemall = imemfld+imemfrm

      endif

      return
      end subroutine field_read


      subroutine field_write(istatC,
     x FieldBaseC,MolecBaseC,FieldSiteC)
     x bind (C, name="field_write_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for writing FIELD file from C-struct
c     equivalent derived types passed from VOTCA DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
      use iso_c_binding
!      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: istatC,iflgC

      TYPE(C_PTR),target, intent(inout) :: FieldBaseC
      TYPE(C_PTR),target, intent(inout) :: MolecBaseC
      TYPE(C_PTR),target, intent(inout) :: FieldSiteC

      TYPE(FieldSpecsT), pointer :: FldBase
      TYPE(MolecSpecsT), pointer :: MolBase(:)
      TYPE(FieldSiteT),  pointer :: FldSite(:)

      character(len=srec) :: fname,header
      character(len=half) :: units

      integer :: natms,nmols,nmolt,matms,mtatms
      integer :: idnode,iflg,istate,lchr,ich,im,ia,ima

      logical, save :: is_mem_new=.false.,lneut=.false.

      fname  = 'FIELD_CG'
      header = ' '
      units  = ' '

      idnode = 0
      istate = 0
      lchr   = 0

      iflgC = istatC
      iflg  = int(iflgC)

      CALL C_F_POINTER(C_LOC(FieldBaseC), FldBase)

      nmolt  = int(FldBase%nmols)
      natms  = int(FldBase%natms)

c AB: allocate all the needed big arrays (see above)

      is_mem_new = (iflg.gt.0)

      if( is_mem_new .or. natms.gt.mxsite .or. nmolt.gt.mxtmls ) then
!         mxatms = natms
         mxsite = natms
         mxtmls = nmolt
      endif

      is_mem_new = is_mem_new .or. (imemfld.lt.1)

      if( is_mem_new ) then
        if ( .not.memfld(.true.) ) then
          write(*,*)'---'
          write(*,*)
     x       'field_write(1): Failed allocating memory!'
          write(*,*)'==='

          istatC = -6_C_INT
!          stop
          return
        else
          write(*,*)
          write(*,*)
     x       'field_write(1): (re)allocating memory - SUCCESS!'
        endif
      endif

      is_mem_new = .false.
      imemall = imemfld+imemfrm

      CALL C_F_POINTER(C_LOC(MolecBaseC), MolBase, shape=[nmolt])
      CALL C_F_POINTER(C_LOC(FieldSiteC), FldSite, shape=[natms])

      matms = mxatms

      lchr = lenchar2(FldBase%units,C_NULL_CHAR,9)
      lchr = min(lchr,lenchar(FldBase%units,8))
      call cpchr2str(FldBase%units,units,lchr)

      lchr = lenchar2(FldBase%title,C_NULL_CHAR,81)
      lchr = min(lchr,lenchar(FldBase%title,80))
      call cpchr2str(FldBase%title,header,lchr)

      lneut = (FldBase%ineut.eq.1)

!      lchr = lenchar2(fnameC,C_NULL_CHAR,81)
!      lchr = min(lchr,lenchar(fnameC,80))
!      call cpchr2str(fnameC,fname,lchr)

ccc AB: next is for testing (comment out)
!      write(*,*)
!      write(*,*)'field_write(): max nums = ',mxatms,mxsite,mxtmls
!      write(*,*)
!      write(*,*)'field_write(): check FieldSpecs header '
!      write(*,*)'"',(FldBase%title(ich),ich=1,lchr),'"'
!      write(*,*)' =?= '
!      write(*,*)'"',trim(header),'"'
!!      write(*,*)'"',(header(ich:ich),ich=1,lchr),'"',lchr
!      write(*,*)
ccc

!      natms = 0
      matms = 0
      ia = 0

      do im=1,mxtmls

         lchr = lenchar2(MolBase(im)%name,C_NULL_CHAR,81)
         lchr = min(lchr,lenchar(MolBase(im)%name,80))
         call cpchr2str(MolBase(im)%name,mtmolname(im),lchr)

         mtnmolecs(im) = MolBase(im)%nrept
         mtnmolatm(im) = MolBase(im)%nsites
         mtnmolgrp(im) = MolBase(im)%ngroups

ccc AB: next is for testing (comment out)
!         write(*,*)'field_write(): check molecule type "',
!     x     (MolBase(im)%name(ich),ich=1,lchr),'" =?= "',
!     x     trim(mtmolname(im)),'", ',MolBase(im)%id,
!     x     mtnmolecs(im),mtnmolatm(im),mtnmolgrp(im)
ccc
         mtatms = MolBase(im)%nsites

         do ima=1,mtatms
!            natms = natms+1
            ia = ia+1

            lchr = lenchar2(FldSite(ia)%name,C_NULL_CHAR,9)
            lchr = min(lchr,lenchar(FldSite(ia)%name,8))
            call cpchr2str(FldSite(ia)%name,mtatmname(ia),lchr)

            lchr = lenchar2(FldSite(ia)%type,C_NULL_CHAR,9)
            lchr = min(lchr,lenchar(FldSite(ia)%type,8))
            call cpchr2str(FldSite(ia)%type,mtatmtype(ia),lchr)

            mtatmmass(ia) = FldSite(ia)%m
            mtatmchrg(ia) = FldSite(ia)%q

            mtatmnrep(ia) = FldSite(ia)%nrept
            mtatmngrp(ia) = FldSite(ia)%idgrp

            matms = matms+mtnmolecs(im)*mtatmnrep(ia)

         end do

      end do

!      if( natms.ne.mxsite ) then
      if( ia.ne.natms ) then
        write(*,*)
        write(*,*)'field_write(-1): check atoms # in moltypes - ',
     x            ia,' =?= ',natms

        istatC = -1_C_INT
        iflg = -1
      endif

      if( matms.ne.mxatms ) then
        write(*,*)
        write(*,*)'field_write(-2): check atoms # in total - ',
     x            matms,' =?= ',mxatms

        istatC = -2_C_INT
        iflg = -2
      endif

c AB: read the necessary bits from FIELD file

      call fldwrite(idnode,fname,header,units,lneut,mxtmls,
     x  mtmolname,mtnmolecs,mtnmolatm,mtatmname,mtatmtype,
     x  mtatmmass,mtatmchrg,mtatmnrep,matms,istate)

      istatC = int(istate,kind=C_INT)

c AB: deallocate all the big arrays (see above)

      if( iflg.lt.0 .or. istate.lt.0 ) then

        if( iflg.eq.-1 .and. .not.memfld(.false.) ) then
          write(*,*)'---'
          write(*,*)
     x    'field_write(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!          stop
!          return
        elseif( iflg.eq.-2 .and. .not.memall(.false.) ) then

          write(*,*)'---'
          write(*,*)
     x    'field_write(-2): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -8_C_INT
!          stop
!          return
        else
          write(*,*)
          write(*,*)
     x    'field_write(-1/-2): deallocating memory - SUCCESS!'
        endif

        imemall = imemfld+imemfrm

      endif

      return
      end subroutine field_write


      subroutine traj_read(istatC,FrameBaseC,FrameSiteC)
     x bind (C, name="traj_read_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for reading HISTORY file into C-struct
c     equivalent derived types and passing to VOTCA DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: iflgC,istatC

      TYPE(C_PTR), target, intent(inout) :: FrameBaseC
      TYPE(C_PTR), target, intent(inout) :: FrameSiteC

      TYPE(FrameSpecsT), pointer :: FrmBase
      TYPE(FrameSiteT),  pointer :: FrmSite(:)

      character(len=srec) :: fname,header
      character(len=half) :: atmstr

      real(8)       :: cell(9)
      real(8), save :: tstep=0.0,energy=0.0

      integer, save :: idnode=0,iflg=0,imcon=0,levtrj=2,natms=0,nstep=0
      integer       :: istate,lchr,ich,im,ia

      logical, save :: is_mem_new=.false.,istrj=.true.
!      logical, save :: is_mem_new,istrj
!      data is_mem_new/.false./,istrj/.true./

      CALL C_F_POINTER(C_LOC(FrameBaseC), FrmBase)

      fname = 'HISTORY'
      header= ' '
      atmstr= ' '

!      idnode = 0
      istate = 0
      lchr   = 0

!      tstep  = 0.d0
!      nstep  = 0

!      imcon  = 0
!      levtrj = 2

      do ich = 1,9
        cell(ich) = 0.d0
      end do

      iflgC = istatC
      iflg  = int(iflgC)

c AB: iflg > 0  => reading HISTORY header only, to get the params (do NOT close it!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)

c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)
c AB: iflg = -1 => reading HISTORY done after next frame read (close it!)
c AB: iflg < -1 => reading HISTORY done, no frame to be read  (close it!)

c AB: istate = 0 => all right
c AB: istate = 1 => all right, levtrj/keytrj reset to 0      (as found in HISTORY)
c AB: istate = 2 => all right, levtrj/keytrj reset to 0 or 1 (as found in HISTORY)
c AB: istate = 3 => all right, imcon reset to PBC/cell spec found in HISTORY
c AB: istate = 4 => all right, number of site entries reset to that found in HISTORY (<=mxatms)
c AB: istate =-1 => all right, end of HISTORY file reached with last frame read in full
c AB: istate =-2 => error: end of HISTORY file reached before a frame record finished
c AB: istate =-3 => error: HISTORY file not found
c AB: istate =-4 => error: number of site entries in HISTORY is greater than expected (>mxatms)

      if( iflg.gt.0 ) then
c AB: only header will be read - before allocating memory!
c AB: allocate all the needed big arrays

        call IniFrmSpecs(FrmBase)

c AB: read only the header (iflg>0) 
        call frmread
     x  (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x   energy,cell,fratmname,fratmmass,fratmchrg,xyz,vel,frc,istate)

        istatC = int(istate,kind=C_INT)

        if( iflg.lt.-2 ) then
c AB: failed to open or read the file

          istatC = -5_C_INT

          if( .not.memfrm(.false.) ) then

            write(*,*)'---'
            write(*,*)'traj_read(1): Failed deallocating memory!'
            write(*,*)'==='

            istatC = -7_C_INT
!            stop
!            return
          endif

          imemall = imemfld+imemfrm

          return

        elseif( imemfrm.lt.1 .or. mxatms.lt.natms ) then

          mxatms = natms

          if( .not.memfrm(.true.) ) then
            write(*,*)'---'
            write(*,*)
     x        'traj_read(1): Failed (re)allocating memory!'
            write(*,*)'==='

            istatC = -6_C_INT

            FrmBase%nsites = int(natms,kind=C_INT)

!            stop
            return
          else
            write(*,*)
            write(*,*)
     x       'traj_read(1): (re)allocating memory - SUCCESS!'
          endif

        endif

        iflg = 0

        imemall = imemfld+imemfrm

        lchr = len_trim(header)
        call cpstr2chr(header,FrmBase%title,lchr)
        FrmBase%title(lchr+1) = C_NULL_CHAR

        FrmBase%nsites = int(natms,kind=C_INT)
        FrmBase%imcon  = int(imcon,kind=C_INT)
        FrmBase%keytrj = int(levtrj,kind=C_INT)

ccc AB: next is for testing (comment out)
!        write(*,*)
!        write(*,*)'traj_read(1): max nums = ',mxatms,mxsite,mxtmls
!        write(*,*)
!        write(*,*)'traj_read(1): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!        write(*,*)
ccc        
        return

!      if( iflg.gt.0 ) then - done
      else

        natms  = int(FrmBase%nsites)
        imcon  = int(FrmBase%imcon,kind=C_INT)
        levtrj = int(FrmBase%keytrj,kind=C_INT)

        if( imemfrm.lt.1 .or. mxatms.lt.natms ) then

          mxatms = natms

          if ( .not.memfrm(.true.) ) then
            write(*,*)'---'
            write(*,*)
     x        'traj_read(0): Failed (re)allocating memory!'
            write(*,*)'==='

            istatC = -6_C_INT

!            stop
            return
          else
            write(*,*)
            write(*,*)
     x        'traj_read(0): (re)allocating memory - SUCCESS!'
          endif

        endif

      endif

c AB: read the next frame (iflg<1)
      call frmread
     x (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x  energy,cell,fratmname,fratmmass,fratmchrg,xyz,vel,frc,istate)

      istatC = int(istate,kind=C_INT)

      if( istate.gt.0 ) then
c AB: update FrameSpecs if anything has changed (except for the header)

!        lchr = len_trim(header)
!        call cpstr2chr(header,FrmBase%title,lchr)
!        FrmBase%title(lchr+1) = C_NULL_CHAR

        FrmBase%nsites = int(natms,kind=C_INT)
        FrmBase%imcon  = int(imcon,kind=C_INT)
        FrmBase%keytrj = int(levtrj,kind=C_INT)

ccc AB: next is for testing (comment out)
!        lchr = lenchar2(FrmBase%title,C_NULL_CHAR,81)
!        lchr = min(lchr,lenchar(FrmBase%title,80))
!        write(*,*)
!        write(*,*)'traj_read(.): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!        write(*,*)
ccc
      endif

      if( istate.gt.-1 ) then

        call SetCell(FrmBase%cell,cell)

!        FrmBase%energy = real(energy,kind=C_DOUBLE)
        FrmBase%tstep  = real(tstep,kind=C_DOUBLE)
        FrmBase%nstep  = int(nstep,kind=C_INT)
        FrmBase%nsites = int(natms,kind=C_INT)
        FrmBase%imcon  = int(imcon,kind=C_INT)
        FrmBase%keytrj = int(levtrj,kind=C_INT)

        CALL C_F_POINTER(C_LOC(FrameSiteC), FrmSite, shape=[natms])

        call IniFrmSites(FrmSite,natms)

        do ia=1,natms

          atmstr = fratmname(ia)
          lchr = len_trim(atmstr)
          call cpstr2chr(atmstr,FrmSite(ia)%name,lchr)
          call cpstr2chr(atmstr,FrmSite(ia)%type,lchr)
          FrmSite(ia)%name(lchr+1) = C_NULL_CHAR
          FrmSite(ia)%type(lchr+1) = C_NULL_CHAR

          FrmSite(ia)%r(1) = real(xyz(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%r(2) = real(xyz(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%r(3) = real(xyz(3,ia),kind=C_DOUBLE)

         if( levtrj.gt.0 ) then
          
          FrmSite(ia)%v(1) = real(vel(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%v(2) = real(vel(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%v(3) = real(vel(3,ia),kind=C_DOUBLE)

         endif

         if( levtrj.gt.1 ) then
          
          FrmSite(ia)%f(1) = real(frc(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%f(2) = real(frc(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%f(3) = real(frc(3,ia),kind=C_DOUBLE)

         endif

          FrmSite(ia)%m = real(fratmmass(ia),kind=C_DOUBLE)
          FrmSite(ia)%q = real(fratmchrg(ia),kind=C_DOUBLE)

          FrmSite(ia)%id = int(ia,kind=C_INT)

c AB: site id's within a molecule and res/charge-group are to be set elsewhere!
!          FrmSite(ia)%im = 0_c_int
!          FrmSite(ia)%ig = 0_c_int

        end do

      endif

      if( iflg.lt.0 .or. istate.lt.0 ) then

!        if( iflg.lt.-2 ) istatC = -5_C_INT

c AB: deallocate all the big arrays (see above)
        if( .not.memfrm(.false.) ) then
          write(*,*)'---'
          write(*,*)'traj_read(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!          stop
!          return
        else
          write(*,*)
          write(*,*)
     x    'traj_read(-1): deallocating memory - SUCCESS!'
        endif

        imemall = imemfld+imemfrm

      endif

      return
      end subroutine traj_read


      subroutine traj_write(istatC,FrameBaseC,FrameSiteC)
     x bind (C, name="traj_write_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for writing HISTORY file from C-struct
c     equivalent derived types passed from VOTCA to this DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: iflgC,istatC

      TYPE(C_PTR), target, intent(inout) :: FrameBaseC
      TYPE(C_PTR), target, intent(inout) :: FrameSiteC

      TYPE(FrameSpecsT), pointer :: FrmBase
      TYPE(FrameSiteT),  pointer :: FrmSite(:)

      character(len=srec) :: fname,header
      character(len=half) :: atmstr

!      real(8), save :: cell(9),tstep,energy
      real(8)       :: cell(9)
      real(8), save :: tstep=0.0,energy=0.0

      integer, save :: idnode=0,iflg=0,imcon=0,levtrj=2,natms=0,nstep=0
      integer :: istate,lchr,ich,im,ia

      logical, save :: is_mem_new=.false.,istrj=.true.
!      logical, save :: is_mem_new,istrj
!      data is_mem_new/.false./,istrj/.true./

      fname = 'HISTORY_CG'
      header= ' '
      atmstr= ' '

!      idnode = 0
      istate = 0
      lchr   = 0

      iflgC = istatC
      iflg  = int(iflgC)

      CALL C_F_POINTER(C_LOC(FrameBaseC), FrmBase)

      lchr = lenchar2(FrmBase%title,C_NULL_CHAR,81)
      lchr = min(lchr,lenchar(FrmBase%title,80))
      call cpchr2str(FrmBase%title,header,lchr)

      do ich=1,9
        cell(ich) = dble(FrmBase%cell(ich))
      end do

      tstep  = dble(FrmBase%tstep)
      energy = dble(FrmBase%energy)

      nstep  = int(FrmBase%nstep)
      natms  = int(FrmBase%nsites)
      imcon  = int(FrmBase%imcon)
      levtrj = int(FrmBase%keytrj)

      CALL C_F_POINTER(C_LOC(FrameSiteC), FrmSite, shape=[natms])

c AB: allocate all the needed big arrays (see above)

!      is_mem_new = (iflgC.gt.0_C_INT) .or. (mxatms.lt.natms)
      is_mem_new = mxatms.lt.natms

      if( is_mem_new ) then
         mxatms = natms
      endif

      is_mem_new = is_mem_new .or. (imemfrm.lt.1)
      if( is_mem_new ) then

        if ( .not.memfrm(.true.) ) then
          write(*,*)'---'
          write(*,*)
     x       'traj_write(1): Failed (re)allocating memory!'
          write(*,*)'==='
!          stop

          istatC = -6_C_INT

          return
        else
          write(*,*)
          write(*,*)
     x       'traj_write(1): (re)allocating memory - SUCCESS!'
        endif

      endif

      is_mem_new = .false.
      imemall = imemfld+imemfrm

ccc AB: next is for testing (comment out)
!      if(nstep.lt.1001) then
!        write(*,*)
!        write(*,*)'traj_write(.): max nums = ',mxatms,mxsite,mxtmls
!        write(*,*)
!        write(*,*)'traj_write(.): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!      endif
ccc

      do ia=1,natms

        lchr = lenchar2(FrmSite(ia)%name,C_NULL_CHAR,9)
        lchr = min(lchr,lenchar(FrmSite(ia)%name,8))
        call cpchr2str(FrmSite(ia)%name,fratmname(ia),lchr)

!        lchr = lenchar2(FrmSite(ia)%type,C_NULL_CHAR,9)
!        lchr = min(lchr,lenchar(FrmSite(ia)%type,8))
!        call cpchr2str(FrmSite(ia)%type,fratmtype(ia),lchr)

        xyz(1,ia) = dble(FrmSite(ia)%r(1))
        xyz(2,ia) = dble(FrmSite(ia)%r(2))
        xyz(3,ia) = dble(FrmSite(ia)%r(3))

        if( levtrj.gt.0 ) then

          vel(1,ia) = dble(FrmSite(ia)%v(1))
          vel(2,ia) = dble(FrmSite(ia)%v(2))
          vel(3,ia) = dble(FrmSite(ia)%v(3))

        endif

        if( levtrj.gt.1 ) then
          
          frc(1,ia) = dble(FrmSite(ia)%f(1))
          frc(2,ia) = dble(FrmSite(ia)%f(2))
          frc(3,ia) = dble(FrmSite(ia)%f(3))

        endif

        fratmmass(ia) = dble(FrmSite(ia)%m)
        fratmchrg(ia) = dble(FrmSite(ia)%q)

c AB: site id's within a molecule and res/charge-group are to be set elsewhere!
!        FrmSite(ia)%im = 0_c_int
!        FrmSite(ia)%ig = 0_c_int

      end do

      call frmwrite
     x (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x  energy,cell,fratmname,fratmchrg,fratmmass,xyz,vel,frc,istate)

      istatC = int(istate,kind=C_INT)

c AB: deallocate all the big arrays (see above)

      if( iflg.lt.0 .or. istate.lt.0 ) then

        if( .not.memfrm(.false.) ) then
          write(*,*)'---'
          write(*,*)
     x    'traj_write(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!          stop
!          return
        else
          write(*,*)
          write(*,*)
     x    'traj_write(-1): deallocating memory - SUCCESS!'
        endif

        imemall = imemfld+imemfrm

      endif

      return
      end subroutine traj_write


      subroutine conf_read(istatC,FrameBaseC,FrameSiteC)
     x bind (C, name="conf_read_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for reading HISTORY file into C-struct
c     equivalent derived types and passing to VOTCA DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: iflgC,istatC

      TYPE(C_PTR), target, intent(inout) :: FrameBaseC
      TYPE(C_PTR), target, intent(inout) :: FrameSiteC

      TYPE(FrameSpecsT), pointer :: FrmBase
      TYPE(FrameSiteT),  pointer :: FrmSite(:)

      character(len=srec) :: fname,header
      character(len=half) :: atmstr

!      real(8), save :: cell(9),tstep,energy
      real(8)       :: cell(9)
      real(8), save :: tstep=0.0,energy=0.0

      integer, save :: idnode=0,iflg=0,imcon=0,levtrj=2,natms=0,nstep=0
      integer       :: istate,lchr,ich,im,ia

      logical, save :: is_mem_new=.false.,istrj=.false.
!      logical, save :: is_mem_new,istrj
!      data is_mem_new/.false./,istrj/.false./

      CALL C_F_POINTER(C_LOC(FrameBaseC), FrmBase)

      fname = 'CONFIG'
      header= ' '
      atmstr= ' '

!      idnode = 0
      istate = 0
      lchr   = 0

      imcon  = 0
      levtrj = 2

!      tstep  = 0.d0
!      nstep  = 0

!      imcon  = 0
!      levtrj = 2

      do ich = 1,9
        cell(ich) = 0.d0
      end do

      iflgC = istatC
      iflg  = int(iflgC)

c AB: iflg > 0  => reading HISTORY header only, to get the params (do NOT close it!)
c AB: iflg = 0  => reading HISTORY for next frame (do NOT close it!)
c AB: iflg = -1 => reading HISTORY done after next frame read (close it!)
c AB: iflg < -1 => reading HISTORY done, no frame to be read  (close it!)

c AB: istate = 0 => all right
c AB: istate = 1 => all right, levtrj/keytrj reset to 0      (as found in HISTORY)
c AB: istate = 2 => all right, levtrj/keytrj reset to 0 or 1 (as found in HISTORY)
c AB: istate = 3 => all right, imcon reset to PBC/cell spec found in HISTORY
c AB: istate = 4 => all right, number of site entries reset to that found in HISTORY (<=mxatms)
c AB: istate =-1 => all right, end of HISTORY file reached with last frame read in full
c AB: istate =-2 => error: end of HISTORY file reached before a frame record finished
c AB: istate =-3 => error: HISTORY file not found
c AB: istate =-4 => error: number of site entries in HISTORY is greater than expected (>mxatms)

c AB: istate =-5 => error: end of HISTORY file reached while reading the header
c AB: istate =-6 => error: memory (re)allocation failure
c AB: istate =-7 => error: memory deallocation failure

!      if( iflgC.gt.0_C_INT ) then
      if( iflg.gt.0 ) then
c AB: only header will be read - before allocating memory!
c AB: allocate all the needed big arrays

        call IniFrmSpecs(FrmBase)

c AB: read only the header (iflg>0)
        call frmread
     x  (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x   energy,cell,fratmname,fratmmass,fratmchrg,xyz,vel,frc,istate)

        istatC = int(istate,kind=C_INT)

        if( iflg.lt.-2 ) then
c AB: failed to open or read the file

          istatC = -5_C_INT

          if( .not.memfrm(.false.) ) then

            write(*,*)'---'
            write(*,*)'conf_read(1): Failed deallocating memory!'
            write(*,*)'==='

            istatC = -7_C_INT
!            stop
!            return
          endif

          imemall = imemfld+imemfrm

          return

        elseif( imemfrm.lt.1 .or. mxatms.lt.natms ) then

          mxatms = natms

          if( .not.memfrm(.true.) ) then
            write(*,*)'---'
            write(*,*)
     x        'conf_read(1): Failed (re)allocating memory!'
            write(*,*)'==='

            istatC = -6_C_INT

            FrmBase%nsites = int(natms,kind=C_INT)

!            stop
            return
          else
            write(*,*)
            write(*,*)
     x       'conf_read(1): (re)allocating memory - SUCCESS!'
          endif

        endif

        iflg = 0

        imemall = imemfld+imemfrm

        lchr = len_trim(header)
        call cpstr2chr(header,FrmBase%title,lchr)
        FrmBase%title(lchr+1) = C_NULL_CHAR

        FrmBase%energy = real(energy,kind=C_DOUBLE)
        FrmBase%nsites = int(natms,kind=C_INT)
        FrmBase%imcon  = int(imcon,kind=C_INT)
        FrmBase%keytrj = int(levtrj,kind=C_INT)

ccc AB: next is for testing (comment out)
!        write(*,*)
!        write(*,*)'conf_read(1): max nums = ',mxatms,mxsite,mxtmls
!        write(*,*)
!        write(*,*)'conf_read(1): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!        write(*,*)
ccc
        return

!      if( iflg.gt.0 ) then - done
      else

        natms  = int(FrmBase%nsites,kind=C_INT)
        imcon  = int(FrmBase%imcon,kind=C_INT)
        levtrj = int(FrmBase%keytrj,kind=C_INT)

        if( imemfrm.lt.1 .or. mxatms.lt.natms ) then

          mxatms = natms

          if ( .not.memfrm(.true.) ) then
            write(*,*)'---'
            write(*,*)
     x        'conf_read(0): Failed (re)allocating memory!'
            write(*,*)'==='

            istatC = -6_C_INT

!            stop
            return
          else
            write(*,*)
            write(*,*)
     x        'conf_read(0): (re)allocating memory - SUCCESS!'
          endif

        endif

      endif

c AB: read the next frame (iflg<1)
      call frmread
     x (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x  energy,cell,fratmname,fratmmass,fratmchrg,xyz,vel,frc,istate)

      istatC = int(istate,kind=C_INT)

      if( istate.gt.0 ) then
c AB: update FrameSpecs if anything has changed (except for the header)

!        lchr = len_trim(header)
!        call cpstr2chr(header,FrmBase%title,lchr)
!        FrmBase%title(lchr+1) = C_NULL_CHAR

        FrmBase%energy = real(energy,kind=C_DOUBLE)
        FrmBase%nsites = int(natms,kind=C_INT)
        FrmBase%imcon  = int(imcon,kind=C_INT)
        FrmBase%keytrj = int(levtrj,kind=C_INT)

      endif

ccc AB: next is for testing (comment out)
!        lchr = lenchar2(FrmBase%title,C_NULL_CHAR,81)
!        lchr = min(lchr,lenchar(FrmBase%title,80))
!        write(*,*)
!        write(*,*)'conf_read(.): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!        write(*,*)
ccc

      if( istate.gt.-1 ) then

        call SetCell(FrmBase%cell,cell)

c AB: CONFIG/REVCON files do not contain nstep & tstep
!!        FrmBase%nstep = int(nstep,kind=C_INT)
!!        FrmBase%tstep = real(tstep,kind=C_DOUBLE)

!        FrmBase%energy = real(energy,kind=C_DOUBLE)
!        FrmBase%nsites = int(natms,kind=C_INT)
!        FrmBase%imcon  = int(imcon,kind=C_INT)
!        FrmBase%keytrj = int(levtrj,kind=C_INT)

        CALL C_F_POINTER(C_LOC(FrameSiteC), FrmSite, shape=[natms])

        call IniFrmSites(FrmSite,natms)

        do ia=1,natms

          atmstr = fratmname(ia)
          lchr = len_trim(atmstr)
          call cpstr2chr(atmstr,FrmSite(ia)%name,lchr)
          call cpstr2chr(atmstr,FrmSite(ia)%type,lchr)
          FrmSite(ia)%name(lchr+1) = C_NULL_CHAR
          FrmSite(ia)%type(lchr+1) = C_NULL_CHAR

          FrmSite(ia)%r(1) = real(xyz(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%r(2) = real(xyz(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%r(3) = real(xyz(3,ia),kind=C_DOUBLE)

         if( levtrj.gt.0 ) then
          
          FrmSite(ia)%v(1) = real(vel(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%v(2) = real(vel(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%v(3) = real(vel(3,ia),kind=C_DOUBLE)

         endif

         if( levtrj.gt.1 ) then
          
          FrmSite(ia)%f(1) = real(frc(1,ia),kind=C_DOUBLE)
          FrmSite(ia)%f(2) = real(frc(2,ia),kind=C_DOUBLE)
          FrmSite(ia)%f(3) = real(frc(3,ia),kind=C_DOUBLE)

         endif

          FrmSite(ia)%id = int(ia,kind=C_INT)

c AB: site mass & charge are to be set in accord with the force-field (topology)!
!          FrmSite(ia)%m = real(fratmmass(ia),kind=C_DOUBLE)
!          FrmSite(ia)%q = real(fratmchrg(ia),kind=C_DOUBLE)

c AB: site id's within a molecule and res/charge-group are to be set elsewhere!
!          FrmSite(ia)%im = 0_c_int
!          FrmSite(ia)%ig = 0_c_int

        end do

      endif

      if( iflg.lt.0 .or. istate.lt.0 ) then

!        if( iflg.lt.-2 ) istatC = -5_C_INT

c AB: deallocate all the big arrays (see above)
        if( .not.memfrm(.false.) ) then
          write(*,*)'---'
          write(*,*)'conf_read(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!          stop
!          return
        else
          write(*,*)
          write(*,*)
     x    'conf_read(-1): deallocating memory - SUCCESS!'
        endif

        imemall = imemfld+imemfrm

      endif

      return
      end subroutine conf_read


      subroutine conf_write(istatC,FrameBaseC,FrameSiteC)
     x bind (C, name="conf_write_")

c***********************************************************************
c     
c     DL_POLY wrapper-routine for writing HISTORY file from C-struct
c     equivalent derived types passed from VOTCA to this DL_POLY interface
c     
c     copyright - daresbury laboratory 2013
c     author    - Andrey Brukhno  July 2013
c     
c***********************************************************************
!      use iso_c_binding
      use, intrinsic :: iso_c_binding

      use dlp_io_module

      implicit none

      integer(kind=C_INT) :: iflgC,istatC

      TYPE(C_PTR), target, intent(inout) :: FrameBaseC
      TYPE(C_PTR), target, intent(inout) :: FrameSiteC

      TYPE(FrameSpecsT), pointer :: FrmBase
      TYPE(FrameSiteT),  pointer :: FrmSite(:)

      character(len=srec) :: fname,title,header
      character(len=half) :: units,atmstr

      real(8), save :: cell(9),tstep,energy

      integer, save :: idnode,iflg,imcon,levtrj,natms,nstep=0
      integer :: istate,lchr,ich,im,ia

      logical, save :: is_mem_new=.false.,istrj=.false.
!      logical, save :: is_mem_new,istrj
!      data is_mem_new/.false./,istrj/.false./

      fname = 'CONFIG_CG'
      title = ' '
      header= ' '
      units = ' '
      atmstr= ' '

      istate = 0
      idnode = 0
      lchr   = 0

      tstep = 0.d0
      nstep = 0

      iflgC = istatC
      iflg  = int(iflgC)

      CALL C_F_POINTER(C_LOC(FrameBaseC), FrmBase)

      lchr = lenchar2(FrmBase%title,C_NULL_CHAR,81)
      lchr = min(lchr,lenchar(FrmBase%title,80))
      call cpchr2str(FrmBase%title,header,lchr)

      do ich=1,9
        cell(ich) = dble(FrmBase%cell(ich))
      end do

c AB: CONFIG/REVCON files do not contain nstep & tstep
!      tstep  = dble(FrmBase%tstep)
      energy = dble(FrmBase%energy)

!      nstep  = int(FrmBase%nstep)
      natms  = int(FrmBase%nsites)
      imcon  = int(FrmBase%imcon)
      levtrj = int(FrmBase%keytrj)

      CALL C_F_POINTER(C_LOC(FrameSiteC), FrmSite, shape=[natms])

c AB: allocate all the needed big arrays (see above)

!      is_mem_new = (iflgC.gt.0_C_INT) .or. (mxatms.lt.natms)
      is_mem_new = mxatms.lt.natms

      if( is_mem_new ) then
         mxatms = natms
      endif

      is_mem_new = is_mem_new .or. (imemfrm.lt.1)
      if( is_mem_new ) then
        if ( .not.memfrm(.true.) ) then
          write(*,*)'---'
          write(*,*)
     x       'conf_write(1): Failed (re)allocating memory!'
          write(*,*)'==='

          istatC = -6_C_INT

!          stop
          return
        else
          write(*,*)
          write(*,*)
     x       'conf_write(1): (re)allocating memory - SUCCESS!'
        endif
      endif

      is_mem_new = .false.
      imemall = imemfld+imemfrm

ccc AB: next is for testing (comment out)
!      if(nstep.lt.1001) then
!        write(*,*)
!        write(*,*)'conf_write(.): max nums = ',mxatms,mxsite,mxtmls
!        write(*,*)
!        write(*,*)'conf_write(.): check FrameSpecs header'
!        write(*,*)'"',(FrmBase%title(ich),ich=1,lchr),'"'
!     x       ,FrmBase%keytrj,FrmBase%imcon,FrmBase%nsites
!        write(*,*)' =?= '
!        write(*,*)'"',trim(header),'"'
!     x       ,levtrj,imcon,natms
!      endif
ccc

      do ia=1,natms

        lchr = lenchar2(FrmSite(ia)%name,C_NULL_CHAR,9)
        lchr = min(lchr,lenchar(FrmSite(ia)%name,8))
        call cpchr2str(FrmSite(ia)%name,fratmname(ia),lchr)

!        lchr = lenchar2(FrmSite(ia)%type,C_NULL_CHAR,9)
!        lchr = min(lchr,lenchar(FrmSite(ia)%type,8))
!        call cpchr2str(FrmSite(ia)%type,fratmtype(ia),lchr)

        xyz(1,ia) = dble(FrmSite(ia)%r(1))
        xyz(2,ia) = dble(FrmSite(ia)%r(2))
        xyz(3,ia) = dble(FrmSite(ia)%r(3))

        if( levtrj.gt.0 ) then

          vel(1,ia) = dble(FrmSite(ia)%v(1))
          vel(2,ia) = dble(FrmSite(ia)%v(2))
          vel(3,ia) = dble(FrmSite(ia)%v(3))

        endif

        if( levtrj.gt.1 ) then
          
          frc(1,ia) = dble(FrmSite(ia)%f(1))
          frc(2,ia) = dble(FrmSite(ia)%f(2))
          frc(3,ia) = dble(FrmSite(ia)%f(3))

        endif

c AB: CONFIG/REVCON files do not contain masses & charges
!        fratmmass(ia) = dble(FrmSite(ia)%m)
!        fratmchrg(ia) = dble(FrmSite(ia)%q)

c AB: site id's within a molecule and res/charge-group are to be set elsewhere!
!        FrmSite(ia)%im = 0_c_int
!        FrmSite(ia)%ig = 0_c_int

      end do

      call frmwrite
     x (idnode,istrj,iflg,fname,header,imcon,levtrj,natms,nstep,tstep,
     x  energy,cell,fratmname,fratmchrg,fratmmass,xyz,vel,frc,istate)

      istatC = int(istate,kind=C_INT)

c AB: deallocate all the big arrays (see above)

      if( iflg.lt.0 .or. istate.lt.0 ) then

        if( .not.memfrm(.false.) ) then
          write(*,*)'---'
          write(*,*)
     x    'conf_write(-1): Failed deallocating memory!'
          write(*,*)'==='

          istatC = -7_C_INT
!         stop
!         return
        else
          write(*,*)
          write(*,*)
     x    'conf_write(-1): deallocating memory - SUCCESS!'
        endif

       imemall = imemfld+imemfrm

      endif

      return
      end subroutine conf_write

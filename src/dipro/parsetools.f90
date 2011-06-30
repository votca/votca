subroutine parse_mols
use orbitals
use strings
use parameters
implicit none
integer stat
character(len=100) buf
character*10 buff, bufff
character(len=20) :: FMT, rFMT
character*20 hstr
character*20 out(6)

real*8, allocatable :: crapbuf(:)

integer iorb, i,j, iact
integer nlA, k, nrA, nlB, nrB, nshift, nlD, nrD, nseg, imin, count, cnt
logical CP

real*8 expo
integer expoint

integer, external :: get_index

      ! Open files for parsing
      open(10,file="dimer.log", status="old")
      open(11,file="../molA/monomer.log", status="old")
      open(12,file="../molB/monomer.log", status="old")

      ! number of basis functions in molA
      do
         read(11,*,iostat=stat) buf, buff
         if (stat /= 0) exit
         if ( index(buf,'NBasis=') .ne. 0 ) then
            read(buff,*) nA
            exit
         end if
      end do
      close(11)

      write(*,*) "Basis functions in A:", nA

      ! number of basis functions in molB
      do
         read(12,*,iostat=stat) buf, buff
         if (stat /= 0) exit
         if ( index(buf,'NBasis=') .ne. 0 ) then
            read(buff,*) nB
            exit
         end if
      end do
      close(12)

      write(*,*) "Basis functions in B:", nB

      ! number of basis functions in dimer
      do
         read(10,*,iostat=stat) buf, buff
         if (stat /= 0) exit
         if ( index(buf,'NBasis=') .ne. 0 ) then
            read(buff,*) nD
            exit
         end if
      end do
      close(10)

      write(*,*) "Basis functions in D:", nD

      ! Consistency checks, check for CP/noCP basis or ERROR
      if (nD .ne. (nA+nB) ) then
         if ( (nD.ne.nA) .or. (nD .ne. nB)) then
            stop "Count of basis functions doesn't match. Failing."
         else
            CP=.true.
            nshift=0
         end if
      else
         CP=.false.
         nshift=nA
      end if


      write(*,*) "Counterpoise basis?", CP

select case(QCP)
   case ('G')

      
      ! read from molecule A
      nAt = homoAsub+lumoAsub
      allocate(molAcoef(nAt,nD), molAener(nAt))
      molAcoef = 0.0d0
      molAener = 0.0d0
      nlA = nA/5
      nrA = nA-nlA*5
      open(8,file="../molA/fort.7",status='old')
      read(8,'(A8)') FMT
      do i = 1, homoAind(1)+lumoAsub
         read(8,"(a)") buf
         CALL PARSE(buf,'=',out,nseg)
         if (i .lt. homoAind(homoAsub) ) then
            do j = 1, nlA
               read(8,*) buf
            end do
            if (nrA .ne. 0 ) read(8,*) buf
         else
            iact = i-homoAind(homoAsub)+1
            read(out(2),*) molAener(iact)
            do j = 1, nlA
               read(8,FMT) molAcoef(iact,(j-1)*5+1:j*5)
            end do
            if (nrA .ne. 0) then
               write(hstr,'(I1)') nrA
               hstr = trim(hstr)
               rFMT = '('//hstr
               rFMT = trim(rFMT)//'D15.8)'
               read(8,rFMT) (molAcoef(iact,k),k=nA-nrA+1,nA)
            end if
         end if
      end do
      close(8)
      molAener = molAener*conv_Hrt_eV


      write(*,*) "Reading A done!"

      ! read from molecule B
      nBt = homoBsub+lumoBsub
      allocate(molBcoef(nBt,nD), molBener(nBt))
      molBcoef = 0.0d0
      molBener = 0.0d0
      nlB = nB/5
      nrB = nB-nlB*5
      open(8,file="../molB/fort.7",status='old')
      read(8,'(A8)') FMT
      do i = 1, homoBind(1)+lumoBsub
         read(8,'(a)') buf
         CALL PARSE(buf,'=',out,nseg)
         if (i .lt. homoBind(homoBsub) ) then
            do j = 1, nlB
               read(8,*) buf
            end do
            if (nrB .ne. 0 ) read(8,*) buf
         else
            iact = i-homoBind(homoBsub)+1
            read(out(2),*) molBener(iact)
            do j = 1, nlB
               read(8,FMT) molBcoef(iact,nshift+(j-1)*5+1:nshift+j*5)
            end do
            if (nrB .ne. 0) then
               write(hstr,'(I1)') nrB
               hstr = trim(hstr)
               rFMT = '('//hstr
               rFMT = trim(rFMT)//'D15.8)'
               read(8,rFMT) (molBcoef(iact,nshift+k),k=nB-nrB+1,nB)
            end if
         end if
      end do
      close(8)
      molBener = molBener*conv_Hrt_eV

      write(*,*) "Reading B done!"
   case ('GW')

      ! read from molecule A
      nAt = n_QP_A
      allocate(molAcoef(nAt,nD), molAener(nAt))
      molAcoef = 0.0d0
      molAener = 0.0d0
      nlA = nA/5
      nrA = nA-nlA*5
      open(8,file="../molA/fort.7",status='old')
      read(8,'(A8)') FMT
      do i = 1, nAt
         read(8,"(a)") buf
         CALL PARSE(buf,'=',out,nseg)
         iact = i
         read(out(2),*) molAener(iact)
         do j = 1, nlA
            read(8,FMT) molAcoef(iact,(j-1)*5+1:j*5)
         end do
         if (nrA .ne. 0) then
            write(hstr,'(I1)') nrA
            hstr = trim(hstr)
            rFMT = '('//hstr
            rFMT = trim(rFMT)//'D15.8)'
            read(8,rFMT) (molAcoef(iact,k),k=nA-nrA+1,nA)
         end if
      end do
      close(8)
      molAener = molAener*conv_Hrt_eV

      write(*,*) "Reading A done (GW)!"

      ! read from molecule B
      nBt = n_QP_B
      allocate(molBcoef(nBt,nD), molBener(nBt))
      molBcoef = 0.0d0
      molBener = 0.0d0
      nlB = nB/5
      nrB = nB-nlB*5
      open(8,file="../molB/fort.7",status='old')
      read(8,'(A8)') FMT
      do i = 1, nBt
         read(8,'(a)') buf
         CALL PARSE(buf,'=',out,nseg)
         iact = i
         read(out(2),*) molBener(iact)
         do j = 1, nlB
            read(8,FMT) molBcoef(iact,nshift+(j-1)*5+1:nshift+j*5)
         end do
         if (nrB .ne. 0) then
            write(hstr,'(I1)') nrB
            hstr = trim(hstr)
            rFMT = '('//hstr
            rFMT = trim(rFMT)//'D15.8)'
            read(8,rFMT) (molBcoef(iact,nshift+k),k=nB-nrB+1,nB)
         end if
      end do
      close(8)
      molBener = molBener*conv_Hrt_eV

      write(*,*) "Reading B done (GW)!"
      
   end select

      ! read from dimer
   call system("cat fort.7 | grep Alpha | wc -l > NORB")
   open(27,file="NORB", status="old")
   read(27,*) nDt
   close(27)
   allocate( dimcoef(nDt,nD), fock_diag(nDt,nDt) )
      dimcoef   = 0.0d0
      fock_diag = 0.0d0
      nlD       = nD/5
      nrD       = nD-nlD*5
      open(8,file="fort.7",status='old')
      read(8,'(A8)') FMT
      do i = 1, nDt
         read(8,'(a)') buf
         CALL PARSE(buf,'=',out,nseg)
         read(out(2),*) fock_diag(i,i)
         do j = 1, nlD
            read(8,FMT) dimcoef(i,(j-1)*5+1:j*5)
         end do
         if (nrD .ne. 0) then
            write(hstr,'(I1)') nrD
            hstr = trim(hstr)
            rFMT = '('//hstr
            rFMT = trim(rFMT)//'D15.8)'
            read(8,rFMT) (dimcoef(i,k),k=nD-nrD+1,nD)
         end if
      end do
      close(8)
      fock_diag = fock_diag*conv_Hrt_eV


      write(*,*) "Reading D done!"

      ! read overlap matrix from dimer LOG-file (this will be nasty!)
      allocate( S(nD,nD) )
      S = 0.0d0
      open(8,file="dimer.log",status='old')
      do 
         read(8,'(a)',iostat=stat) buf
         if (stat /= 0) stop "Error reading overlap matrix."
         if ( index(buf,'Overlap') .ne. 0 ) exit 
      end do

      j=0
      imin=1
      ! go through all cols
      do while (j .le. nD)
         read(8,'(a)',iostat=stat) buf
         ! reset line counter
         i=imin
         ! go through line
         do while ( i .le. nD )
            read(8,'(a)',iostat=stat) buf
            CALL PARSE(buf,' ',out,nseg)
            read(out(2:nseg),*) (S(i,j+k),k=1,nseg-1)
            read(out(2:nseg),*) (S(j+k,i),k=1,nseg-1)
            i=i+1
         end do
      
         j = j+5
         imin = imin +5

      end do

      
      close(8)




      write(*,*) "Reading S done!"
      write(111,*) S
!   end select

end subroutine parse_mols

integer function get_index(i,j)
implicit none
integer i,j, iact, jact, div

if (j .gt. i ) then
   iact=j
   jact=i
else
   iact = i
   jact =j
end if

div = iact/2

if (mod(iact,2).eq.0) then
   get_index = iact*(div-1)+div+jact
else
   get_index = iact*(div)+jact
end if

end function get_index



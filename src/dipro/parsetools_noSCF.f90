subroutine parse_mols_TNS
use orbitals
use strings
use parameters
use m_mrgrnk
implicit none
integer stat
character(len=100) buf
character*10 buff, bufff
character(len=20) :: FMT, rFMT
character*20 hstr
character*20 out(6)

real*8, allocatable :: crapbuf(:)
integer, allocatable :: ImolA(:)
integer iorb, i,j, iact
integer nlA, k, nrA, nlB, nrB, nshift, nlD, nrD, nseg, imin, count, cnt
logical CP

real*8 expo
integer expoint

integer, external :: get_index

write(*,*) "TURBOMOLE noSCF"
      
! Open files for parsing
open(10,file="mos", status="old")
!!$open(11,file="../molA/mos", status="old")
!!$open(12,file="../molB/mos", status="old")
!!$
!!$      ! get basis set dimensions for molA
!!$      do
!!$         read(11,'(a)',iostat=stat) buf
!!$         if (stat /= 0) exit
!!$         if ( index(buf,'nsaos=') .ne. 0 ) then
!!$            read( buf(index(buf,'nsaos=')+6:),*) nA
!!$            exit
!!$         end if
!!$      end do
!!$      rewind(11)
!!$      write(*,*) nA
!!$
!!$      ! get basis set dimensions for molB
!!$      do
!!$         read(12,'(a)',iostat=stat) buf
!!$         if (stat /= 0) exit
!!$         if ( index(buf,'nsaos=') .ne. 0 ) then
!!$            read( buf(index(buf,'nsaos=')+6:),*) nB
!!$            exit
!!$         end if
!!$      end do
!!$      rewind(12)
!!$      write(*,*) nB

! get basis set dimensions for dimer
do
   read(10,'(a)',iostat=stat) buf
   if (stat /= 0) exit
   if ( index(buf,'nsaos=') .ne. 0 ) then
      read( buf(index(buf,'nsaos=')+6:),*) nD
      exit
   end if
end do
rewind(10)
write(*,*) nD


!!$      ! read orbitals from molA
!!$      nAt = 2+homoAsub+lumoAsub
!!$      allocate(molAcoef(nAt,nA), molAener(nAt))
!!$      molAcoef = 0.0d0
!!$      molAener = 0.0d0
!!$      nlA = nA/4
!!$      nrA = nA-nlA*4
!!$      ! read header
!!$      read(11,'(a)',iostat=stat) buf 
!!$      read(11,'(a)',iostat=stat) buf 
!!$      read(11,'(a)',iostat=stat) buf 
!!$      do i = 1, homoA+1+lumoAsub
!!$         if (i .lt. homoA-homoAsub) then
!!$            read(11,'(a)',iostat=stat) buf 
!!$            do j = 1, nlA
!!$               read(11,'(a)',iostat=stat) buf
!!$            end do
!!$            if (nrA .ne. 0 ) read(11,*) buf
!!$         else
!!$            cnt =0
!!$            iact = i-homoA+homoAsub+1
!!$            read(11,'(a)',iostat=stat) buf 
!!$            if ( index(buf,'ue=') .ne. 0 ) then
!!$               read( buf(index(buf,'ue=')+3:index(buf,'ue=')+16),*) molAener(iact)
!!$            end if
!!$            do j = 1, nlA
!!$               read(11,'(a)',iostat=stat) buf
!!$               do k = 1, 4
!!$                  cnt = cnt +1
!!$                  read( buf((k-1)*20+1:k*20),*) molAcoef(iact,cnt)
!!$               end do
!!$            end do
!!$            if (nrA .ne. 0 ) then
!!$               read(11,'(a)',iostat=stat) buf
!!$               do k = 1, nrA
!!$                  cnt = cnt +1
!!$                  read( buf((k-1)*20+1:k*20),*) molAcoef(iact,cnt)
!!$               end do
!!$            end if
!!$         end if
!!$      end do
!!$      close(11)
!!$      molAener = molAener*conv_Hrt_eV
!!$
!!$      ! read orbitals from molA
!!$      nBt = 2+homoBsub+lumoBsub
!!$      allocate(molBcoef(nBt,nB), molBener(nBt))
!!$      molBcoef = 0.0d0
!!$      molBener = 0.0d0
!!$      nlB = nB/4
!!$      nrB = nB-nlB*4
!!$      ! read header
!!$      read(12,'(a)',iostat=stat) buf
!!$      read(12,'(a)',iostat=stat) buf
!!$      read(12,'(a)',iostat=stat) buf
!!$      do i = 1, homoB+1+lumoBsub
!!$         if (i .lt. homoB-homoBsub) then
!!$            read(12,'(a)',iostat=stat) buf 
!!$            do j = 1, nlB
!!$               read(12,'(a)',iostat=stat) buf
!!$            end do
!!$            if (nrB .ne. 0 ) read(12,*) buf
!!$         else
!!$            cnt =0
!!$            iact = i-homoB+homoBsub+1
!!$            read(12,'(a)',iostat=stat) buf 
!!$            if ( index(buf,'ue=') .ne. 0 ) then
!!$               read( buf(index(buf,'ue=')+3:index(buf,'ue=')+16),*) molBener(iact)
!!$            end if
!!$            do j = 1, nlB
!!$               read(12,'(a)',iostat=stat) buf
!!$               do k = 1, 4
!!$                  cnt = cnt +1
!!$                  read( buf((k-1)*20+1:k*20),*) molBcoef(iact,cnt)
!!$               end do
!!$            end do
!!$            if (nrB .ne. 0 ) then
!!$               read(12,'(a)',iostat=stat) buf
!!$               do k = 1, nrB
!!$                  cnt = cnt +1
!!$                  read( buf((k-1)*20+1:k*20),*) molBcoef(iact,cnt)
!!$               end do
!!$            end if
!!$         end if
!!$      end do
!!$      close(12)
!!$      molBener = molBener * conv_Hrt_eV
!!$

! read from dimer
allocate( dimcoef(nD,nD), dimen(nD) )
dimcoef   = 0.0d0
nlD = nD/4
nrD = nD-nlD*4  
read(10,'(a)',iostat=stat) buf 
do i = 1, nD
   cnt =0
   read(10,'(a)',iostat=stat) buf 
   if ( index(buf,'ue=') .ne. 0 ) then
      read( buf(index(buf,'ue=')+3:index(buf,'ue=')+18),*) dimen(i)
      read( buf(index(buf,'ue=')+20:index(buf,'ue=')+22),*) expoint
      dimen(i) = dimen(i)*10.0d0**dble(expoint)*conv_Hrt_eV
   end if
   do j = 1, nlD
      read(10,'(a)',iostat=stat) buf
      do k = 1, 4
         cnt = cnt +1
         read( buf((k-1)*20+1:k*20),*) dimcoef(i,cnt)
      end do
   end do
   if (nrD .ne. 0 ) then
      read(10,'(a)',iostat=stat) buf
      do k = 1, nrD
         cnt = cnt +1
         read( buf((k-1)*20+1:k*20),*) dimcoef(i,cnt)
      end do
   end if
end do
close(10)
      
! read overlap matrix from dimer LOG-file (this will be nasty!)
allocate( fock_diag(nD,nD) )
fock_diag = 0.0d0
open(8,file="dimer.out",status='old')
do 
   read(8,'(a)',iostat=stat) buf
   if (stat /= 0) stop "Error reading Fock matrix."
   if ( index(buf,'closed-shell fock matrix') .ne. 0 ) exit 
end do

nlD = ((nD*nD-nD)/2+nD)/3
nrD = ((nD*nD-nD)/2+nD)-nlD*3  


! Upper diagonal of overlap matrix is stored in ONE (!) single line
allocate(crapbuf((nD*nD-nD)/2+nD  ))
cnt=0
! rewrite vector into overlap matrix
do i = 1, nlD
   read(8,'(a)',iostat=stat) buf
   do k = 1, 3
      cnt=cnt+1
      read( buf( (k-1)*24+5:(k-1)*24+24 ), *) crapbuf(cnt) 
   end do
end do
if (nrD .ne. 0 ) then
   read(8,'(a)',iostat=stat) buf
   do k = 1, nrD
      cnt = cnt +1
      read( buf( (k-1)*24+5:(k-1)*24+24 ), *) crapbuf(cnt) 
   end do
end if

do i = 1, nD
   do j = 1, nD
      fock_diag(i,j) = crapbuf(get_index(i,j))
   end do
end do
close(8)
fock_diag = fock_diag*conv_Hrt_eV


end subroutine parse_mols_TNS

!!$integer function get_index(i,j)
!!$implicit none
!!$integer i,j, iact, jact, div
!!$
!!$if (j .gt. i ) then
!!$   iact=j
!!$   jact=i
!!$else
!!$   iact = i
!!$   jact =j
!!$end if
!!$
!!$div = iact/2
!!$
!!$if (mod(iact,2).eq.0) then
!!$   get_index = iact*(div-1)+div+jact
!!$else
!!$   get_index = iact*(div)+jact
!!$end if
!!$
!!$end function get_index

! bloody different fockmatrix output for dscf...
subroutine parse_mols_DNS
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

write(*,*) "TURBOMOLE noSCF"
      
! Open files for parsing
open(10,file="mos", status="old")
open(11,file="../molA/mos", status="old")
open(12,file="../molB/mos", status="old")

      ! get basis set dimensions for molA
      do
         read(11,'(a)',iostat=stat) buf
         if (stat /= 0) exit
         if ( index(buf,'nsaos=') .ne. 0 ) then
            read( buf(index(buf,'nsaos=')+6:),*) nA
            exit
         end if
      end do
      rewind(11)
      write(*,*) nA

      ! get basis set dimensions for molB
      do
         read(12,'(a)',iostat=stat) buf
         if (stat /= 0) exit
         if ( index(buf,'nsaos=') .ne. 0 ) then
            read( buf(index(buf,'nsaos=')+6:),*) nB
            exit
         end if
      end do
      rewind(12)
      write(*,*) nB

! get basis set dimensions for dimer
do
   read(10,'(a)',iostat=stat) buf
   if (stat /= 0) exit
   if ( index(buf,'nsaos=') .ne. 0 ) then
      read( buf(index(buf,'nsaos=')+6:),*) nD
      exit
   end if
end do
rewind(10)
write(*,*) nD


!!$      ! read orbitals from molA
      nAt = 2+homoAsub+lumoAsub
      allocate(molAcoef(nAt,nA), molAener(nAt))
      molAcoef = 0.0d0
      molAener = 0.0d0
      nlA = nA/4
      nrA = nA-nlA*4
      ! read header
      read(11,'(a)',iostat=stat) buf 
      read(11,'(a)',iostat=stat) buf 
      read(11,'(a)',iostat=stat) buf 
      do i = 1, homoA+1+lumoAsub
         if (i .lt. homoA-homoAsub) then
            read(11,'(a)',iostat=stat) buf 
            do j = 1, nlA
               read(11,'(a)',iostat=stat) buf
            end do
            if (nrA .ne. 0 ) read(11,*) buf
         else
            cnt =0
            iact = i-homoA+homoAsub+1
            read(11,'(a)',iostat=stat) buf 
            if ( index(buf,'ue=') .ne. 0 ) then
               read( buf(index(buf,'ue=')+3:index(buf,'ue=')+16),*) molAener(iact)
            end if
            do j = 1, nlA
               read(11,'(a)',iostat=stat) buf
               do k = 1, 4
                  cnt = cnt +1
                  read( buf((k-1)*20+1:k*20),*) molAcoef(iact,cnt)
               end do
            end do
            if (nrA .ne. 0 ) then
               read(11,'(a)',iostat=stat) buf
               do k = 1, nrA
                  cnt = cnt +1
                  read( buf((k-1)*20+1:k*20),*) molAcoef(iact,cnt)
               end do
            end if
         end if
      end do
      close(11)
      molAener = molAener*conv_Hrt_eV


      ! read orbitals from molA
      nBt = 2+homoBsub+lumoBsub
      allocate(molBcoef(nBt,nB), molBener(nBt))
      molBcoef = 0.0d0
      molBener = 0.0d0
      nlB = nB/4
      nrB = nB-nlB*4
      ! read header
      read(12,'(a)',iostat=stat) buf
      read(12,'(a)',iostat=stat) buf
      read(12,'(a)',iostat=stat) buf
      do i = 1, homoB+1+lumoBsub
         if (i .lt. homoB-homoBsub) then
            read(12,'(a)',iostat=stat) buf 
            do j = 1, nlB
               read(12,'(a)',iostat=stat) buf
            end do
            if (nrB .ne. 0 ) read(12,*) buf
         else
            cnt =0
            iact = i-homoB+homoBsub+1
            read(12,'(a)',iostat=stat) buf 
            if ( index(buf,'ue=') .ne. 0 ) then
               read( buf(index(buf,'ue=')+3:index(buf,'ue=')+16),*) molBener(iact)
            end if
            do j = 1, nlB
               read(12,'(a)',iostat=stat) buf
               do k = 1, 4
                  cnt = cnt +1
                  read( buf((k-1)*20+1:k*20),*) molBcoef(iact,cnt)
               end do
            end do
            if (nrB .ne. 0 ) then
               read(12,'(a)',iostat=stat) buf
               do k = 1, nrB
                  cnt = cnt +1
                  read( buf((k-1)*20+1:k*20),*) molBcoef(iact,cnt)
               end do
            end if
         end if
      end do
      close(12)
      molBener = molBener * conv_Hrt_eV


! read from dimer
allocate( dimcoef(nD,nD) )
dimcoef   = 0.0d0
nlD = nD/4
nrD = nD-nlD*4  
read(10,'(a)',iostat=stat) buf 
do i = 1, nD
   cnt =0
   read(10,'(a)',iostat=stat) buf 
!   if ( index(buf,'ue=') .ne. 0 ) then
!      read( buf(index(buf,'ue=')+3:index(buf,'ue=')+18),*) fock_diag(i,i)
!      read( buf(index(buf,'ue=')+20:index(buf,'ue=')+22),*) expoint
!      fock_diag(i,i) = fock_diag(i,i)*10.0d0**dble(expoint)
!   end if
   do j = 1, nlD
      read(10,'(a)',iostat=stat) buf
      do k = 1, 4
         cnt = cnt +1
         read( buf((k-1)*20+1:k*20),*) dimcoef(i,cnt)
      end do
   end do
   if (nrD .ne. 0 ) then
      read(10,'(a)',iostat=stat) buf
      do k = 1, nrD
         cnt = cnt +1
         read( buf((k-1)*20+1:k*20),*) dimcoef(i,cnt)
      end do
   end if
end do
close(10)
      
! read overlap matrix from dimer LOG-file (this will be nasty!)
allocate( fock_diag(nD,nD) )
fock_diag = 0.0d0
open(8,file="dimer.out",status='old')
do 
   read(8,'(a)',iostat=stat) buf
   if (stat /= 0) stop "Error reading Fock matrix."
   if ( index(buf,'closed-shell fock matrix') .ne. 0 ) exit 
end do

nlD = ((nD*nD-nD)/2+nD)/4
nrD = ((nD*nD-nD)/2+nD)-nlD*4  


! Upper diagonal of overlap matrix is stored in ONE (!) single line
allocate(crapbuf((nD*nD-nD)/2+nD  ))
cnt=0
! rewrite vector into overlap matrix
do i = 1, nlD
   read(8,'(a)',iostat=stat) buf
   do k = 1, 4
      cnt=cnt+1
      read( buf( (k-1)*25+10:(k-1)*25+25 ), *) crapbuf(cnt) 
   end do
end do
if (nrD .ne. 0 ) then
   read(8,'(a)',iostat=stat) buf
   do k = 1, nrD
      cnt = cnt +1
      read( buf( (k-1)*25+10:(k-1)*25+25 ), *) crapbuf(cnt) 
   end do
end if

do i = 1, nD
   do j = 1, nD
      fock_diag(i,j) = crapbuf(get_index(i,j))
   end do
end do
close(8)
fock_diag = fock_diag*conv_Hrt_eV


end subroutine parse_mols_DNS

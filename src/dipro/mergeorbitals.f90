!
! Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
PROGRAM merge_orbitals
use strings
use m_mrgrnk
implicit none

character(len=20) :: FMT, a_ch, b_ch, c_ch, d_ch, e_ch, nul
integer bla
character(len=20), external :: modnum
character*20 out(6)
integer nseg

real*8 a,b,c,d,e
real*8, allocatable :: molb(:), merged_orbitals(:,:), merged_energies(:), merged_orbitals_ranked(:,:)

integer n_basis_a, n_lines_a, n_rest_a
integer n_basis_b, n_lines_b, n_rest_b
integer n_orbitals_a, n_orbitals_b, n_basis, n_orbitals
integer, allocatable :: rank(:)
integer i,j, cnt
CHARACTER*100 LINEARG  ,fileA, pathname, fileB, fileC, buf, fileD, fileE
CHARACTER*2 CP

! prepare NULL
nul = modnum(0.0d0)

! number of basis functions and orbitals read from command line
CALL getarg(1,LINEARG)
READ(LINEARG,*) n_basis_a
CALL getarg(2,LINEARG)
READ(LINEARG,*) n_orbitals_a
CALL getarg(3,LINEARG)
READ(LINEARG,*) n_basis_b
CALL getarg(4,LINEARG)
READ(LINEARG,*) n_orbitals_b
CALL getarg(5,LINEARG)
READ(LINEARG,*) pathname
CALL getarg(6,LINEARG)
READ(LINEARG,*) CP

CP=trim(CP)

fileA=trim(trim(pathname)//'/molA/fort.7')
fileB=trim(trim(pathname)//'/molB/fort.7')
fileC=trim(trim(pathname)//'/dim/dimer.com')
fileD=trim(trim(pathname)//'/dim/dimer.fchkin')
fileE=trim(trim(pathname)//'/dim/fort.22')

n_orbitals = n_orbitals_a + n_orbitals_b
if ( CP .eq. "CP") then
   n_basis = n_basis_a
else
   n_basis    = n_basis_a + n_basis_b
end if

! read orbitals into matrix
allocate( merged_orbitals(n_orbitals,n_basis),merged_energies(n_orbitals) )
allocate(rank(n_orbitals))
merged_orbitals = 0.0d0
merged_energies = 0.0d0

! how many basis function lines for molecule A
n_lines_a    = n_basis_a/5
n_rest_a     = n_basis_a - 5*n_lines_a
!if (n_rest_a .ne. 0) then
!   n_basis_b  = n_basis_b - (5-n_rest_a)
!end if
n_lines_b  = n_basis_b/5
n_rest_b   = n_basis_b - 5*n_lines_b

write(*,*) 'basis per mol:', n_basis_a
write(*,*) 'lines from molA:', n_lines_a
write(*,*) 'rest from molA:', n_rest_a
write(*,*) ' basis molB:', n_basis_b
write(*,*) ' lines molB:', n_lines_b
write(*,*) ' rest molB:', n_rest_b



allocate( molB(n_basis))

select case(CP)
            
   case DEFAULT

      ! read from molecule A
      open(8,file=fileA,status='old')
      open(9,file=fileC,status='old',position='append')
      read(8,'(A8)') FMT
      write(9,'(A8)') FMT
      do i = 1, n_orbitals_a
         read(8,"(a)") buf
         CALL PARSE(buf,'=',out,nseg)
         read(out(2),*) merged_energies(i)
         cnt =0
         do j = 1, n_lines_a

            read(8,FMT) merged_orbitals(i,cnt+1), &
                        merged_orbitals(i,cnt+2), &
                        merged_orbitals(i,cnt+3), &
                        merged_orbitals(i,cnt+4), &
                        merged_orbitals(i,cnt+5)

            cnt = cnt +5
         end do
         if (n_rest_a == 1) then
            read(8,'(1D15.8)') merged_orbitals(i,cnt+1)
         else if (n_rest_a == 2) then
            read(8,'(2D15.8)') merged_orbitals(i,cnt+1), &  
                               merged_orbitals(i,cnt+2)
         else if (n_rest_a == 3) then
            read(8,'(3D15.8)') merged_orbitals(i,cnt+1), &
                               merged_orbitals(i,cnt+2), &
                               merged_orbitals(i,cnt+3)
         else if (n_rest_a == 4) then
            read(8,'(4D15.8)') merged_orbitals(i,cnt+1), &
                               merged_orbitals(i,cnt+2), &
                               merged_orbitals(i,cnt+3), &
                               merged_orbitals(i,cnt+4)
         end if
      end do
      close(8)

      ! now molecule B
      open(8,file=fileB,status='old')

      read(8,*) FMT
      do i = 1, n_orbitals_b
         read(8,"(a)") buf
         CALL PARSE(buf,'=',out,nseg)
         read(out(2),*) merged_energies(n_orbitals_a+i)
         cnt = n_basis_a
         do j = 1, n_lines_b

            read(8,FMT) merged_orbitals(n_orbitals_a+i,cnt+1), &
                        merged_orbitals(n_orbitals_a+i,cnt+2), &
                        merged_orbitals(n_orbitals_a+i,cnt+3), &
                        merged_orbitals(n_orbitals_a+i,cnt+4), &
                        merged_orbitals(n_orbitals_a+i,cnt+5)

            cnt = cnt +5
         end do
         if (n_rest_b == 1) then
            read(8,'(1D15.8)') merged_orbitals(n_orbitals_a+i,cnt+1)
         else if (n_rest_b == 2) then
            read(8,'(2D15.8)') merged_orbitals(n_orbitals_a+i,cnt+1), &  
                               merged_orbitals(n_orbitals_a+i,cnt+2)
         else if (n_rest_b == 3) then
            read(8,'(3D15.8)') merged_orbitals(n_orbitals_a+i,cnt+1), &
                               merged_orbitals(n_orbitals_a+i,cnt+2), &
                               merged_orbitals(n_orbitals_a+i,cnt+3)
         else if (n_rest_b == 4) then
            read(8,'(4D15.8)') merged_orbitals(n_orbitals_a+i,cnt+1), &
                               merged_orbitals(n_orbitals_a+i,cnt+2), &
                               merged_orbitals(n_orbitals_a+i,cnt+3), &
                               merged_orbitals(n_orbitals_a+i,cnt+4)
         end if

      end do

      close(8)



      CALL  mrgrnk (merged_energies, rank)                  
allocate(merged_orbitals_ranked(n_orbitals,n_basis))
      do i = 1, n_orbitals
         write(9,'(I5,A,D15.8)') i," Alpha MO OE=", merged_energies(rank(i))
         write(9,FMT) merged_orbitals(rank(i),:)
		merged_orbitals_ranked(i,:) = merged_orbitals(rank(i),:)
      end do
      
      close(9)





	  open(9,file=fileD,status='old',position='append')
        write(9,'(5E16.8)') merged_orbitals_ranked
	close(9)

	open(9,file=fileE,status='replace')
	do i = 1, n_orbitals
	write(9,*) i,  merged_energies(rank(i))
	do j = 1, n_basis
		write(9,*) i,j,merged_orbitals_ranked(i,j)
    end do
end do
	close(9)
   end select




 end PROGRAM merge_orbitals

character(len=20) function modnum(x)
implicit none
real*8 x

write(modnum,'(D15.8)') x
modnum = adjustl(modnum)
if (modnum(1:1) == '0') then
   modnum(1:1) = ' '
else if (modnum(1:1) == '-') then
   if (modnum(2:2) == '0') then
      modnum(2:len(modnum)-1) = modnum(3:len(modnum))
   end if
end if

end function modnum
!!$
!!$   case('CP')
!!$      open(8,file=fileA,status='old')
!!$      open(9,file=fileC,status='old',position='append')
!!$      read(8,'(A8)') FMT
!!$      write(9,'(A8)') FMT
!!$      do i = 1, n_orbitals
!!$         read(8,'(I5)') bla
!!$         write(9,'(I5)') bla
!!$         do j = 1, n_lines
!!$            read(8,FMT) a,b,c,d,e
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            d_ch = modnum(d)
!!$            e_ch = modnum(e)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch),' ',trim(d_ch),' ',trim(e_ch)
!!$         end do
!!$         if (n_rest == 1) then
!!$            read(8,'(1D15.8)') a
!!$            a_ch = modnum(a)
!!$            write(9,*) trim(a_ch)
!!$         else if (n_rest == 2) then
!!$            read(8,'(2D15.8)') a, b
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch)
!!$         else if (n_rest == 3) then
!!$            read(8,'(3D15.8)') a, b, c
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch)
!!$         else if (n_rest == 4) then
!!$            read(8,'(4D15.8)') a, b, c, d
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            d_ch = modnum(d)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch),' ',trim(d_ch)
!!$         end if
!!$      end do
!!$      close(8)
!!$
!!$      open(8,file=fileB,status='old')
!!$      read(8,*) FMT
!!$      do i = 1, n_orbitals
!!$         read(8,'(I5)') bla
!!$         write(9,'(I5)') bla+n_orbitals
!!$         do j = 1, n_lines
!!$            read(8,FMT) a,b,c,d,e
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            d_ch = modnum(d)
!!$            e_ch = modnum(e)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch),' ',trim(d_ch),' ',trim(e_ch)
!!$         end do
!!$         if (n_rest == 1) then
!!$            read(8,'(1D15.8)') a
!!$            a_ch = modnum(a)
!!$            write(9,*) trim(a_ch)
!!$         else if (n_rest == 2) then
!!$            read(8,'(2D15.8)') a, b
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch)
!!$         else if (n_rest == 3) then
!!$            read(8,'(3D15.8)') a, b, c
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch)
!!$         else if (n_rest == 4) then
!!$            read(8,'(4D15.8)') a, b, c, d
!!$            a_ch = modnum(a)
!!$            b_ch = modnum(b)
!!$            c_ch = modnum(c)
!!$            d_ch = modnum(d)
!!$            write(9,*) trim(a_ch),' ',trim(b_ch),' ',trim(c_ch),' ',trim(d_ch)
!!$         end if
!!$      end do
!!$      close(8)

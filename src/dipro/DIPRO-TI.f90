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
module orbitals
real*8, allocatable :: molAcoef(:,:), molBcoef(:,:), dimcoef(:,:)
real*8, allocatable :: molAener(:), molBener(:), dimen(:),fock_diag(:,:), S(:,:)
character*3 QCP
character*100 plabel
integer homoA, homoB, homoAsub, homoBsub
integer lumoA, lumoB, lumoAsub, lumoBsub
integer nA, nB, nD, nAt, nBt, nDt
integer n_QP_A, n_QP_B, n_QP_D
real*8, allocatable :: J_hole(:,:)
real*8, allocatable :: J_elec(:,:)
integer,allocatable :: homoAind(:), homoBind(:), lumoAind(:), lumoBind(:)
real*8, allocatable :: homoAboltz(:),homoBboltz(:),lumoAboltz(:),lumoBboltz(:)
real*8 J_sq_degen_h, J_sq_degen_e, J_sq_boltz_h, J_sq_boltz_e

real*8, allocatable :: QPS_molA(:,:), QPS_molB(:,:), QPS_dim(:,:)
real*8, allocatable :: QPE_molA(:), QPE_molB(:), QPE_dim(:,:)
real*8, allocatable :: QP_dimer_ProA(:,:), QP_dimer_ProB(:,:)
real*8, allocatable :: QPA_DimBS(:,:), QPB_DimBS(:,:)
integer n_dim, n_dim_A, n_dim_B

end module orbitals

module parameters
real*8, parameter :: conv_Hrt_eV=27.21138386d0
end module parameters



PROGRAM DIPRO_TI
use orbitals
IMPLICIT NONE

character*100 linearg
real*8, allocatable :: dimer_Pro(:,:) 
real*8, allocatable :: PsiA_DimBS(:,:) 
real*8, allocatable :: PsiB_DimBS(:,:) 
real*8, allocatable :: helpmat(:,:)
real*8, allocatable :: JMAT(:,:)
real*8, allocatable :: SMAT(:,:)
real*8, allocatable :: diagS(:,:), eigvec_t(:,:), transmat(:,:), JMAT_eff(:,:)
real*8, allocatable :: eigval(:),work(:),eigval_sq(:)
real*8 TI_hole, TI_elec, sum, tread
integer i,  info,lwork, j, arg,k,cnt




! number of basis functions and orbitals read from command line
arg=1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) homoAsub

if (homoAsub.ne.0) then
   allocate(homoAind(homoAsub))
   do i = 1, homoAsub
      arg=arg+1
      CALL getarg(arg,LINEARG)
      READ(LINEARG,*) homoAind(i)
   end do
end if

arg=arg+1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) homoBsub
if (homoBsub.ne.0) then
   allocate(homoBind(homoBsub))
   do i = 1, homoBsub
      arg=arg+1
      CALL getarg(arg,LINEARG)
      READ(LINEARG,*) homoBind(i)
   end do
end if


arg=arg+1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) lumoAsub
if (lumoAsub.ne.0) then
   allocate(lumoAind(lumoAsub))
   do i = 1, lumoAsub
      arg=arg+1
      CALL getarg(arg,LINEARG)
      READ(LINEARG,*) lumoAind(i)
   end do
end if


arg=arg+1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) lumoBsub
if (lumoBsub.ne.0) then
   allocate(lumoBind(lumoBsub))
   do i = 1, lumoBsub
      arg=arg+1
      CALL getarg(arg,LINEARG)
      READ(LINEARG,*) lumoBind(i)
   end do
end if

arg=arg+1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) QCP

arg=arg+1
CALL getarg(arg,LINEARG)
READ(LINEARG,*) plabel

QCP=trim(QCP)
plabel=trim(plabel)

! Check validity of QCP parameter
if (  (QCP .ne. "T") .and. (QCP .ne. "G") .and. (QCP .ne. 'GW') ) then
   write(*,*) "Quantum-chemistry program not supported!"
   write(*,*) "Specify G for Gaussian09,"
   write(*,*) "        T for Turbomole"
   stop
end if

! Read orbital coefficients and orbital energies from QC software output

if (QCP .eq. "T") then
   CALL parse_mols_TNS
   
   allocate(J_hole(homoAsub,homoBsub),J_elec(lumoAsub,lumoBsub))


   do i = 1, homoAsub
      do j = 1, homoBsub
         J_hole(i,j) = dot_product(dimcoef(homoAind(i),:),matmul(fock_diag,dimcoef(homoBind(j),:)))
      end do
   end do

   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_elec(i,j) = dot_product(dimcoef(lumoAind(i),:),matmul(fock_diag,dimcoef(lumoBind(j),:)))
      end do
   end do

   allocate(homoAboltz(homoAsub),homoBboltz(homoBsub)) 
   allocate(lumoAboltz(lumoAsub),lumoBboltz(lumoBsub)) 

   sum = 0.0d0
   do i = 1, homoAsub
      homoAboltz(i)=exp((dimen(homoAind(i))-dimen(homoAind(1)))/0.0256d0)
      sum = sum + homoAboltz(i)*homoAboltz(i)
   end do
   homoAboltz=homoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, homoBsub
      homoBboltz(i)=exp((dimen(homoBind(i))-dimen(homoBind(1)))/0.0256d0)
      sum = sum + homoBboltz(i)*homoBboltz(i)
   end do
   homoBboltz=homoBboltz/sqrt(sum)


   sum = 0.0d0
   do i = 1, lumoAsub
      lumoAboltz(i)=exp(-(dimen(lumoAind(i))-dimen(lumoAind(1)))/0.0256d0)
      sum = sum + lumoAboltz(i)*lumoAboltz(i)
   end do
   lumoAboltz=lumoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, lumoBsub
      lumoBboltz(i)=exp(-(dimen(lumoBind(i))-dimen(lumoBind(1)))/0.0256d0)
      sum = sum + lumoBboltz(i)*lumoBboltz(i)
   end do
   lumoBboltz=lumoBboltz/sqrt(sum)


   J_sq_degen_h = 0.0d0
   J_sq_boltz_h = 0.0d0
   do i = 1, homoAsub
      do j = 1, homoBsub
         J_sq_degen_h = J_sq_degen_h + J_hole(i,j)**2.0d0
         J_sq_boltz_h = J_sq_boltz_h + (homoAboltz(i)*J_hole(i,j)*homoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_h = J_sq_degen_h/dble(homoAsub*homoBsub)

   J_sq_degen_e = 0.0d0
   J_sq_boltz_e = 0.0d0
   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_sq_degen_e = J_sq_degen_e + J_elec(i,j)**2.0d0
         J_sq_boltz_e = J_sq_boltz_e + (lumoAboltz(i)*J_elec(i,j)*lumoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_e = J_sq_degen_e/dble(lumoAsub*lumoBsub)


   ! Call data to XML output function (TO BE COMPLETED!)
   CALL XMLoutput

else if ( (QCP .eq. "G") ) then
   CALL parse_mols
   write(*,*) "Reading done!"

   ! Take mo - the alpha electrons to get matrix in correct order
   allocate( dimer_Pro(nD,nDt) )
   dimer_Pro = transpose(matmul(dimcoef,S))

  !  DIPRO.x 8 8 7 6 5 4 3 2 1 8 8 7 6 5 4 3 2 1 6 9 10 11 12 13 14 6 9 10 11 12 13 14 G ma
  ! SMatrix of a dimer
  do i = 1, nD
  do j = 1, nD
    write(100,*) i,j, S(i,j)
  end do	
  end do


   ! Determine orbitals of monomers in the basis set of the dimer
   allocate( PsiA_DimBS(nAt,nDt), PsiB_DimBS(nBt,nDt) )
   PsiA_DimBS = matmul(molAcoef, dimer_Pro)
   PsiB_DimBS = matmul(molBcoef, dimer_Pro)


  do i = 1, nD
  do j = 1, nDt
    write(102,*) i,j, dimer_Pro(i,j)
  end do	
  end do

  do i = 1, nAt
  do j = 1, nDt
    write(101,*) i,j, PsiA_DimBS(i,j)
  end do	
  end do


   allocate(helpmat(nBt,nDt))
   helpmat=matmul( PsiB_DimBS, fock_diag)

   ! Determine transfer matrix JAB and site energies JAA, JBB
   allocate( JMAT(nAt+nBt,nAt+nBt), SMAT(nAt+nBt,nAt+nBt) )
   JMAT(1:nAt,1:nAt) = matmul(matmul( PsiA_DimBS, fock_diag) , transpose(PsiA_DimBS) )
   JMAT(nAt+1:nAt+nBt,1:nAt) =   matmul(helpmat , transpose(PsiA_DimBS) )
   JMAT(1:nAt,nAt+1:nAt+nBt) = transpose( JMAT(nAt+1:nAt+nBt,1:nAt))
   JMAT(nAt+1:nAt+nBt,nAt+1:nAt+nBt) = matmul(helpmat , transpose(PsiB_DimBS) )


  do i = 1, nAt+nbT
  do j = 1, nAt+nbT
    write(122,*) i,j, JMAT(i,j)
  end do	
  end do

   ! Determine overlap analogous to J
   SMAT(1:nAt,1:nAt) = matmul(PsiA_DimBS , transpose(PsiA_DimBS) )
   SMAT(nAt+1:nAt+nBt,nAt+1:nAt+nBt)= matmul(PsiB_DimBS , transpose(PsiB_DimBS) )
   SMAT(nAt+1:nAt+nBt,1:nAt) = matmul(PsiB_DimBS , transpose(PsiA_DimBS) )
   SMAT(1:nAt,nAt+1:nAt+nBt) = transpose( SMAT(nAt+1:nAt+nBt,1:nAt))

  do i = 1, nAt+nbT
  do j = 1, nAt+nbT
    write(120,*) i,j, SMAT(i,j)
  end do
  end do

   ! Determine eigenvales (evalu) and eigenvectors (evec) of overlap matrix
   lwork=3*(nAt+nBt)
   allocate(eigval(nAt+nBt),work(lwork),eigval_sq(nAt+nBt))
   call DSYEV('V','U',nAt+nBt,SMAT,nAt+nBt,eigval,work,lwork,info)

   ! Construct orthogonalization matrix (transmat)
   allocate( diagS(nAt+nBt,nAt+nBt), eigvec_t(nAt+nBt,nAt+nBt), transmat(nAt+nBt,nAt+nBt) )
   !eigval_sq     = sqrt(eigval)
   diagS = 0.0d0
   do i = 1, nAt+nBt
      diagS(i,i) = 1.0d0/sqrt(eigval(i))
   end do
   eigvec_t = transpose(SMAT)
   transmat = matmul(SMAT, matmul(diagS,eigvec_t))


   ! Calculate effective JMAT
   allocate(JMAT_eff(nAt+nBt,nAt+nBt))
   JMAT_eff = matmul( transmat, matmul(JMAT,transmat) )
   
   allocate(J_hole(homoAsub,homoBsub),J_elec(lumoAsub,lumoBsub))

   J_hole = JMAT_eff(1:homoAsub,nAt+1:nAt+homoBsub)
   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_elec(i,j) = JMAT_eff(homoAsub+i,nAt+homoBsub+j)
      end do
   end do

   ! Degenercay runs 
   allocate(homoAboltz(homoAsub),homoBboltz(homoBsub)) 
   allocate(lumoAboltz(lumoAsub),lumoBboltz(lumoBsub)) 

   sum = 0.0d0
   do i = 1, homoAsub
      homoAboltz(i)=exp((molAener(i)-molAener(1))/0.0256d0)
      sum = sum + homoAboltz(i)*homoAboltz(i)
   end do
   homoAboltz=homoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, homoBsub
      homoBboltz(i)=exp((molBener(i)-molBener(1))/0.0256d0)
      sum = sum + homoBboltz(i)*homoBboltz(i)
   end do
   homoBboltz=homoBboltz/sqrt(sum)


   sum = 0.0d0
   do i = 1, lumoAsub
      lumoAboltz(i)=exp(-(molAener(i+homoAsub)-molAener(1+homoAsub))/0.0256d0)
      sum = sum + lumoAboltz(i)*lumoAboltz(i)
   end do
   lumoAboltz=lumoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, lumoBsub
      lumoBboltz(i)=exp(-(molBener(i+homoBsub)-molBener(1+homoBsub))/0.0256d0)
      sum = sum + lumoBboltz(i)*lumoBboltz(i)
   end do
   lumoBboltz=lumoBboltz/sqrt(sum)


   J_sq_degen_h = 0.0d0
   J_sq_boltz_h = 0.0d0
   do i = 1, homoAsub
      do j = 1, homoBsub
         J_sq_degen_h = J_sq_degen_h + J_hole(i,j)**2.0d0
         J_sq_boltz_h = J_sq_boltz_h + (homoAboltz(i)*J_hole(i,j)*homoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_h = J_sq_degen_h/dble(homoAsub*homoBsub)

   J_sq_degen_e = 0.0d0
   J_sq_boltz_e = 0.0d0
   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_sq_degen_e = J_sq_degen_e + J_elec(i,j)**2.0d0
         J_sq_boltz_e = J_sq_boltz_e + (lumoAboltz(i)*J_elec(i,j)*lumoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_e = J_sq_degen_e/dble(lumoAsub*lumoBsub)

   ! Call data to XML output function (TO BE COMPLETED!)
   CALL XMLoutput

else if (QCP.eq."GW") then

   ! Now, we have to read the energies and expansion coefficients of th
   ! quasiparticle states from the respective diag_qp.out files.

   ! Open files for parsing
   open(10,file="diag_QP.dat", status="old")
   open(11,file="../molA/diag_QP.dat", status="old")
   open(12,file="../molB/diag_QP.dat", status="old")

   write(*,*) "Entering GW extension"

   ! molA, only read the QP states we need 
   read(11,*) n_QP_A
   allocate( QPS_molA(homoAsub+lumoAsub,n_QP_A), QPE_molA(n_QP_A) )
   cnt=0
   do i = 1, maxval(lumoAind)
      read(11,*) k, tread
      if ( i .ge. minval(homoAind)) then
         cnt = cnt+1
         QPE_molA(cnt) = tread
      end if
      do j = 1, n_QP_A
         read(11,*) k, tread
         if (i .ge. minval(homoAind)) then
            QPS_molA(cnt,j) = tread
         end if
      end do
   end do
   close(11)

   QPE_molA = QPE_molA*13.605d0
   write(*,*) "read from molA", n_QP_A

   ! molB
   read(12,*) n_QP_B
   allocate( QPS_molB(homoBsub+lumoBsub,n_QP_B), QPE_molB(n_QP_B) )
   cnt=0
   do i = 1, maxval(lumoBind)
      read(12,*) k, tread
      if (i .ge. minval(homoBind)) then
         cnt=cnt+1
         QPE_molB(cnt) = tread
      end if
      do j = 1, n_QP_B
         read(12,*) k, tread
         if (i.ge.minval(homoBind)) then
            QPS_molB(cnt,j) = tread
         end if
      end do
   end do
   close(12)
   QPE_molB = QPE_molB*13.605d0

   write(*,*) "read from molB", n_QP_B

   ! dimer
   read(10,*) n_QP_D
   allocate( QPS_dim(n_QP_D,n_QP_D), QPE_dim(n_QP_D,n_QP_D) )
   QPS_dim = 0.0d0
   QPE_dim = 0.0d0
   do i = 1, n_QP_D
      read(10,*) k, QPE_dim(i,i)
      do j = 1, n_QP_D
         read(10,*) k, QPS_dim(i,j)
      end do
   end do
   close(10)
   QPE_dim = QPE_dim*13.605d0

   write(*,*) "read from dimer", n_QP_D

   
   ! NOW, the parsing of the KS data must be adjusted to fit the GW states
   CALL parse_mols
   write(*,*) "Reading done!"

   ! Take mo - the alpha electrons to get matrix in correct order
   allocate( dimer_Pro(nD,nDt) )
   dimer_Pro = transpose(matmul(dimcoef,S))

   ! Determine orbitals of monomers in the basis set of the dimer
   allocate( PsiA_DimBS(nAt,nDt), PsiB_DimBS(nBt,nDt) )
   PsiA_DimBS = matmul(molAcoef, dimer_Pro)
   PsiB_DimBS = matmul(molBcoef, dimer_Pro)




   ! Take dimer QP states and multiply by overlap matrix <KS_A|KS_D>
   allocate( QP_dimer_ProA(n_QP_A,n_QP_D), QP_dimer_ProB(n_QP_B,n_QP_D) )
   QP_dimer_ProA = matmul(PsiA_DimBS(:,1:n_QP_D),transpose(QPS_dim))
   QP_dimer_ProB = matmul(PsiB_DimBS(:,1:n_QP_D),transpose(QPS_dim))
   ! QP_dimer_ProX contain the terms (mat(S)xvec(d)_i) on the columns

   write(*,*) "calculated QP_diemr_ProX"

!   write(*,*) QP_dimer_ProB


   ! Determine orbitals of monomers in the basis set of the dimer
   allocate( QPA_DimBS(homoAsub+lumoAsub,n_QP_D), QPB_DimBS(homoBsub+lumoBsub,n_QP_D) )
   QPA_DimBS = matmul(QPS_molA, QP_dimer_ProA)
   QPB_DimBS = matmul(QPS_molB, QP_dimer_ProB)
   ! QPX_DimBS(i,j) contain the projection of the i-th QP state of molX
   !                on the j-th QP state of the dimer

   write(*,*) "calculated QPX_DimBS"

!   write(*,*) QPB_DimBS

   allocate(helpmat(homoBsub+lumoBsub,n_QP_D))
   helpmat=matmul( QPB_DimBS, QPE_dim)

   write(*,*) "calculated helpmat"
!   write(*,*) helpmat

   ! Determine transfer matrix JAB and site energies JAA, JBB
!   allocate( JMAT(n_QP_A+n_QP_B,n_QP_A+n_QP_B), SMAT(n_QP_A+n_QP_B,n_QP_A+n_QP_B))
   n_dim   = homoAsub+lumoAsub+homoBsub+lumoBsub
   n_dim_A = homoAsub+lumoAsub
   n_dim_B = homoBsub+lumoBsub

   allocate( JMAT(n_dim,n_dim), SMAT(n_dim,n_dim))

   JMAT(1:n_dim_A,1:n_dim_A) = matmul(matmul( QPA_DimBS, QPE_dim),transpose(QPA_DimBS) )
   JMAT(n_dim_A+1:n_dim_A+n_dim_B,1:n_dim_A) = matmul(helpmat , transpose(QPA_DimBS) )
   JMAT(1:n_dim_A,n_dim_A+1:n_dim_A+n_dim_B) = transpose( JMAT(n_dim_A+1:n_dim_A+n_dim_B,1:n_dim_A))
   JMAT(n_dim_A+1:n_dim_A+n_dim_B,n_dim_A+1:n_dim_A+n_dim_B) = matmul(helpmat , transpose(QPB_DimBS) )

   write(*,*) "calculated JMAT"
!   write(*,*) JMAT

   ! Determine overlap analogous to J
   SMAT(1:n_dim_A,1:n_dim_A) = matmul(QPA_DimBS , transpose(QPA_DimBS) )
   SMAT(n_dim_A+1:n_dim_A+n_dim_B,n_dim_A+1:n_dim_A+n_dim_B)= matmul(QPB_DimBS , transpose(QPB_DimBS) )
   SMAT(n_dim_A+1:n_dim_A+n_dim_B,1:n_dim_A) = matmul(QPB_DimBS , transpose(QPA_DimBS) )
   SMAT(1:n_dim_A,n_dim_A+1:n_dim_A+n_dim_B) = transpose( SMAT(n_dim_A+1:n_dim_A+n_dim_B,1:n_dim_A))

   write(*,*) "calculated SMAT"
!   write(*,*) SMAT

   ! Determine eigenvales (evalu) and eigenvectors (evec) of overlap matrix
   lwork=3*(n_dim)
   allocate(eigval(n_dim),work(lwork),eigval_sq(n_dim))
   call DSYEV('V','U',n_dim,SMAT,n_dim,eigval,work,lwork,info)

!   write(*,*) info

   ! Construct orthogonalization matrix (transmat)
   allocate( diagS(n_dim,n_dim), eigvec_t(n_dim,n_dim), transmat(n_dim,n_dim) )
   !eigval_sq     = sqrt(eigval)
   diagS = 0.0d0
   do i = 1, n_dim
      diagS(i,i) = 1.0d0/sqrt(eigval(i))
   end do
   eigvec_t = transpose(SMAT)
   transmat = matmul(SMAT, matmul(diagS,eigvec_t))

   ! Calculate effective JMAT
   allocate(JMAT_eff(n_dim,n_dim))
   JMAT_eff = matmul( transmat, matmul(JMAT,transmat) )
   
   allocate(J_hole(homoAsub,homoBsub),J_elec(lumoAsub,lumoBsub))

   J_hole = JMAT_eff(1:homoAsub,n_dim_A+1:n_dim_A+homoBsub)
   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_elec(i,j) = JMAT_eff(homoAsub+i,n_dim_A+homoBsub+j)
      end do
   end do

   write(*,*) "holes A->B:", JMAT_eff(homoAsub,n_dim_A+homoBsub)
   write(*,*) "holes B->A:", JMAT_eff(n_dim_A+homoBsub,homoAsub)
   write(*,*) "elecs A->B:", JMAT_eff(homoAsub+1,n_dim_A+homoBsub+1)
   write(*,*) "elecs B->A:", JMAT_eff(n_dim_A+homoBsub+1,homoAsub+1)

   ! Degenercay runs 
   allocate(homoAboltz(homoAsub),homoBboltz(homoBsub)) 
   allocate(lumoAboltz(lumoAsub),lumoBboltz(lumoBsub)) 

   sum = 0.0d0
   do i = 1, homoAsub
      homoAboltz(i)=exp((QPE_molA(i)-QPE_molA(1))/0.0256d0)
      sum = sum + homoAboltz(i)*homoAboltz(i)
   end do
   homoAboltz=homoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, homoBsub
      homoBboltz(i)=exp((QPE_molB(i)-QPE_molB(1))/0.0256d0)
      sum = sum + homoBboltz(i)*homoBboltz(i)
   end do
   homoBboltz=homoBboltz/sqrt(sum)


   sum = 0.0d0
   do i = 1, lumoAsub
      lumoAboltz(i)=exp(-(QPE_molA(i+homoAsub)-QPE_molA(1+homoAsub))/0.0256d0)
      sum = sum + lumoAboltz(i)*lumoAboltz(i)
   end do
   lumoAboltz=lumoAboltz/sqrt(sum)

   sum = 0.0d0
   do i = 1, lumoBsub
      lumoBboltz(i)=exp(-(QPE_molB(i+homoBsub)-QPE_molB(1+homoBsub))/0.0256d0)
      sum = sum + lumoBboltz(i)*lumoBboltz(i)
   end do
   lumoBboltz=lumoBboltz/sqrt(sum)


   J_sq_degen_h = 0.0d0
   J_sq_boltz_h = 0.0d0
   do i = 1, homoAsub
      do j = 1, homoBsub
         J_sq_degen_h = J_sq_degen_h + J_hole(i,j)**2.0d0
         J_sq_boltz_h = J_sq_boltz_h + (homoAboltz(i)*J_hole(i,j)*homoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_h = J_sq_degen_h/dble(homoAsub*homoBsub)

   J_sq_degen_e = 0.0d0
   J_sq_boltz_e = 0.0d0
   do i = 1, lumoAsub
      do j = 1, lumoBsub
         J_sq_degen_e = J_sq_degen_e + J_elec(i,j)**2.0d0
         J_sq_boltz_e = J_sq_boltz_e + (lumoAboltz(i)*J_elec(i,j)*lumoBboltz(j))**2.0d0
      end do
   end do
   J_sq_degen_e = J_sq_degen_e/dble(lumoAsub*lumoBsub)


   call XMLoutput


end if
END PROGRAM DIPRO_TI



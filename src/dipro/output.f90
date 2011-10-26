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
subroutine XMLoutput
use strings
use orbitals
implicit none
character*100 printval
integer n, i, j 
character*40 out(20), pval, fileN, DLAB, ival, jval
logical ex
real*8 dist(4)

open(16,file="../TI.xml", status="replace")

printval='<pair name="'//trim(plabel)//'">'
write(16,*) trim(printval)

! workaround for distances (-> rates)
fileN = '../'//trim(plabel)//'.dist'
inquire(file=fileN, exist=ex) 
if ( ex ) then
   open(15,file=fileN,status='old')
   read(15,*) DLAB, dist(:)
   close(15)
else
   dist = 0.0d0
end if

write(16,*) '   <parameters>'
write(pval,*) homoAind(1)
CALL REMOVESP(pval)
printval='      <HOMO_A>'// trim(pval) // '</HOMO_A>'
write(16,*) trim(printval)

write(pval,*) homoAsub
CALL REMOVESP(pval)
printval='      <NoccA>' // trim(pval) // '</NoccA>'
write(16,*) trim(printval)

write(pval,*) lumoAind(1)
CALL REMOVESP(pval)
printval='      <LUMO_A>'// trim(pval) // '</LUMO_A>'
write(16,*) trim(printval)

write(pval,*) lumoAsub
CALL REMOVESP(pval)
printval='      <NvirtA>' // trim(pval) // '</NvirtA>'
write(16,*) trim(printval)

write(pval,*) homoBind(1)
CALL REMOVESP(pval)
printval='      <HOMO_B>' // trim(pval) // '</HOMO_B>'
write(16,*) trim(printval)

write(pval,*) homoBsub
CALL REMOVESP(pval)
printval='      <NoccB>' // trim(pval) // '</NoccB>'
write(16,*) trim(printval)

write(pval,*) lumoBind(1)
CALL REMOVESP(pval)
printval='      <LUMO_B>' // trim(pval) // '</LUMO_B>'
write(16,*) trim(printval)

write(pval,*) lumoBsub
CALL REMOVESP(pval)
printval='      <NvirtB>' // trim(pval) // '</NvirtB>'
write(16,*) trim(printval)


write(16,*) '   </parameters>'
write(16,*) '    <transport name="hole">'
write(16,*) '        <channel name="single">'
if (QCP.eq."T") then
   write(pval,*) J_hole(1,1)
else
   write(pval,*) J_hole(homoAsub,homoBsub)
end if
CALL REMOVESP(pval)
printval='            <J>' // trim(pval) // '</J>'
write(16,*) trim(printval)

if ( QCP .eq. "T" ) then
   write(pval,*) dimen(homoAind(1))
else if ( QCP .eq. "G" ) then   
   write(pval,*) molAener(homoAsub)
else if (QCP.eq."GW")  then
   write(pval,*) QPE_molA(homoAsub)
end if

CALL REMOVESP(pval)
printval='            <e_A>' // trim(pval) // '</e_A>'
write(16,*) trim(printval)
if ( QCP .eq. "T" ) then
   write(pval,*) dimen(homoBind(1))
else if (QCP .eq. "G") then
   write(pval,*) molBener(homoBsub)
else if (QCP .eq. "GW") then
   write(pval,*) QPE_molB(homoBsub)
end if
CALL REMOVESP(pval)
printval='            <e_B>' // trim(pval) // '</e_B>'
write(16,*) trim(printval)
write(16,*) '        </channel>'
write(16,*) '        <channel name="multi">'
write(16,*) '           <molecule name="A">'
do i = 1, homoAsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   if ( QCP .eq. "T" ) then
      write(pval,*) dimen(homoAind(i))
   else if (QCP .eq. "G") then
      write(pval,*) molAener(homoAsub+1-i)
   elseif (QCP.eq."GW") then
      write(pval,*) QPE_molA(homoAsub+1-i)
   end if
   CALL REMOVESP(pval)
   printval='               <e_HOMOm'//trim(ival)//'>'//trim(pval) // '</e_HOMOm'//trim(ival)//'>'
   write(16,*) trim(printval)
end do
write(16,*) '           </molecule>'
write(16,*) '           <molecule name="B">'
do i = 1, homoBsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   if ( QCP .eq. "T" ) then 
      write(pval,*) dimen(homoBind(i))
   else if (QCP.eq."G") then
      write(pval,*) molBener(homoBsub+1-i)
   else if (QCP.eq."GW") then
      write(pval,*) QPE_molB(homoBsub+1-i)
   end if
   CALL REMOVESP(pval)
   printval='               <e_HOMOm'//trim(ival)//'>'//trim(pval) // '</e_HOMOm'//trim(ival)//'>'
   write(16,*) trim(printval)
end do
write(16,*) '           </molecule>'

write(16,*) '               <dimer name="integrals">'
do i = 1, homoAsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   do j = 1, homoBsub   
      write(jval,*) j-1
      CALL REMOVESP(jval)
      if (QCP.eq."T") then
         write(pval,*) J_hole(i,j)
      else
         write(pval,*) J_hole(homoAsub+1-i,homoBsub+1-j)
      end if
      CALL REMOVESP(pval)
      printval='                    <T_'//trim(ival)//trim(jval)//'>'//trim(pval) // '</T_'//trim(ival)//trim(jval)//'>'
      write(16,*) trim(printval)
   end do
end do
write(pval,*) J_sq_degen_h
CALL REMOVESP(pval)
printval='                    <J_sq_degen>'//trim(pval) // '</J_sq_degen>'
write(16,*) trim(printval)
write(pval,*) J_sq_boltz_h
CALL REMOVESP(pval)
printval='                    <J_sq_boltz>'//trim(pval) // '</J_sq_boltz>'
write(16,*) trim(printval)
write(16,*) '               </dimer>'
write(16,*) '        </channel>'
write(16,*) '    </transport>'
!!$
!!$
write(16,*) '    <transport name="electron">'
write(16,*) '        <channel name="single">'

write(pval,*) J_elec(1,1)
CALL REMOVESP(pval)
printval='            <J>' // trim(pval) // '</J>'
write(16,*) trim(printval)
if (QCP .eq. "T" ) then
   write(pval,*) dimen(lumoAind(1))
else if (QCP.eq."G") then
   write(pval,*) molAener(homoAsub+1)
else if (QCP.eq."GW") then
   write(pval,*) QPE_molA(homoAsub+1)
end if
CALL REMOVESP(pval)
printval='            <e_A>' // trim(pval) // '</e_A>'
write(16,*) trim(printval)
if (QCP .eq."T" ) then
   write(pval,*) dimen(lumoBind(1))
else if (QCP.eq."G") then
   write(pval,*) molBener(homoBsub+1)
else if (QCP.eq."GW") then
   write(pval,*) QPE_molB(homoBsub+1)
end if
CALL REMOVESP(pval)
printval='            <e_B>' // trim(pval) // '</e_B>'
write(16,*) trim(printval)
write(16,*) '        </channel>'
write(16,*) '        <channel name="multi">'
write(16,*) '               <molecule name="A">'

do i = 1, lumoAsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   if (QCP .eq. "T") then
      write(pval,*) dimen(lumoAind(i))
   else if (QCP.eq."G")then
      write(pval,*) molAener(homoAsub+i)
   elseif (QCP.eq."GW") then
      write(pval,*) QPE_molA(homoAsub+i)
   end if
   CALL REMOVESP(pval)
   printval='                    <e_LUMOp'//trim(ival)//'>'//trim(pval)//'</e_LUMOp'//trim(ival)//'>'
   write(16,*) trim(printval)
end do

write(16,*) '               </molecule>'
write(16,*) '               <molecule name="B">'

do i = 1, lumoBsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   if (QCP .eq. "T" ) then
      write(pval,*) dimen(lumoBind(i))
   else if (QCP.eq."G") then
      write(pval,*) molBener(homoBsub+i)
   else if (QCP.eq."GW") then
      write(pval,*) QPE_molB(homoBsub+i)
   end if
   CALL REMOVESP(pval)
   printval='                    <e_LUMOp'//trim(ival)//'>'//trim(pval)//'</e_LUMOp'//trim(ival)//'>'
   write(16,*) trim(printval)
end do
write(16,*) '               </molecule>'
write(16,*) '               <dimer name="integrals">'
do i = 1, lumoAsub
   write(ival,*) i-1
   CALL REMOVESP(ival)
   do j = 1, lumoBsub   
      write(jval,*) j-1
      CALL REMOVESP(jval)
      write(pval,*) J_elec(i,j)
      CALL REMOVESP(pval)
      printval='                    <T_'//trim(ival)//trim(jval)//'>'//trim(pval) // '</T_'//trim(ival)//trim(jval)//'>'
      write(16,*) trim(printval)
   end do
end do      
write(pval,*) J_sq_degen_e
CALL REMOVESP(pval)
printval='                    <J_sq_degen>'//trim(pval) // '</J_sq_degen>'
write(16,*) trim(printval)
write(pval,*) J_sq_boltz_e
CALL REMOVESP(pval)
printval='                    <J_sq_boltz>'//trim(pval) // '</J_sq_boltz>'
write(16,*) trim(printval)
write(16,*) '               </dimer>'
write(16,*) '        </channel>'
write(16,*) '    </transport>'
write(16,*) '</pair>'
!!$
!!$	
!!$
end subroutine XMLoutput


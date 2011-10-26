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
module precision

! Real kinds

integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

! Integer kinds

integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

!Complex kinds

integer, parameter :: kc4 = kr4                            ! single precision complex
integer, parameter :: kc8 = kr8                            ! double precision complex

end module precision

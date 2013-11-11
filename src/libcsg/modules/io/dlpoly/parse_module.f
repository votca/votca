      module parse_module

c***********************************************************************
c     
c     dl_poly module for defining parsing arrays
c
c copyright (c) 2004, daresbury laboratory
c     author    - w. smith    jan 2004
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
c     
c***********************************************************************

      integer, parameter :: lenrec=255
      character*1 record(lenrec)
      save record

      contains

      subroutine getrec(safe,idnode,ifile)

c*********************************************************************
c     
c     dl_poly subroutine to read a character string on one node
c     and broadcast it to all other nodes
c     
c     copyright daresbury laboratory 1994
c     author w.smith december 1994
c     
c*********************************************************************

      implicit none
      
      logical safe

      character*(lenrec) line
      integer export,import,idnode,ifile,i
      dimension export(lenrec),import(lenrec)
      
      safe=.true.
      
      call gsync()
      
      if(idnode.eq.0)then
        
        read(ifile,'(a150)',end=100)line

        do i=1,lenrec

          record(i)=line(i:i)
          export(i)=ichar(line(i:i))
          
        enddo
        
        call gstate(safe)
        call gisum(export,lenrec,import)
        
        return
        
  100   safe=.false.
        
        call gstate(safe)
        
      else
        
        call gstate(safe)
        if(.not.safe)return

        do i=1,lenrec

          export(i)=0

        enddo

        call gisum(export,lenrec,import)
        
        do i=1,lenrec
          
          record(i)=char(export(i))
          
        enddo
        
        return
        
      endif
      
      end subroutine getrec

      integer function intstr(word,len,lst)

c***********************************************************************
c     
c     dl_poly function for extracting integers from a 
c     character string
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     integer string
c     
c***********************************************************************
      
      implicit none

      logical flag,count,final
      character*1 n,word,ksn
      integer lst,len,j,isn

      dimension n(0:9),word(len)
      data n/'0','1','2','3','4','5','6','7','8','9'/

      isn=1
      lst=0
      ksn='+'
      intstr=0
      flag=.false.
      final=.false.
      count=.false.
      
      do while(lst.lt.len.and.(.not.final))

        lst=lst+1
        flag=.false.

        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            intstr=10*intstr+j
            count=.true.
            flag=.true.
            
          endif
          
        enddo

        if(count.and.(.not.flag))final=.true.
        if(flag.and.ksn.eq.'-')isn=-1
        ksn=word(lst)

      enddo

      intstr=isn*intstr

      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      return
      end function intstr

      real(8) function dblstr(word,len,lst)

c***********************************************************************
c     
c     dl_poly function for extracting double precisions from a 
c     character string. 
c     modified from dl_poly function intstr
c     
c     copyright - daresbury laboratory 1994
c     author    - w. smith may 1994.
c     modified  - t. forester april 1994
c     
c     parameters:
c     word   - input character string
c     len    - working length of character string
c     lst    - location of space character at end of
c     double precision string
c     
c***********************************************************************
      
      implicit none
      
      character*1 n,word,ksn,dot,d,e
      logical flag,ldot,start,final
      integer len,lst,iexp,idum,i,j,fail
      real(8) sn,ten,one
      dimension n(0:9),word(len)
      character*1, allocatable :: work(:)

      data n/'0','1','2','3','4','5','6','7','8','9'/
      data dot/'.'/
      data d/'d'/
      data e/'e'/
      
      allocate(work(len),stat=fail)

      lst=0
      sn=1.d0
      ksn='+'
      ten=10.d0
      one=1.d0
      
      dblstr=0.d0
      iexp=0
      idum=0
      start=.false.
      ldot=.false.
      final=.false.

      do while(lst.lt.len.and.(.not.final))
        
        lst=lst+1
        flag=.false.
        
        do j=0,9
          
          if(n(j).eq.word(lst))then
            
            dblstr=ten*dblstr+one*dble(j)
            flag=.true.
            start=.true.
            
          endif
          
        enddo
        
        if(dot.eq.word(lst))then
          
          flag=.true.
          ten=1.d0
          ldot=.true.
          start=.true.
          
        endif

        if(flag.and.ksn.eq.'-') sn=-1.d0
        if(ldot) one=one/10.d0
        ksn=word(lst)
        if(ksn.eq."D")ksn="d"
        if(ksn.eq."E")ksn="e"
        
        if(start)then
          
          if(d.eq.ksn.or.e.eq.ksn)then
            
            do i=1,len-lst
              work(i)=word(i+lst)
            enddo
            iexp=intstr(work,len-lst,idum)
            final=.true.

          endif

          if(.not.flag)final=.true.
          
        endif
        
      enddo
      
      dblstr=sn*dblstr*(10.d0**iexp)
      lst=lst+idum
      
      do j=lst,len
        word(j-lst+1)=word(j)
      enddo
      do j=len-lst+2,len
        word(j)=' '
      enddo

      deallocate(work,stat=idum)

      return
      end function dblstr

      subroutine strip(string,imax)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c     
c***********************************************************************

      implicit none

      integer i,imax,j
      character*1 string(imax)
      
      do i=1,imax

        if(string(1).eq.' ')then

          do j=1,imax-1

            string(j)=string(j+1)

          enddo

          string(imax)=' '

        endif

      enddo

      return
      end subroutine strip

      subroutine lowcase(string,length)

c***********************************************************************
c     
c     DL_POLY routine to lowercase a string of up to 255 characters.
c     Transportable to non-ASCII machines
c     
c     copyright daresbury laboratory 1993
c     author    t. forester     july 1993
c     
c***********************************************************************

      implicit none

      character*1 string(*)
      character*1 letter
      integer i,length

      do i=1,min(255,length)

        letter=string(i)

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

        string(i)=letter

      enddo

      return
      end subroutine lowcase

      subroutine copystring(oldstr,newstr,length)

c***********************************************************************
c     
c     DL_POLY routine to copy one string into another
c     
c     copyright daresbury laboratory
c     author    w. smith    jan 2004
c     
c***********************************************************************

      implicit none

      character*1 newstr(*),oldstr(*)
      integer i,length

      do i=1,length

        newstr(i)=oldstr(i)

      enddo

      return
      end subroutine copystring

      logical function findstring(seek,string,here)

c***********************************************************************
c     
c     DL_POLY routine to find an explicit string in an input record
c     note: variable `seek' is a character string while variable
c    `string' is a character*1 array i.e. code is application specific
c
c     copyright daresbury laboratory
c     author    w.smith   jan   2004
c     
c***********************************************************************

      implicit none

      integer i,n,m,here
      character*(*) seek
      character*1 string(lenrec)

      m=lenrec
      n=len(seek)
      findstring=.false.

      here=0
      do while(here.le.m-n.and.(.not.findstring))

        findstring=.true.

        do i=1,n
          if(seek(i:i).ne.string(here+i))findstring=.false.
        enddo

        here=here+1

      enddo

      return
      end function findstring

      subroutine striptext(string,length,nwords)

c***********************************************************************
c     
c     DL_POLY routine to strip leading text from a data record
c     
c     copyright daresbury laboratory
c     author   w.smith jan 2004
c     
c***********************************************************************

      implicit none

      logical final
      integer length,nwords,i,j,k
      character*1 string(length)
      
      do k=1,nwords

        i=0
        final=.false.
        
        do while(.not.final.and.i.lt.length)
          
          i=i+1
          
          if(string(1).eq.' ')then
            
            final=.true.
            
          else
            
            do j=1,length-1
              
              string(j)=string(j+1)
              
            enddo
            
            string(length)=' '
            
          endif
          
        enddo
        
      enddo

      return
      end subroutine striptext

      subroutine getword(word,string,len1,len2)

c***********************************************************************
c     
c     DL_POLY routine to fetch an 8 character word from a string
c     while ignoring leading blanks
c
c     copyright daresbury laboratory
c     author   w.smith jan 2004
c     
c***********************************************************************

      implicit none

      logical final
      character*8 word
      integer len1,len2,i,j,k
      character*1 wrdseq(len1),string(len2)
      
      do i=1,len1
        wrdseq(i)=' '
      enddo

      i=0
      k=0
      final=.false.
      
      do while(.not.final.and.i.lt.len2)
        
        i=i+1
        
        if(string(1).eq.' ')then
          
          if(k.gt.0)final=.true.
          
        else
          
          k=k+1
          wrdseq(k)=string(1)
          if(k.eq.len1)final=.true.

        endif
        
        do j=1,len2-1
          
          string(j)=string(j+1)
          
        enddo
        
        string(len2)=' '
          
      enddo
      
      word=mkwd8(wrdseq)

      return
      end subroutine getword

      character*8 function mkwd8(string)

c***********************************************************************
c     
c     DL_POLY routine to make an 8 character word from a string
c
c     copyright daresbury laboratory
c     author   w.smith nov 2006
c     
c***********************************************************************

      implicit none

      integer i
      character*1 string(*)
      
      do i=1,8
         mkwd8(i:i)=string(i)
      enddo
      
      return
      end function mkwd8
      
      end module parse_module



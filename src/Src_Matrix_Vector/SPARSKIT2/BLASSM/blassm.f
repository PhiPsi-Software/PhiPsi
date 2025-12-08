       subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                  c,jc,ic,nzmax,iw,ierr) 
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
      real*8 scal 
      logical values
      values = (job .ne. 0) 
      len = 0
      ic(1) = 1 
      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
      do 500 ii=1, nrow 
         do 200 ka=ia(ii), ia(ii+1)-1 
            if (values) scal = a(ka)
            jj   = ja(ka)
          do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  endif
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               endif
 100          continue
 200     continue
         do 201 k=ic(ii), len
          iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
      end
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
      logical values
      values = (job .ne. 0) 
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
      do 500 ii=1, nrow
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            if (values) c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
         iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
      end
      subroutine aplb1(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
      logical values
      values = (job .ne. 0) 
      ierr = 0
      kc = 1
      ic(1) = kc 
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1 
 5       continue 
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncol+1
         endif
         if (kb .le. kbmax) then 
            j2 = jb(kb)         
         else 
            j2 = ncol+1
         endif
         if (kc .gt. nzmax) goto 999
         if (j1 .eq. j2) then 
            if (values) c(kc) = a(ka)+b(kb)
            jc(kc) = j1
            ka = ka+1
            kb = kb+1
            kc = kc+1
         else if (j1 .lt. j2) then
            jc(kc) = j1
            if (values) c(kc) = a(ka)
            ka = ka+1
            kc = kc+1
         else if (j1 .gt. j2) then
            jc(kc) = j2
            if (values) c(kc) = b(kb)
            kb = kb+1
            kc = kc+1
         endif
         if (ka .le. kamax .or. kb .le. kbmax) goto 5
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i 
      return
      end
      subroutine aplsb (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,
     *     iw,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
      ierr = 0
      len = 0
      ic(1) = 1 
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
      do 500 ii=1, nrow
         do 200 ka=ia(ii), ia(ii+1)-1 
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol 
            c(len)  = a(ka) 
            iw(jcol)= len
 200     continue
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               c(len)  = s*b(kb)
               iw(jcol)= len
            else
               c(jpos) = c(jpos) + s*b(kb)
            endif
 300     continue
         do 301 k=ic(ii), len
         iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
      end
      subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
      ierr = 0
      kc = 1
      ic(1) = kc 
      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1 
 5       continue 
         if (ka .le. kamax .or. kb .le. kbmax) then 
            if (ka .le. kamax) then
               j1 = ja(ka)
            else
               j1 = ncol+1
            endif
            if (kb .le. kbmax) then 
               j2 = jb(kb)         
            else 
               j2 = ncol+1
            endif
            if (j1 .eq. j2) then 
               c(kc) = a(ka)+s*b(kb)
               jc(kc) = j1
               ka = ka+1
               kb = kb+1
               kc = kc+1
            else if (j1 .lt. j2) then
               jc(kc) = j1
               c(kc) = a(ka)
               ka = ka+1
               kc = kc+1
            else if (j1 .gt. j2) then
               jc(kc) = j2
               c(kc) = s*b(kb)
               kb = kb+1
               kc = kc+1
            endif
            if (kc .gt. nzmax) goto 999
            goto 5
         endif
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i 
      return
      end
      subroutine apmbt (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*) 
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*) 
      logical values
      values = (job .ne. 0) 
      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      endif
      ljob = 0
      if (values) ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
      if (job .eq. -1) then
         do 2 k=1,len
         c(k) = -c(k)
 2       continue
      endif
      do 500 ii=1, nrow
         do 200 k = ic(ii),ic(ii+1)-1
            iw(jc(k)) = k
 200     continue
         do 300 ka = ia(ii), ia(ii+1)-1 
            jcol = ja(ka)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               
               ic(len) = ii
               if (values) c(len)  = a(ka)
            else
               if (values) c(jpos) = c(jpos) + a(ka)
            endif
 300     continue
         do 301 k=ic(ii), ic(ii+1)-1
            iw(jc(k)) = 0
 301     continue
 500  continue
      ljob = 2
      if (values) ljob = 3
      do 501 i=1, nrow+1
         iw(i) = ic(i) 
 501  continue
      call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
      ljob = 0
      if (values) ljob = 1
      call coicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
      subroutine aplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)
      real*8 a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),iw(*)
      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue
      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      endif
      ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic) 
      do 2 k=1,len
 2       c(k) = c(k)*s
         do 500 ii=1, nrow
            do 200 k = ic(ii),ic(ii+1)-1
               iw(jc(k)) = k
 200        continue
            do 300 ka = ia(ii), ia(ii+1)-1 
               jcol = ja(ka)
               jpos = iw(jcol)
           if (jpos .eq. 0) then
              len = len+1
              if (len .gt. nzmax) goto 999
              jc(len) = jcol              
              ic(len) = ii
              c(len)  = a(ka)
           else
              c(jpos) = c(jpos) + a(ka)
           endif
 300    continue
        do 301 k=ic(ii), ic(ii+1)-1
           iw(jc(k)) = 0
 301    continue
 500  continue
      ljob = 3
      do 501 i=1, nrow+1
         iw(i) = ic(i) 
 501  continue
      call csrcoo (nrow,ljob,nnzb,c,jc,iw,nnzb,c,ic,jc,ierr)
      ljob = 1
      call coicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii) 
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue
      if (job .eq. 0) return
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
      end 
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1) 
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k)) 
 2       continue
 1    continue
      if (job .eq. 0) return
      do 3 ii=1, nrow+1
         ib(ii) = ia(ii)
 3    continue
      do 31 k=ia(1), ia(nrow+1) -1 
         jb(k) = ja(k)
 31   continue
      return
      end 
      subroutine aplsca (nrow, a, ja, ia, scal,iw) 
      real*8 a(*), scal
      integer ja(*), ia(nrow+1),iw(*)
      logical test
      call diapos (nrow,ja,ia,iw)
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            a(iw(j)) = a(iw(j)) + scal 
         endif
 1    continue
      if (icount .eq. 0) return
      ko = ia(nrow+1)+icount
      do 5 ii=nrow, 1, -1 
         k1 = ia(ii)
         k2 = ia(ii+1)-1 
         ia(ii+1) = ko
         test = (iw(ii) .eq. 0) 
         do 4 k = k2,k1,-1 
            j = ja(k)
            if (test .and. (j .lt. ii)) then 
               test = .false. 
               ko = ko - 1
               a(ko) = scal 
               ja(ko) = ii
               iw(ii) = ko
            endif
            ko = ko-1
            a(ko) = a(k) 
            ja(ko) = j
 4       continue
         if (test) then
            ko = ko-1
            a(ko) = scal 
            ja(ko) = ii
            iw(ii) = ko
         endif
 5    continue
      ia(1) = ko 
      return
      end
      subroutine apldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw) 
      real*8 a(*), b(*), diag(nrow) 
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*)
      logical test
      if (job .ne. 0) then 
         nnz = ia(nrow+1)-1
         do 2  k=1, nnz
            jb(k) = ja(k)
            b(k)  = a(k) 
 2       continue
         do 3 k=1, nrow+1
            ib(k) = ia(k)
 3       continue
      endif 
      call diapos (nrow,ja,ia,iw)
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            b(iw(j)) = a(iw(j)) + diag(j) 
         endif
 1    continue
      if (icount .eq. 0) return
      ko = ib(nrow+1)+icount
      do 5 ii=nrow, 1, -1 
         k1 = ib(ii)
         k2 = ib(ii+1)-1 
         ib(ii+1) = ko
         test = (iw(ii) .eq. 0) 
         do 4 k = k2,k1,-1 
            j = jb(k)
            if (test .and. (j .lt. ii)) then 
               test = .false. 
               ko = ko - 1
               b(ko) = diag(ii) 
               jb(ko) = ii
               iw(ii) = ko
            endif
            ko = ko-1
            b(ko) = a(k) 
            jb(ko) = j
 4       continue
         if (test) then
            ko = ko-1
            b(ko) =  diag(ii) 
            jb(ko) = ii
            iw(ii) = ko
         endif
 5    continue
      ib(1) = ko 
      return
      end

       

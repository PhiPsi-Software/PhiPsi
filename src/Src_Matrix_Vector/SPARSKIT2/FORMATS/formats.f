      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr) 
      real*8 dns(ndns,*),a(*)
      integer ja(*),ia(*)
      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
            dns(i,j) = 0.0d0
 2       continue
 1    continue
      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k) 
            if (j .gt. ncol) then
               ierr = i
               return
            endif
            dns(i,j) = a(k)
 3       continue           
 4    continue
      return
      end
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
      real*8 dns(ndns,*),a(*)
      integer ia(*),ja(*)
      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1, ncol 
            if (dns(i,j) .eq. 0.0d0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            endif
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3       continue           
         ia(i+1) = next
 4    continue
      return
      end
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
      subroutine coicsr (n,nnz,job,a,ja,ia,iwk)
      integer ia(*),ja(nnz),iwk(n+1) 
      real*8 a(*)
      real*8 t,tnext 
      logical values
      t = 0.0
      tnext = 0.0
      values = (job .eq. 1) 
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue 
      iwk(1) = 1 
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue 
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
 6    k = k+1                 
      ipos = iwk(i)
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
      if (values) a(ipos)  = t
      ja(ipos) = j
      iwk(i) = ipos+1
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext 
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
      goto 5
 70   do 80 i=1,n 
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
      end
      subroutine csrcoo(nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
      real*8 a(*),ao(*) 
      integer ir(*),jc(*),ja(*),ia(nrow+1) 
      integer,intent(in):: nnz
      integer calculated_nnz
      ierr = 0
      calculated_nnz = ia(nrow+1)-1
      if (calculated_nnz .gt. nzmax) then
         ierr = 1
         return
      endif
      goto (3,2,1) job
 1    do 10 k=1,calculated_nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,calculated_nnz
         jc(k) = ja(k)
 11   continue
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
      end
      subroutine csrssr (nrow,a,ja,ia,nzmax,ao,jao,iao,ierr)
      real*8 a(*), ao(*), t
      integer ia(*), ja(*), iao(*), jao(*)
      ierr = 0
      ko = 0
      do  7 i=1, nrow
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            if (ko .gt. nzmax) then
               ierr = i
               return
            endif
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t
         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
      iao(nrow+1) = ko+1
      return
      end
      subroutine ssrcsr(job, value2, nrow, a, ja, ia, nzmax,
     &                  ao, jao, iao, indu, iwk, ierr)
      integer            ierr, job, nrow, nzmax, value2
      integer            ia(nrow+1), iao(nrow+1), indu(nrow),
     &                   iwk(nrow+1), ja(*), jao(nzmax)
      real*8             a(*), ao(nzmax)
      integer            i, ipos, j, k, kfirst, klast, ko, kosav, nnz
      real*8             tmp
      ierr = 0
      do 10 i = 1, nrow
         indu(i) = 0
         iwk(i) = 0
 10   continue
      iwk(nrow+1) = 0
      do 30 i = 1, nrow
         do 20 k = ia(i), ia(i+1) - 1
            j = ja(k)
            if (j.ne.i)
     &         iwk(j+1) = iwk(j+1) + 1
 20      continue
 30   continue
      iwk(1) = 1
      do 40 i = 1, nrow
         indu(i) = iwk(i) + ia(i+1) - ia(i)
         iwk(i+1) = iwk(i+1) + indu(i)
         indu(i) = indu(i) - 1
 40   continue
      nnz = iwk(nrow+1) - 1
      if (nnz.gt.nzmax) then
         ierr = nnz
         return
      endif
      kosav = iwk(nrow+1)
      do 60 i = nrow, 1, -1
         klast = ia(i+1) - 1
         kfirst = ia(i)
         iao(i+1) = kosav
         kosav = iwk(i)
         ko = iwk(i) - kfirst
         iwk(i) = ko + klast + 1
         do 50 k = klast, kfirst, -1
            if (value2.ne.0)
     &         ao(k+ko) = a(k)
            jao(k+ko) = ja(k)
 50      continue
 60   continue
      iao(1) = 1
      do 80 i = 1, nrow
         do 70 k = iao(i), indu(i)
            j = jao(k)
            if (j.ne.i) then
               ipos = iwk(j)
               if (value2.ne.0)
     &            ao(ipos) = ao(k)
               jao(ipos) = i
               iwk(j) = ipos + 1
            endif
 70      continue
 80   continue
      if (job.le.0) return
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = iao(i)
 90   continue
      iwk(nrow+1) = iao(nrow+1)
      k = 1
      do 120 i = 1, nrow
         iao(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = jao(ipos)
            if (indu(j).eq.0) then
               if (value2.ne.0) then
                  if (ao(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     jao(k) = jao(ipos)
                     ao(k) = ao(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  jao(k) = jao(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
               ao(indu(j)) = ao(indu(j)) + ao(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
         do 110 ipos = iao(i), k - 1
            indu(jao(ipos)) = 0
 110     continue
 120  continue
      iao(nrow+1) = k
      if (job.le.1) return
      do 140 i = 1, nrow
         klast = iao(i+1) - 1
         kfirst = iao(i)
 130     if (klast.gt.kfirst) then
            if (jao(klast).lt.i .and. jao(kfirst).ge.i) then
               j = jao(klast)
               jao(klast) = jao(kfirst)
               jao(kfirst) = j
               if (value2.ne.0) then
                  tmp = ao(klast)
                  ao(klast) = ao(kfirst)
                  ao(kfirst) = tmp
               endif
            endif
            if (jao(klast).ge.i)
     &         klast = klast - 1
            if (jao(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
         if (jao(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
      do 190 i = 1, nrow
         do 160 ipos = iao(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), iao(i+1)-1
            do 170 j = iao(i+1)-1, ipos+1, -1
               k = j - 1
               if (jao(k).gt.jao(j)) then
                  ko = jao(k)
                  jao(k) = jao(j)
                  jao(j) = ko
                  if (value2.ne.0) then
                     tmp = ao(k)
                     ao(k) = ao(j)
                     ao(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
      return
      end
      subroutine xssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      real*8 a(*),ao(nzmax)
      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0     
 1    continue
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue 
 3    continue
      indu(1) = 1 
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
      nnz = indu(nrow+1)-1 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i) 
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
            ko = ko+1
 5       continue
         indu(i) = ko 
 6    continue
      iao(1) = 1
      do 8 i=1,nrow
         do 9 k=iao(i), iao(i+1)-1 
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1 
 9       continue
 8    continue
      return
      end
      subroutine csrell (nrow,a,ja,ia,maxcol,coef,jcoef,ncoef,
     *                   ndiag,ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)  
      real*8 a(*), coef(ncoef,1)
      ierr = 0
      ndiag = 0
      do 3 i=1, nrow
         k = ia(i+1)-ia(i)
         ndiag = max0(ndiag,k) 
 3    continue
      if (ndiag .gt. maxcol) then
         ierr = 1 
         return
      endif
      do 4 j=1,ndiag 
         do 41 i=1,nrow
            coef(i,j) = 0.0d0
            jcoef(i,j) = i
 41      continue
 4    continue
      do 6 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1)-1
         do 5 k=k1,k2
            coef(i,k-k1+1) = a(k)
            jcoef(i,k-k1+1) = ja(k)
 5       continue
 6    continue
      return
      end
      subroutine ellcsr(nrow,coef,jcoef,ncoef,ndiag,a,ja,ia,nzmax,ierr)
      integer ia(nrow+1), ja(*), jcoef(ncoef,1) 
      real*8 a(*), coef(ncoef,1)
      ierr = 0
      kpos = 1
      ia(1) = kpos
      do 6 i=1, nrow
         do 5 k=1,ndiag
            if (coef(i,k) .ne. 0.0d0) then
               if (kpos .gt. nzmax) then
                  ierr = kpos
                  return
               endif
               a(kpos) = coef(i,k)
               ja(kpos) = jcoef(i,k)
               kpos = kpos+1
            endif
 5       continue
         ia(i+1) = kpos
 6    continue        
      return
      end
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)
      real*8 a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)
      icount = 0
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1 
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
      iptr = n + ia(n+1) - icount
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr-1
            endif
 100     continue
 500  continue
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i) 
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return        
      end
      subroutine msrcsr (n,a,ja,ao,jao,iao,wk,iwk)
      real*8 a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1),iwk(n+1)
      logical added
      do 1 i=1,n
         wk(i) = a(i)
         iwk(i) = ja(i)
 1    continue
      iwk(n+1) = ja(n+1)
      iao(1) = 1
      iptr = 1
      do 500 ii=1,n 
         added = .false.
         idiag = iptr + (iwk(ii+1)-iwk(ii)) 
         do 100 k=iwk(ii),iwk(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            else 
               idiag = iptr
               iptr = iptr+1
               added = .true.
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr+1
            endif
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr 
 500  continue
      return    
      end
      subroutine csrcsc(n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
      call csrcsc2 (n,n,job,ipos,a,ja,ia,ao,jao,iao)
      end
      subroutine csrcsc2 (n,n2,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n2+1),ja(*),jao(*)
      real*8  a(*),ao(*)
      do 1 i=1,n2+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
      iao(1) = ipos 
      do 4 i=1,n2
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
      do 7 i=n2,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
      end
      subroutine csrlnk (n,ia,link) 
      integer n,  ia(n+1), link(*)
      integer i, k
      do 100 i =1, n
         istart = ia(i) 
         iend = ia(i+1)-1
         if (iend .gt. istart) then
            do 99  k=istart, iend-1 
               link(k) = k+1
 99         continue
            link(iend) = 0
         else
            ia(i) = 0
         endif
 100  continue
      return
      end
      subroutine lnkcsr (n, a, jcol, istart, link, ao, jao, iao) 
      real*8 a(*), ao(*) 
      integer n, jcol(*), istart(n), link(*), jao(*), iao(*) 
      integer irow, ipos, next
      ipos = 1
      iao(1) = ipos
      do 100 irow =1, n
         next = istart(irow)
 10      if (next .eq. 0) goto 99
         jao(ipos) = jcol(next)
         ao(ipos)  = a(next)
         ipos = ipos+1
         next = link(next) 
         goto 10
 99      iao(irow+1) = ipos 
 100  continue
      return
      end
      subroutine csrdia (n,idiag,job,a,ja,ia,ndiag,
     *                   diag,ioff,ao,jao,iao,ind)
      real*8 diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)
      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia(n,ja,ia,ind,idum)
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      endif
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0d0
 54      continue
 55   continue
      ko = 1
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1 
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
      end
      subroutine diacsr (n,job,idiag,diag,ndiag,ioff,a,ja,ia)
      real*8 diag(ndiag,idiag), a(*), t
      integer ia(*), ja(*), ioff(*)
      ia(1) = 1
      ko = 1
      do 80 i=1, n
         do 70 jj = 1, idiag
            j = i+ioff(jj) 
            if (j .lt. 1 .or. j .gt. n) goto 70
            t = diag(i,jj) 
            if (job .eq. 0 .and. t .eq. 0.0d0) goto 70
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 70      continue
         ia(i+1) = ko
 80   continue
      return
      end
      subroutine bsrcsr (job, n, m, na, a, ja, ia, ao, jao, iao)
      implicit none
      integer job, n, m, na, ia(*), ja(*), jao(*), iao(n+1)
      real*8 a(na,*), ao(*)
      integer i, i1, i2, ij, ii, irow, j, jstart, k, krow
      logical val
      val = (job.ne.0)
      irow = 1        
      krow = 1
      iao(irow) = 1      
      do 2 ii=1, n
         i1 = ia(ii)
         i2 = ia(ii+1)-1
         do 23 i=1,m
            do 21 k=i1, i2
               jstart = m*(ja(k)-1)
               do 22  j=1,m
                  ij = (j-1)*m + i
                  if (val) ao(krow) = a(ij,k) 
                  jao(krow) = jstart+j
                  krow = krow+1
 22            continue            
 21         continue
            irow = irow+1 
            iao(irow) = krow
 23      continue
 2    continue
      return
      end
      subroutine csrbsr (job,nrow,m,na,a,ja,ia,ao,jao,iao,iw,ierr)
      implicit none
      integer job,ierr,nrow,m,na,ia(nrow+1),ja(*),jao(na),iao(*),iw(*)
      real*8 a(*),ao(na,*)
      integer nr, m2, io, ko, ii, len, k, jpos, j, i, ij, jr, irow   
      logical vals  
      ierr = 0 
      if (m*m .gt. na) ierr = 2 
      if (m .eq. 0) ierr = 1 
      if (ierr .ne. 0) return
      vals = (job .gt. 0) 
      nr = 1 + (nrow-1) / m
      m2 = m*m 
      ko = 1 
      io = 1 
      iao(io) = 1 
      len = 0 
         do j=1, nr
            iw(j) = 0
         enddo
      do ii=1, nrow, m
         irow = 0
         do while (ii+irow .le. nrow .and. irow .le. m-1) 
            do k=ia(ii+irow),ia(ii+irow+1)-1
               j = ja(k)-1 
               jr = j/m + 1                
               j = j - (jr-1)*m 
               jpos = iw(jr) 
               if (jpos .eq. 0) then
                  iw(jr) = ko 
                  jao(ko) = jr 
                  if (vals) then
                     do i=1, m2
                        ao(i,ko) = 0.0d0
                     enddo
                     ij = j*m + irow + 1 
                     ao(ij,ko) = a(k) 
                  endif
                  ko = ko+1
               else
                  jao(jpos) = jr 
                  ij = j*m + irow + 1 
                  if (vals) ao(ij,jpos) = a(k) 
               endif 
            enddo  
            irow = irow+1
         enddo
         do j = iao(io),ko-1 
            iw(jao(j)) = 0
         enddo
         if (job .eq. -1) then
            len = len + ko-1
            ko = 1
         else
            io = io+1 
            iao(io) = ko
         endif
      enddo
      if (job .eq. -1) iao(1) = len
      return
      end
      subroutine csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)
      real*8 a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
      ierr = 0
      if (job .eq. 1) call getbwd(n,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return
      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
            abd(ii,j) = 0.0d0
 10      continue
 15   continue
      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k) 
 20      continue
 30   continue
      return
      end
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)
      real*8 a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)
      ierr = 0
      if (lowd .gt. nabd .or. lowd .le. 0) then 
         ierr = -1
         return
      endif
      ko = 1
      ia(1) = 1
      do 30 irow=1,n
         i = lowd 
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j) 
             if (t .eq. 0.0d0) goto 19
             if (ko .gt. len) then 
               ierr = irow 
               return
            endif
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
 21      ia(irow+1) = ko
 30   continue
      return
      end
      subroutine csrssk (n,imod,a,ja,ia,asky,isky,nzmax,ierr)
      real*8 a(*),asky(nzmax) 
      integer n, imod, nzmax, ierr, ia(n+1), isky(n+1), ja(*)
      ierr = 0
      isky(1) = 0
      do 3 i=1,n
         ml = 0
         do 31 k=ia(i),ia(i+1)-1 
            ml = max(ml,i-ja(k)+1) 
 31      continue
         isky(i+1) = isky(i)+ml
 3    continue
      nnz = isky(n+1) 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
      do 1 k=1, nnz 
         asky(k) = 0.0d0
 1    continue
      do 4 i=1,n
         kend = isky(i+1) 
         do 41 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .le. i) asky(kend+j-i) = a(k)
 41      continue
 4    continue
      if (imod .eq. 0) return
      if (imod .eq. 1) then 
         do 50 k=1, n+1
            isky(k) = isky(k)+1
 50      continue
      endif
      if (imod .eq. 2) then
         do 60 k=1, n
            isky(k) = isky(k+1) 
 60      continue
      endif
      return
      end
      subroutine sskssr (n,imod,asky,isky,ao,jao,iao,nzmax,ierr)
      real*8 asky(*),ao(nzmax) 
      integer n, imod,nzmax,ierr, isky(n+1),iao(n+1),jao(nzmax) 
      integer next, kend, kstart, i, j 
      ierr = 0
      if (imod.ne.0 .and. imod.ne.1 .and. imod .ne. 2) then
         ierr =-1
         return
      endif 
      next = 1
      kend = 0 
      if (imod .eq. 1) kend = isky(1)-1
      if (imod .eq. 0) kend = isky(1) 
      do 50 i=1,n
         iao(i) = next
         kstart = kend+1
         if (imod .eq. 0) kend = isky(i+1)
         if (imod .eq. 1) kend = isky(i+1)-1
         if (imod .eq. 2) kend = isky(i) 
         do 40 k=kstart,kend
            if (asky(k) .eq. 0.0d0) goto 40
            j = i-(kend-k) 
            jao(next) = j
            ao(next)  = asky(k)
            next=next+1
            if (next .gt. nzmax+1) then
               ierr = i
               return
            endif 
 40      continue
 50    continue
      iao(n+1) = next
      return
      end
      subroutine csrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow) 
      real*8 a(*), ao(*)
      idiag = 0
      ilo = nrow 
      do 10 j=1, nrow
         iperm(j) = j 
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len) 
         idiag = max(idiag,len) 
         jao(j) = len
 10   continue 
      call dcsort (jao, nrow, iao, iperm, ilo, idiag) 
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k)) 
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i) 
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
      end
      subroutine jadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao) 
      integer ja(*), jao(*), ia(idiag+1), iperm(nrow), iao(nrow+1) 
      real*8 a(*), ao(*)
      do 137 j=1,nrow
         jao(j) = 0
 137  continue
      do 140 i=1, idiag
         len = ia(i+1)-ia(i) 
         do 138 k=1,len
            jao(iperm(k)) = jao(iperm(k))+1
 138     continue
 140  continue
      kpos = 1
      iao(1) = 1
      do 141 i=1, nrow 
         kpos = kpos+jao(i) 
         iao(i+1) = kpos
 141  continue
      do 200 jj = 1, idiag
         k1 = ia(jj)-1
         len = ia(jj+1)-k1-1 
         do 160 k=1,len
            kpos = iao(iperm(k))
            ao(kpos) = a(k1+k) 
            jao(kpos) = ja(k1+k) 
            iao(iperm(k)) = kpos+1
 160     continue
 200  continue
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
      subroutine dcsort(ival, n, icnt, index, ilo, ihi)
      integer n, ilo, ihi, ival(n), icnt(ilo:ihi), index(n)
      integer i, j, ivalj
      do 10 i = ilo, ihi
        icnt(i) = 0
 10   continue
      do 20 i = 1, n
        icnt(ival(i)) = icnt(ival(i)) + 1
 20   continue
      do 30 i = ihi-1,ilo,-1
        icnt(i) = icnt(i) + icnt(i+1)
 30   continue
      do 40 j = n, 1, -1
        ivalj = ival(j)
        index(icnt(ivalj)) = j
        icnt(ivalj) = icnt(ivalj) - 1
 40   continue
      return
      end
      subroutine cooell(job,n,nnz,a,ja,ia,ao,jao,lda,ncmax,nc,ierr)
      implicit none
      integer job,n,nnz,lda,ncmax,nc,ierr
      integer ja(nnz),ia(nnz),jao(lda,ncmax)
      real*8  a(nnz),ao(lda,ncmax)
      integer i,j,k,ip
      real*8  zero
      logical copyval
      parameter (zero=0.0D0)
      copyval = (job.ne.0)
      if (lda .lt. n) then
         ierr = -1
         return
      endif
      do i = 1, n
         jao(i,ncmax) = 0
      enddo
      nc = 0
      do k = 1, nnz
         i = ia(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
      enddo
      nc = 0
      do i = 1, n
         if (nc.lt.jao(i,ncmax)) nc = jao(i,ncmax)
         jao(i,ncmax) = 0
      enddo
      if (nc.gt.ncmax) then
         ierr = nc
         return
      endif
      do k = 1, nnz
         i = ia(k)
         j = ja(k)
         jao(i,ncmax) = jao(i,ncmax) + 1
         ip = jao(i,ncmax)
         if (ip.gt.nc) nc = ip
         if (copyval) ao(i,ip) = a(k)
         jao(i,ip) = j
      enddo
      do i = 1, n
         do j = ia(i+1)-ia(i)+1, nc
            jao(i,j)=i
            if(copyval) ao(i,j) = zero
         enddo
      enddo
      ierr = 0
      return
      end
      subroutine xcooell(n,nnz,a,ja,ia,ac,jac,nac,ner,ncmax,ierr)
      real*8 a(nnz), ac(nac,ner)
      integer ja(nnz), ia(nnz), jac(nac,ner), ierr, ncmax, icount
      ierr = 0
      do 4 in = 1,ner
         do 4 innz =1,n
            jac(innz,in) = n
            ac(innz,in) = 0.0d0
 4    continue
      ncmax = 1
      do 10 is = 1,n
         k = 0
         do 30 ii = 1,nnz
            if(ia(ii).eq.is)then
               k = k + 1
               if (k .le. ner) then
                  ac(is,k) = a(ii)
                  jac(is,k) = ja(ii)
               endif 
            endif
 30      continue
         if (k.ge.ncmax) ncmax = k
 10   continue
      if (ncmax.eq.ner) ierr = 0
      if (ncmax.gt.ner) then
         ierr = -1
         return
      endif
      do 45 in = 1,ncmax
         icount = 0
         do 44 inn =1,n
            if (ac(inn,in).ne.0.0d0) icount = 1
 44      continue
         if (icount.eq.0) then
            ierr = 1
            return
         endif
 45   continue
      do 55 inn = 1,n
         icount = 0
         do 54 in =1,ncmax
            if (ac(inn,in).ne.0.0d0) icount = 1
 54      continue
         if (icount.eq.0) then
            ierr = 2
            return
         endif
 55   continue
      return
      end
      subroutine csruss (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer nrow,ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),
     *     iau(nrow+1)
      integer i, j, k, kl, ku 
      do 1 i=1,nrow+1
         iau(i) = 0
 1    continue
      do 3 i=1, nrow
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)
            if (j .gt. i) iau(j+1) = iau(j+1)+1
 2       continue 
 3    continue
      iau(1) = 1
      do 4 i=1,nrow
         iau(i+1) = iau(i)+iau(i+1)
         ial(i+1) = ial(i)+ial(i+1)
 4    continue
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
         do 71 k = ia(i), ia(i+1)-1
            j = ja(k) 
            if (j  .gt. i) then
               ku = iau(j) 
               au(ku) = a(k)
               jau(ku) = i
               iau(j) = ku+1
            elseif (j  .eq. i) then
               diag(i) = a(k) 
            elseif (j .lt. i) then
               al(kl) = a(k)
               jal(kl) = j
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
      do 8 i=nrow,1,-1
         iau(i+1) = iau(i)
 8    continue
      iau(1) = 1
      end 
      subroutine usscsr (nrow,a,ja,ia,diag,al,jal,ial,au,jau,iau) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1),jau(*),iau(nrow+1)
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
      do 3 i=1, nrow
         do 2 k=iau(i), iau(i+1)-1 
            j = jau(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ja(ka) = i
         ia(i) = ka+1
 6    continue
      do 8 i=1, nrow
         do 7 k=iau(i), iau(i+1)-1
            jak = jau(k) 
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
      end 
      subroutine csrsss (nrow,a,ja,ia,sorted,diag,al,jal,ial,au)
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1)
      logical sorted 
      kl = 1
      ial(1) = kl
      do  7 i=1, nrow
         do 71 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .eq. i) then
               diag(i) = a(k) 
            elseif (jak .lt. i) then
               al(kl) = a(k)
               jal(kl) = jak
               kl = kl+1
            endif
 71      continue
         ial(i+1) = kl 
 7    continue
      if (.not. sorted) then
         call csort (nrow, al, jal, ial, .true.) 
      endif
      do  8 i=1, nrow
         do 81 k = ia(i), ia(i+1)-1
            jak = ja(k) 
            if (jak  .gt. i) then
               ku = ial(jak) 
               au(ku) = a(k)
               ial(jak) = ku+1
            endif
 81      continue
 8    continue
      do 9 i=nrow,1,-1
         ial(i+1) = ial(i)
 9    continue
      ial(1) = 1
      end 
      subroutine ssscsr (nrow,a,ja,ia,diag,al,jal,ial,au) 
      real*8 a(*),al(*),diag(*),au(*) 
      integer ja(*),ia(nrow+1),jal(*),ial(nrow+1) 
      do 1 i=1, nrow
         ia(i+1) = ial(i+1)-ial(i)+1
 1    continue
      do 3 i=1, nrow
         do 2 k=ial(i), ial(i+1)-1 
            j = jal(k)
            ia(j+1) = ia(j+1)+1
 2       continue 
 3    continue
      ia(1) = 1
      do 4 i=1,nrow
         ia(i+1) = ia(i)+ia(i+1)
 4    continue
      do 6 i=1, nrow
         ka = ia(i) 
         do 5 k=ial(i), ial(i+1)-1
            a(ka) = al(k) 
            ja(ka) = jal(k) 
            ka = ka+1
 5       continue
         a(ka) = diag(i) 
         ia(i) = ka+1
 6    continue
      do 8 i=1, nrow
         do 7 k=ial(i), ial(i+1)-1
            jak = jal(k) 
            ka = ia(jak) 
            a(ka) = au(k) 
            ja(ka) = i
            ia(jak) = ka+1
 7       continue
 8    continue
      do 9 i=nrow,1,-1
         ia(i+1) = ia(i)
 9    continue
      ia(1) = 1
      end
      subroutine csrvbr(n,ia,ja,a,nr,nc,kvstr,kvstc,ib,jb,kb,
     &     b, job, iwk, nkmax, nzmax, ierr )
      integer n, ia(n+1), ja(*), nr, nc, ib(*), jb(nkmax-1), kb(nkmax)
      integer kvstr(*), kvstc(*), job, iwk(*), nkmax, nzmax, ierr
      real*8  a(*), b(nzmax)
      integer ncol, nb, neqr, numc, a0, b0, b1, k0, i, ii, j, jj, jnew
      logical sorted
      ierr = 0
      call csorted(n, ia, ja, sorted)
      if (.not. sorted) then
         call csort (n, a, ja, ia, .true.)
      endif
      if (job .eq. 1 .or. job .eq. 2) then
         ncol = 0
         do i = 2, n
            ncol = max0(ncol, ja(ia(i)-1))
         enddo
         do i = 1, ncol
            iwk(i) = 0
         enddo
         call csrkvstr(n, ia, ja, nr, kvstr)
         call csrkvstc(n, ia, ja, nc, kvstc, iwk)
      endif
      if (job .eq. 2) then
         if (kvstr(nr+1) .ne. kvstc(nc+1)) then
            ierr = 3
            return
         endif
         call kvstmerge(nr, kvstr, nc, kvstc, nb, iwk)
         nr = nb
         nc = nb
         do i = 1, nb+1
            kvstr(i) = iwk(i)
            kvstc(i) = iwk(i)
         enddo
      endif
      do i = 1, nc
         do j = kvstc(i), kvstc(i+1)-1
            iwk(j) = i
         enddo
      enddo
      ncol = kvstc(nc+1)-1
      if (job .eq. 0) goto 400
      a0 = 1
      b0 = 1
      k0 = 1
      kb(1) = 1
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
         j = 0
         do jj = ia(kvstr(i)), ia(kvstr(i)+1)-1
            jnew = iwk(ja(jj))
            if (jnew .ne. j) then
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
               j = jnew
               b0 = b0 + neqr * (kvstc(j+1) - kvstc(j))
               kb(k0+1) = b0
               jb(k0) = j
               k0 = k0 + 1
            endif
         enddo
         do ii = 0, neqr-1
            b1 = kb(ib(i))+ii
            do jj = 1, numc
               if (b1 .gt. nzmax) then
                  ierr = 2
                  write (*,*) 'csrvbr: no space in b for block row ', i
                  return
               endif
               b(b1) = a(a0)
               b1 = b1 + neqr
               a0 = a0 + 1
            enddo
         enddo
      enddo
      ib(nr+1) = k0
      return
 400  continue
      do i = ncol+1, ncol+nc
         iwk(i) = 0
      enddo
      k0 = 1
      kb(1) = 1
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ia(kvstr(i)+1) - ia(kvstr(i))
         ib(i) = k0
         do jj = ia(kvstr(i)), ia(kvstr(i+1))-1
            iwk(iwk(ja(jj))+ncol) = 1
         enddo
         do j = 1, nc
            if (iwk(j+ncol) .ne. 0) then
               if (k0+1 .gt. nkmax) then
                  ierr = 1
                  write (*,*) 'csrvbr: no space in kb for block row ', i
                  return
               endif
               kb(k0+1) = kb(k0) + neqr * (kvstc(j+1) - kvstc(j))
               jb(k0) = j
               k0 = k0 + 1
               iwk(j+ncol) = 0
            endif
         enddo
      enddo
      ib(nr+1) = k0
      a0 = 1
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         do ii = 0, neqr-1
            b0 = kb(ib(i)) + ii
            do j = ib(i), ib(i+1)-1
               do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
                  if (b0 .gt. nzmax) then
                     ierr = 2
                     write (*,*)'csrvbr: no space in b for blk row',i
                     return
                  endif
                  if (a0 .ge. ia(kvstr(i)+ii+1)) then
                     b(b0) = 0.d0
                  else
                     if (jj .eq. ja(a0)) then
                        b(b0) = a(a0)
                        a0 = a0 + 1
                     else
                        b(b0) = 0.d0
                     endif
                  endif
                  b0 = b0 + neqr
               enddo
            enddo
          continue
         enddo
      enddo
      return
      end
      subroutine vbrcsr(ia, ja, a, nr, kvstr, kvstc, ib, jb, kb,
     &   b, nzmax, ierr)
      integer ia(*), ja(*), nr, ib(nr+1), jb(*), kb(*)
      integer kvstr(nr+1), kvstc(*), nzmax, ierr
      real*8  a(*), b(nzmax)
      integer neqr, numc, a0, b0, i, ii, j, jj
      ierr = 0
      a0 = 1
      b0 = 1
      do i = 1, nr
         neqr = kvstr(i+1) - kvstr(i)
         numc = ( kb(ib(i+1)) - kb(ib(i)) ) / neqr
         do j = ib(i), ib(i+1)-1
            do jj = kvstc(jb(j)), kvstc(jb(j)+1)-1
               ja(a0) = jj
               a0 = a0 + 1
            enddo
         enddo
         do ii = 1, neqr-1
            do j = 1, numc
               ja(a0) = ja(a0-numc)
               a0 = a0 + 1
            enddo
         enddo
         a0 = kb(ib(i))
         do ii = 0, neqr-1
            ia(kvstr(i)+ii) = a0
            b0 = kb(ib(i)) + ii
            do jj = 1, numc
               if (a0 .gt. nzmax) then
                  ierr = -(kvstr(i)+ii)
                  write (*,*) 'vbrcsr: no space for row ', -ierr
                  return
               endif
               a(a0) = b(b0)
               a0 = a0 + 1
               b0 = b0 + neqr
            enddo
         enddo
      enddo
      ia(kvstr(nr+1)) = a0
      return
      end
      subroutine csorted(n, ia, ja, sorted)
      integer n, ia(n+1), ja(*)
      logical sorted
      integer i,j
      do i = 1, n
         do j = ia(i)+1, ia(i+1)-1
            if (ja(j-1) .ge. ja(j)) then
               sorted = .false.
               return
            endif
         enddo
      enddo
      sorted = .true.
      return
      end
      


      subroutine submat (job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)
      integer job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      real*8 a(*),ao(*) 
      nr = i2-i1+1
      nc = j2-j1+1
      if ( nr .le. 0 .or. nc .le. 0) return
      klen = 0
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1
         do 60 k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) = j - j1+1
            endif
 60      continue
 100  continue
      iao(nr+1) = klen+1
      return
      end
      subroutine filter(n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)
      real*8 a(*),b(*),drptol
      integer ja(*),jb(*),ia(*),ib(*),n,job,len,ierr
      real*8 norm,loctol
      integer index,row,k,k1,k2 
      index = 1
      do 10 row= 1,n
         k1 = ia(row)
         k2 = ia(row+1) - 1
         ib(row) = index
         goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = 0.0d0
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = 0.0d0
         do 23 k = k1,k2
            if( abs(a(k))  .gt. norm) then
               norm = abs(a(k))
            endif
 23      continue
 400     loctol = drptol * norm
         do 30 k = k1,k2
            if( abs(a(k)) .gt. loctol)then 
               if (index .gt. len) then
               ierr = row 
               return
            endif
            b(index) =  a(k)
            jb(index) = ja(k)
            index = index + 1
         endif
 30   continue
 10   continue
      ib(n+1) = index
      return
      end
      subroutine filterm (n,job,drop,a,ja,b,jb,len,ierr)
      real*8 a(*),b(*),drop
      integer ja(*),jb(*),n,job,len,ierr
      real*8 norm,loctol
      integer index,row,k,k1,k2 
      index = n+2
      do 10 row= 1,n
         k1 = ja(row)
         k2 = ja(row+1) - 1
         jb(row) = index
         goto (100,200,300) job
 100     norm = 1.0d0
         goto 400
 200     norm = a(row)**2 
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = abs(a(row)) 
         do 23 k = k1,k2
            norm = max(abs(a(k)),norm) 
 23      continue
 400     loctol = drop * norm
         do 30 k = k1,k2
            if( abs(a(k)) .gt. loctol)then 
               if (index .gt. len) then
                  ierr = row 
                  return
               endif
               b(index) =  a(k)
               jb(index) = ja(k)
               index = index + 1
            endif
 30      continue
 10   continue
      jb(n+1) = index
      return
      end
      subroutine csort (n,a,ja,ia,values) 
      logical values
      integer n, ja(*), ia(n+1)
      real*8 a(*) 
      integer row, i, k, j
      real*8 rj
      rj = 0.0 
      do row=1, n
         do k=ia(row)+1, ia(row+1)-1 
            j = ja(k)
            if (values) rj = a(k)
            i = k-1;
            do while ((i .ge. ia(row)) .and. (j < ja(i)) ) 
               ja(i+1) = ja(i)
               if (values) a(i+1) = a(i)
               i = i-1
               if (i .eq. 0) exit
            enddo
            ja(i+1) = j
            if (values) a(i+1) = rj
         enddo
      enddo
      return 
      end
      subroutine clncsr(job,value2,nrow,a,ja,ia,indu,iwk)
      integer job, nrow, value2
      integer ia(nrow+1),indu(nrow),iwk(nrow+1),ja(*)
      real*8  a(*)
      integer i,j,k,ko,ipos,kfirst,klast
      real*8  tmp
      if (job.le.0) return
      do 90 i = 1, nrow
         indu(i) = 0
         iwk(i) = ia(i)
 90   continue
      iwk(nrow+1) = ia(nrow+1)
      k = 1
      do 120 i = 1, nrow
         ia(i) = k
         ipos = iwk(i)
         klast = iwk(i+1)
 100     if (ipos.lt.klast) then
            j = ja(ipos)
            if (indu(j).eq.0) then
               if (value2.ne.0) then
                  if (a(ipos) .ne. 0.0D0) then
                     indu(j) = k
                     ja(k) = ja(ipos)
                     a(k) = a(ipos)
                     k = k + 1
                  endif
               else
                  indu(j) = k
                  ja(k) = ja(ipos)
                  k = k + 1
               endif
            else if (value2.ne.0) then
               a(indu(j)) = a(indu(j)) + a(ipos)
            endif
            ipos = ipos + 1
            go to 100
         endif
         do 110 ipos = ia(i), k - 1
            indu(ja(ipos)) = 0
 110     continue
 120  continue
      ia(nrow+1) = k
      if (job.le.1) return
      do 140 i = 1, nrow
         klast = ia(i+1) - 1
         kfirst = ia(i)
 130     if (klast.gt.kfirst) then
            if (ja(klast).lt.i .and. ja(kfirst).ge.i) then
               j = ja(klast)
               ja(klast) = ja(kfirst)
               ja(kfirst) = j
               if (value2.ne.0) then
                  tmp = a(klast)
                  a(klast) = a(kfirst)
                  a(kfirst) = tmp
               endif
            endif
            if (ja(klast).ge.i)
     &         klast = klast - 1
            if (ja(kfirst).lt.i)
     &         kfirst = kfirst + 1
            go to 130
         endif
         if (ja(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
 140  continue
      if (job.le.2) return
      do 190 i = 1, nrow
         do 160 ipos = ia(i), indu(i)-1
            do 150 j = indu(i)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 150        continue
 160     continue
         do 180 ipos = indu(i), ia(i+1)-1
            do 170 j = ia(i+1)-1, ipos+1, -1
               k = j - 1
               if (ja(k).gt.ja(j)) then
                  ko = ja(k)
                  ja(k) = ja(j)
                  ja(j) = ko
                  if (value2.ne.0) then
                     tmp = a(k)
                     a(k) = a(j)
                     a(j) = tmp
                  endif
               endif
 170        continue
 180     continue
 190  continue
      return
      end
      subroutine copmat (nrow,a,ja,ia,ao,jao,iao,ipos,job) 
      real*8 a(*),ao(*) 
      integer nrow, ia(*),ja(*),jao(*),iao(*), ipos, job 
      integer kst, i, k 
      kst    = ipos -ia(1) 
      do 100 i = 1, nrow+1
         iao(i) = ia(i) + kst
 100  continue
      do 200 k=ia(1), ia(nrow+1)-1
         jao(kst+k)= ja(k)
 200  continue
      if (job .ne. 1) return
      do 201 k=ia(1), ia(nrow+1)-1
         ao(kst+k) = a(k)
 201  continue
      return
      end
      subroutine msrcop (nrow,a,ja,ao,jao,job) 
      real*8 a(*),ao(*) 
      integer nrow, ja(*),jao(*), job 
      integer i, k 
      do 100 i = 1, nrow+1
         jao(i) = ja(i) 
 100  continue
      do 200 k=ja(1), ja(nrow+1)-1
         jao(k)= ja(k)
 200  continue
      if (job .ne. 1) return
      do 201 k=ja(1), ja(nrow+1)-1
         ao(k) = a(k)
 201  continue
      do 202 k=1,nrow
         ao(k) = a(k)
 202  continue
      return
      end
      double precision function getelm (i,j,a,ja,ia,iadd,sorted) 
      integer i, ia(*), iadd, j, ja(*)
      double precision a(*)
      logical sorted 
      integer ibeg, iend, imid, k
      iadd = 0 
      getelm = 0.0
      ibeg = ia(i)
      iend = ia(i+1)-1
      if (.not. sorted) then 
         do 5  k=ibeg, iend
            if (ja(k) .eq.  j) then
               iadd = k 
               goto 20 
            endif
 5       continue
      else
 10      imid = ( ibeg + iend ) / 2
         if (ja(imid).eq.j) then
            iadd = imid 
            goto 20
         endif
         if (ibeg .ge. iend) goto 20
         if (ja(imid).gt.j) then
            iend = imid -1
         else 
            ibeg = imid +1
         endif
         goto 10  
      endif
 20   if (iadd .ne. 0) getelm = a(iadd) 
      return
      end 
      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
      real*8 diag(*),a(*)
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
      integer istart, max, iend, i, kold, k, kdiag, ko
      istart = max(0,-ioff)
      iend = min(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
         diag(i) = 0.0d0 
 1    continue
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            endif
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
      ko = 0
      do  7 i=1, nrow 
         kold = ko
         kdiag = idiag(i) 
         do 71 k= ia(i), ia(i+1)-1 
            if (k .ne. kdiag) then
               ko = ko+1
               a(ko) = a(k)
               ja(ko) = ja(k)
            endif
 71      continue
         ia(i) = kold+1
 7    continue
      ia(nrow+1) = ko+1
      return
      end
      subroutine transp (nrow,ncol,a,ja,ia,iwk,ierr)
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      real*8 a(*) 
      real*8 t, t1
      ierr = 0
      nnz = ia(nrow+1)-1
      jcol = 0
      do 1 k=1, nnz
         jcol = max(jcol,ja(k))
 1    continue
      if (jcol .gt. ncol) then
         ierr = jcol
         return
      endif
      ncol = jcol
      do 3 i=1,nrow
         do 2 k=ia(i),ia(i+1)-1
            iwk(k) = i
 2       continue 
 3    continue
      do 35 i=1,ncol+1
         ia(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ja(k)
         ia(i+1) = ia(i+1)+1
 4    continue 
      ia(1) = 1 
      do 44 i=1,ncol
         ia(i+1) = ia(i) + ia(i+1)
 44   continue 
      init = 1
      k = 0
 5    t = a(init)
      i = ja(init)
      j = iwk(init)
      iwk(init) = -1
 6    k = k+1                 
      l = ia(i)
      t1 = a(l)
      inext = ja(l)
      a(l)  = t
      ja(l) = j
      ia(i) = l+1
      if (iwk(l) .lt. 0) goto 65
      t = t1
      i = inext
      j = iwk(l)
      iwk(l) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (iwk(init) .lt. 0) goto 65
      goto 5
 70   continue
      do 80 i=ncol,1,-1 
         ia(i+1) = ia(i)
 80   continue
      ia(1) = 1
      return
      end 
      subroutine getl (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*8 a(*), ao(*)
      real*8 t
      integer ko, kold, kdiag, k, i
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
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
      iao(n+1) = ko+1
      return
      end
      subroutine getu (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*8 a(*), ao(*)
      real*8 t 
      integer ko, k, i, kdiag, kfirst 
      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t
         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
      iao(n+1) = ko+1
      return
      end
      subroutine levels (n, jal, ial, nlev, lev, ilev, levnum)
      integer jal(*),ial(*), levnum(*), ilev(*), lev(*) 
      do 10 i = 1, n
         levnum(i) = 0
 10   continue
      nlev = 0
      do 20 i = 1, n
         levi = 0
         do 15 j = ial(i), ial(i+1) - 1
            levi = max (levi, levnum(jal(j)))
 15      continue
         levi = levi+1 
         levnum(i) = levi 
         nlev = max(nlev,levi) 
 20   continue
      do 21 j=1, nlev+1
         ilev(j) = 0
 21   continue
      do 22 j=1, n
         i = levnum(j)+1
         ilev(i) = ilev(i)+1
 22   continue
      ilev(1) = 1
      do 23 j=1, nlev
         ilev(j+1) = ilev(j)+ilev(j+1)
 23   continue
      do 30 j=1,n
         i = levnum(j)
         lev(ilev(i)) = j
         ilev(i) = ilev(i)+1
 30   continue
      do 35 j=nlev, 1, -1
         ilev(j+1) = ilev(j) 
 35   continue
      ilev(1) = 1
      return
      end
      subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,
     *                  c,jc,ic,iw,nzmax,ierr)
      real*8 a(*),c(*) 
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*),imask(nrow+1) 
      logical iw(ncol)
      ierr = 0
      len = 0
      do 1 j=1, ncol
         iw(j) = .false.
 1    continue
      do 100 ii=1, nrow
         do 2 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
 2       continue
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do 200 k=k1,k2 
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               endif
               jc(len) = j
               c(len) = a(k)
            endif
 200     continue              
         do 3 k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
 3       continue
 100  continue          
      ic(nrow+1)=len+1
      return
      end
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      real*8 a(*),ao(*) 
      logical values
      values = (job .eq. 1) 
      do 50 j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
 50   continue
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
      do 100 ii=1,nrow
         ko = iao(perm(ii)) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
      return
      end
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job) 
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      real*8 a(*), ao(*) 
      integer k, i, nnz
      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k)) 
 100  continue
      if (job .ne. 1) return
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue
      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue
      return
      end
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      real*8 a(*),ao(*) 
      integer locjob, mod
      locjob = mod(job,2) 
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
      locjob = 0
      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob) 
      else 
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob) 
      endif 
      return
      end
      subroutine dperm1 (i1,i2,a,ja,ia,b,jb,ib,perm,ipos,job)
      integer i1,i2,job,ja(*),ia(*),jb(*),ib(*),perm(*)
      real*8 a(*),b(*)
      integer ko,irow,k 
      logical values
      values = (job .eq. 1) 
      ko = ipos 
      ib(1) = ko
      do 900 i=i1,i2
         irow = perm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = ja(k)
            ko=ko+1
 800     continue
         ib(i-i1+2) = ko
 900  continue
      return
      end
      subroutine dperm2 (i1,i2,a,ja,ia,b,jb,ib,cperm,rperm,istart,
     *        ipos,job)
      integer i1,i2,job,istart,ja(*),ia(*),jb(*),ib(*),cperm(*),rperm(*) 
      real*8 a(*),b(*)
      integer ko,irow,k 
      logical values
      values = (job .eq. 1) 
      ko = ipos 
      ib(istart) = ko
      do 900 i=i1,i2
         irow = rperm(i) 
         do 800 k=ia(irow),ia(irow+1)-1
            if (values) b(ko) = a(k)
            jb(ko) = cperm(ja(k))
            ko=ko+1
 800     continue
         ib(istart+i-i1+1) = ko
 900  continue
      return
      end 
      subroutine dmperm (nrow,a,ja,ao,jao,perm,job)
      integer nrow,ja(*),jao(*),perm(nrow),job
      real*8 a(*),ao(*) 
      integer n1, n2
      n1 = nrow+1
      n2 = n1+1
      call dperm (nrow,a,ja,ja,ao(n2),jao(n2),jao,perm,perm,job) 
      jao(1) = n2
      do 101 j=1, nrow 
         ao(perm(j)) = a(j) 
         jao(j+1) = jao(j+1)+n1
 101  continue
      return
      end
      subroutine dvperm (n, x, perm) 
      integer n, perm(n) 
      real*8 x(n)
      real*8 tmp, tmp1
      init      = 1
      tmp       = x(init)        
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
 6    k = k+1
      tmp1      = x(ii) 
      x(ii)     = tmp
      next          = perm(ii) 
      if (next .lt. 0 ) goto 65
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
      goto 6
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp        = x(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
      return
      end
      subroutine ivperm (n, ix, perm) 
      integer n, perm(n), ix(n)
      integer tmp, tmp1
      init      = 1
      tmp        = ix(init)        
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
 6    k = k+1
      tmp1          = ix(ii) 
      ix(ii)     = tmp
      next          = perm(ii) 
      if (next .lt. 0 ) goto 65
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next 
      goto 6
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp        = ix(init)
      ii        = perm(init)
      perm(init)=-perm(init)
      goto 6
 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue 
      return
      end
      subroutine retmx (n,a,ja,ia,dd)
      real*8 a(*),dd(*)
      integer n,ia(*),ja(*)
      integer k2, i, k1, k
      real*8 t, t1, t2
      t = 0.0
      t2 = 0.0
      t1 = 0.0
      k2 = 1
      do 11 i=1,n
         k1 = k2
         k2 = ia(i+1) - 1
         t = 0.0d0
         do 101  k=k1,k2
            t1 = abs(a(k))
            if (t1 .gt. t) t = t1
            if (ja(k) .eq. i) then 
               if (a(k) .ge. 0.0) then 
                  t2 = a(k) 
               else 
                  t2 = - a(k)
               endif
            endif
 101     continue                
         dd(i) =  t2*t
 11   continue
      return
      end
      subroutine diapos  (n,ja,ia,idiag) 
      integer ia(n+1), ja(*), idiag(n) 
      do 1 i=1, n 
         idiag(i) = 0
 1    continue
      do  6 i=1,n
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k) .eq. i) idiag(i) = k
 51      continue
 6    continue
      return
      end
      subroutine dscaldg (n,a,ja,ia,diag,job)
      real*8 a(*), diag(*),t
      integer ia(*),ja(*)
      goto (12,11,10) job+1
 10   do 110 j=1,n
         k1= ia(j)
         k2 = ia(j+1)-1
         t = 0.0d0
         do 111 k = k1,k2
 111        t = t+a(k)*a(k)
 110        diag(j) = sqrt(t)
            goto 12
 11   continue
      call retmx (n,a,ja,ia,diag)
 12   do 1 j=1,n
         if (diag(j) .ne. 0.0d0) then 
            diag(j) = 1.0d0/diag(j)
         else 
            diag(j) = 1.0d0
         endif
 1    continue
      do 2 i=1,n
         t = diag(i)
         do 21 k=ia(i),ia(i+1) -1
            a(k) = a(k)*t
 21      continue
 2    continue
      return 
      end
      subroutine extbdg (n,a,ja,ia,bdiag,nblk,ao,jao,iao)
      implicit real*8 (a-h,o-z)
      real*8 bdiag(*),a(*),ao(*)
      integer ia(*),ja(*),jao(*),iao(*) 
      m = 1 + (n-1)/nblk
      ltr =  ((nblk-1)*nblk)/2 
      l = m * ltr
      do 1 i=1,l
         bdiag(i) = 0.0d0
 1    continue
      ko = 0
      kb = 1
      iao(1) = 1
      do 11 jj = 1,m
         j1 = (jj-1)*nblk+1
         j2 =  min0 (n,j1+nblk-1)
         do 12 j=j1,j2
            do 13 i=ia(j),ia(j+1) -1
               k = ja(i)
               if (k .lt. j1) then
                  ko = ko+1
                  ao(ko) = a(i)
                  jao(ko) = k
               else if (k .lt. j) then
                  bdiag(kb+k-j1) = a(i)
               endif
 13         continue         
            kb = kb + j-j1
            iao(j+1) = ko+1
 12      continue
 11   continue
      return
      end
      subroutine getbwd(n,ja,ia,ml,mu)
      integer ja(*),ia(n+1),ml,mu
      integer ldist,i,k 
      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1 
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
      end
      subroutine blkfnd (nrow,ja,ia,nblk)
      integer ia(nrow+1),ja(*) 
      minlen = ia(2)-ia(1)
      irow   = 1
      do 1 i=2,nrow
         len = ia(i+1)-ia(i)
         if (len .lt. minlen) then
            minlen = len 
            irow = i
         endif
 1    continue
      nblk = 1
      if (minlen .le. 1) return
      do 99 iblk = minlen, 1, -1
         if (mod(minlen,iblk) .ne. 0) goto 99
         len = ia(2) - ia(1)
         len0 = len
         jfirst = ja(1) 
         jlast = ja(ia(2)-1)
         do 10 jrow = irow+1,irow+nblk-1
            i1 = ia(jrow)
            i2 = ia(jrow+1)-1
            len = i2+1-i1
            jf = ja(i1)
            jl = ja(i2) 
            if (len .ne. len0 .or. jf .ne. jfirst .or. 
     *           jl .ne. jlast) goto 99
 10      continue
         call blkchk (nrow,ja,ia,iblk,imsg)               
         if (imsg .eq. 0) then 
            nblk = iblk
            return
         endif 
 99   continue               
      end
      subroutine blkchk (nrow,ja,ia,nblk,imsg)
      integer ia(nrow+1),ja(*) 
      imsg = 0
      if (nblk .le. 1) return
      nr = nrow/nblk
      if (nr*nblk .ne. nrow) goto 101
      irow = 1
      do 20 ii=1, nr
         i1 = ia(irow)
         j2 = i1
         lena = ia(irow+1)-i1
         len = lena/nblk
         if (len* nblk .ne. lena) goto 103
         do 6 i = 1, nblk
            irow = irow + 1
            if (ia(irow)-ia(irow-1) .ne. lena ) goto 104
            do 7 k=0, len-1
               jstart = ja(i1+nblk*k)-1
               if ( (jstart/nblk)*nblk .ne. jstart) goto 102
               do 5 j=1, nblk
                  if (jstart+j .ne. ja(j2) )  goto 104
                  j2 = j2+1 
 5             continue
 7          continue
 6       continue
 20   continue
      return
 101  imsg = -1
      return
 102  imsg = -2
      return
 103  imsg = -3
      return
 104  imsg = -4
      return
      end
      subroutine infdia (n,ja,ia,ind,idiag) 
      integer ia(*), ind(*), ja(*)
      n2= n+n-1
      do 1 i=1,n2
         ind(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            ind(n+j-i) = ind(n+j-i) +1
 2       continue 
 3    continue
      idiag = 0 
      do 41 k=1, n2
         if (ind(k) .ne. 0) idiag = idiag+1
 41   continue
      return
      end
      subroutine amubdg (nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw) 
      integer ja(*),jb(*),ia(nrow+1),ib(ncol+1),ndegr(nrow),iw(ncolb) 
      do 1 k=1, ncolb 
         iw(k) = 0 
 1    continue
      
      do 2 k=1, nrow
         ndegr(k) = 0 
 2    continue
      do 7 ii=1,nrow 
         ldg = 0 
         last = -1 
         do 6 j = ia(ii),ia(ii+1)-1 
            jr = ja(j) 
            do 5 k=ib(jr),ib(jr+1)-1
               jc = jb(k) 
               if (iw(jc) .eq. 0) then 
                  ldg = ldg + 1
                  iw(jc) = last 
                  last = jc
               endif
 5          continue
 6       continue
         ndegr(ii) = ldg
         do 61 k=1,ldg 
            j = iw(last) 
            iw(last) = 0
            last = j
 61      continue
 7    continue
      nnz = 0
      do 8 ii=1, nrow 
         nnz = nnz+ndegr(ii) 
 8    continue
      return
      end
      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw) 
      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow) 
      do 1 k=1, ncol 
         iw(k) = 0 
 1    continue
      do 2 k=1, nrow
         ndegr(k) = 0 
 2    continue
      do 7 ii=1,nrow 
         ldg = 0 
         last = -1 
         do 5 j = ia(ii),ia(ii+1)-1 
            jr = ja(j) 
            ldg = ldg + 1
            iw(jr) = last 
            last = jr
 5       continue
         do 6 j=ib(ii),ib(ii+1)-1
            jc = jb(j)
            if (iw(jc) .eq. 0) then 
               ldg = ldg + 1
               iw(jc) = last 
               last = jc
            endif
 6       continue
         ndegr(ii) = ldg
         do 61 k=1,ldg 
            j = iw(last) 
            iw(last) = 0
            last = j
 61      continue
 7    continue
      nnz = 0
      do 8 ii=1, nrow 
         nnz = nnz+ndegr(ii) 
 8    continue
      return
      end
      subroutine rnrms (nrow, nrm, a, ia, diag) 
      real*8 a(*), diag(nrow), scal 
      integer ia(nrow+1) 
      do 1 ii=1,nrow
         scal = 0.0d0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) ) 
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) ) 
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         endif 
         if (nrm .eq. 2) scal = sqrt(scal) 
         diag(ii) = scal
 1    continue
      return
      end 
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag) 
      real*8 a(*), diag(nrow) 
      integer ja(*), ia(nrow+1) 
      do 10 k=1, nrow 
         diag(k) = 0.0d0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k) 
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) ) 
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) ) 
            else
               diag(j) = diag(j)+a(k)**2
            endif 
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
      end 
      subroutine roscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      real*8 a(*), b(*), diag(nrow) 
      integer nrow,job,nrm,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
      call rnrms (nrow,nrm,a,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0d0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
      end
      subroutine coscal(nrow,job,nrm,a,ja,ia,diag,b,jb,ib,ierr) 
      real*8 a(*),b(*),diag(nrow) 
      integer nrow,job,ja(*),jb(*),ia(nrow+1),ib(nrow+1),ierr 
      call cnrms (nrow,nrm,a,ja,ia,diag)
      ierr = 0
      do 1 j=1, nrow
         if (diag(j) .eq. 0.0) then
            ierr = j 
            return
         else
            diag(j) = 1.0d0/diag(j)
         endif
 1    continue
      call amudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
      end
      subroutine addblk(nrowa, ncola, a, ja, ia, ipos, jpos, job,
     & nrowb, ncolb, b, jb, ib, nrowc, ncolc, c, jc, ic, nzmx, ierr)
      integer nrowa, nrowb, nrowc, ncola, ncolb, ncolc, ipos, jpos
      integer nzmx, ierr, job
      integer ja(1:*), ia(1:*), jb(1:*), ib(1:*), jc(1:*), ic(1:*)
      real*8 a(1:*), b(1:*), c(1:*)
      logical values
      integer i,j1,j2,ka,kb,kc,kamax,kbmax
      kamax = 0 
      
      values = (job .ne. 0) 
      ierr = 0
      nrowc = max(nrowa, nrowb+ipos-1)
      ncolc = max(ncola, ncolb+jpos-1)
      kc = 1
      kbmax = 0
      ic(1) = kc
      do 10 i=1, nrowc
         if (i.le.nrowa) then
            ka = ia(i)
            kamax = ia(i+1)-1
         else
            ka = ia(nrowa+1)
         end if
         if ((i.ge.ipos).and.((i-ipos).le.nrowb)) then
            kb = ib(i-ipos+1)
            kbmax = ib(i-ipos+2)-1 
         else
            kb = ib(nrowb+1)
         end if
 20      continue 
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncolc+1
         endif
         if (kb .le. kbmax) then 
            j2 = jb(kb) + jpos - 1
         else 
            j2 = ncolc+1
         endif
         if ((ka .le. kamax .or. kb .le. kbmax) .and.
     &        (j1 .le. ncolc .or. j2 .le. ncolc)) then
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
            if (kc .gt. nzmx) goto 999
            goto 20
         end if
         ic(i+1) = kc
 10   continue
      return
 999  ierr = i 
      return
      end
      subroutine get1up (n,ja,ia,ju)
      integer  n, ja(*),ia(*),ju(*)
      integer i, k 
      do 5 i=1, n
         ju(i) = 0
         k = ia(i) 
 1       continue
         if (ja(k) .ge. i) then
            ju(i) = k
            goto 5
         elseif (k .lt. ia(i+1) -1) then
            k=k+1
            goto 1
         endif 
 5    continue
      return
      end 
      subroutine xtrows (i1,i2,a,ja,ia,ao,jao,iao,iperm,job)
      integer i1,i2,ja(*),ia(*),jao(*),iao(*),iperm(*),job
      real*8 a(*),ao(*) 
      logical values
      values = (job .eq. 1) 
      ko = 1
      iao(1) = ko
      do 100 j=i1,i2 
         ii = iperm(j) 
         do 60 k=ia(ii), ia(ii+1)-1 
            jao(ko) = ja(k) 
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
         iao(j-i1+2) = ko
 100  continue
      return
      end
      subroutine csrkvstr(n, ia, ja, nr, kvstr)
      integer n, ia(n+1), ja(*), nr, kvstr(*)
      integer i, j, jdiff
      nr = 1
      kvstr(1) = 1
      do i = 2, n
         jdiff = ia(i+1)-ia(i)
         if (jdiff .eq. ia(i)-ia(i-1)) then
            do j = ia(i), ia(i+1)-1
               if (ja(j) .ne. ja(j-jdiff)) then
                  nr = nr + 1
                  kvstr(nr) = i
                  goto 299
               endif
            enddo
 299        continue
         else
            nr = nr + 1
            kvstr(nr) = i
         endif
      enddo
      kvstr(nr+1) = n+1
      return
      end
      subroutine csrkvstc(n, ia, ja, nc, kvstc, iwk)
      integer n, ia(n+1), ja(*), nc, kvstc(*), iwk(*)
      integer i, j, k, ncol
      ncol = 0
      do i = 1, n
         if (ia(i) .lt. ia(i+1)) then
            j = ja(ia(i))
            iwk(j) = 1
            do k = ia(i)+1, ia(i+1)-1
               j = ja(k)
               if (ja(k-1).ne.j-1) then
                  iwk(j) = 1
                  iwk(ja(k-1)+1) = 1
               endif
            enddo
            iwk(j+1) = 1
            ncol = max0(ncol, j)
         endif
      enddo
      nc = 1
      kvstc(1) = 1
      do i = 2, ncol+1
         if (iwk(i).ne.0) then
            nc = nc + 1
            kvstc(nc) = i
            iwk(i) = 0
         endif
      enddo
      nc = nc - 1
      return
      end
      subroutine kvstmerge(nr, kvstr, nc, kvstc, n, kvst)
      integer nr, kvstr(nr+1), nc, kvstc(nc+1), n, kvst(*)
      integer i,j
      if (kvstr(nr+1) .ne. kvstc(nc+1)) return
      i = 1
      j = 1
      n = 1
  200 if (i .gt. nr+1) then
         kvst(n) = kvstc(j)
         j = j + 1
      elseif (j .gt. nc+1) then
         kvst(n) = kvstr(i)
         i = i + 1
      elseif (kvstc(j) .eq. kvstr(i)) then
         kvst(n) = kvstc(j)
         j = j + 1
         i = i + 1
      elseif (kvstc(j) .lt. kvstr(i)) then
         kvst(n) = kvstc(j)
         j = j + 1
      else
         kvst(n) = kvstr(i)
         i = i + 1
      endif
      n = n + 1
      if (i.le.nr+1 .or. j.le.nc+1) goto 200
      n = n - 2
      return
      end

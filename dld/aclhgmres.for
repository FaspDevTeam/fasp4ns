ccccc   Last edit time:  202203171451
      subroutine acln(neq,numel,numcol,nt1)
      implicit real*8 (a-h,o-z)
        dimension numcol(*),lm(1000)
6       format(1x,10i7)
      if(numel.eq.0) go to 1000
      do 100 n=1,neq
        numcol(n)=1
100     continue
        do 500 ke=1,numel
      read(nt1) nd,(lm(i),i=1,nd)
c       write (*,*) 'nd= ',nd, (lm(i),i=1,nd)
      do 400 i=1,nd
      ni = lm(i)
      do 300 j=1,nd
      nj = lm(j)
        if (nj.eq.ni) goto 300
        numcol(ni) = numcol(ni)+1
300     continue
400     continue
500     continue
c       write(*,*) 'numcol ='
c       write(*,6) (numcol(n),n=1,neq)
c        do 600 i=1,neq
c        numcol(neq-i+2)=numcol(neq-i+1)
c600     continue
c        numcol(1)=0
c        do 700 i=1,neq
c        numcol(i+1)=numcol(i+1)+numcol(i)
c700     continue
      do 800 n=1,neq-1
800     numcol(n+1) = numcol(n+1)+numcol(n)
      do 850 n=1,neq
850     numcol(neq-n+2) = numcol(neq-n+1)
      numcol(1) = 0
        rewind (nt1)
1000    return
      end
 
        subroutine aclh(neq,numel,jnz,numcol,mht,nt1,maxt)
      implicit real*8 (a-h,o-z)
        dimension mht(maxt),jnz(*),numcol(*),lm(1000)
6       format(1x,10i7)
      if(numel.eq.0) go to 1000
        call  acln(neq,numel,jnz,nt1)
        if (jnz(neq+1).gt.maxt) then
        write(*,*) 'exceet core memory ....',jnz(neq+1),' > ',maxt
        stop 1111
        endif
        do 100 n=1,neq
      numcol(n)=0
100     continue
      do 450 ke=1,numel
      read(nt1) nd,(lm(i),i=1,nd)
c       write (*,*) 'nd= ',nd, (lm(i),i=1,nd)
      do 400 i=1,nd
      ni = lm(i)
      do 300 j=1,nd
c       if (j.eq.i) goto 300
      nj = lm(j)
      do 200 k=1,numcol(ni)
        if (nj.eq.mht(k+jnz(ni))) goto 300
200     continue
      numcol(ni) = numcol(ni)+1
        if (ni.gt.neq) then
        write(*,*) 'ni ====',ni
        stop 1234
        endif
        mht(numcol(ni)+jnz(ni)) = nj
300     continue
400     continue
450     continue
c       write(*,*) 'numcol ='
c       write(*,6) (numcol(n),n=1,neq)
c       write(*,*) 'mht ='
c       do 2000 n=1,neq
c2000   write(*,6) (mht(i+jnz(n)),i=1,numcol(n))
      do 600 n=1,neq
      nd = numcol(n)
      do 500 i=1,nd
500     lm(i) = mht(i+jnz(n))
      call order(nd,lm)
      do 550 i=1,nd
550     mht(i+jnz(n)) = lm(i)
600     continue
1000    return
      end
 
        subroutine bclh(neq,numel,jnz,numcol,na,maxt)
      implicit real*8 (a-h,o-z)
        dimension jnz(*),numcol(*),na(*)
      if(numel.eq.0) go to 1000
cccc        call aclh(neq,numel,jnz,numcol,na,nt1,maxt)
        nsum = 0
      do 700 n=1,neq
        nn=jnz(n)
      do 700 i=1,numcol(n)
      nsum = nsum+1
      na(nsum) = na(nn+i)
700     continue
      do 800 n=1,neq-1
800     numcol(n+1) = numcol(n+1)+numcol(n)
      do 850 n=1,neq
850     numcol(neq-n+2) = numcol(neq-n+1)
      numcol(1) = 0
c       write(*,*) 'nsum,numcol(neq+1) =',nsum,numcol(neq+1)
c       write(*,*) 'numcol ='
c       write(*,6) (numcol(n),n=1,neq+1)
c       write(*,*) 'na ='
c       write(*,6) (na(i),i=1,nsum)
1000    return
6       format(1x,5i15)
      end
 
      subroutine order(nd,lm)
      implicit real*8 (a-h,o-z)
        dimension lm(*)
c       write(*,*) '**** order ****'
c       write(*,*) (lm(i),i=1,nd)
      do 200 i=1,nd
      ls=lm(i)+1
      do 100 j=i,nd
      if (lm(j).gt.ls) goto 100
      ls=lm(j)
      j0=j
100     continue
      lm(j0)=lm(i)
      lm(i)=ls
200     continue
c       write(*,*) (lm(i),i=1,nd)
c       write(*,*) '-----------------'
      return
      end
      
      
      
      subroutine adda(na,a,numcol,estif,neq,maxa,numel,nt1,nt2)
      implicit real*8 (a-h,o-z)
        dimension a(maxa),na(maxa),numcol(*),lm(1000),estif(200,200)
      if(numel.eq.0) go to 500
      do 100 i=1,maxa
100   a(i) = 0.0
      do 400 ke=1,numel
      read(nt1) nd,(lm(i),i=1,nd)
      read(nt2) ((estif(i,j),j=1,nd),i=1,nd)
c     write (*,*) nercr,nlrec,lrd,nd, (lm(i),i=1,nd)
c     write (*,*) ((estif(i,j),j=1,i),i=1,nd)
      do 300 i=1,nd
      ii = lm(i)
      n0 = numcol(ii)+1
      n1 = numcol(ii+1)
c     write(*,*) 'n0,n1 =',n0,n1
      do 280 j=1,nd
      jj = lm(j)
      do 240 k=n0,n1
      if (na(k).eq.jj) goto 280
240   continue
280   a(k) = a(k) + estif(i,j)
300   continue
400   continue
500   return
      end
 


      subroutine getname(name,IT)
      implicit real*8 (a-h,o-z)
      character name*50,ch3*3
c     IF (IT.LT.10) WRITE(UNIT=CH3,FMT='(I1)') IT
c     IF (IT.GE.10) WRITE(UNIT=CH3,FMT='(I2)') IT
c     IF (IT.GE.100) WRITE(UNIT=CH3,FMT='(I3)') IT
      call getext(it,ch3)
c     write(*,*) 'name =',name,'++++ CH3 =',CH3
      do 10 i=1,50
      if (name(i:i).eq.' ') then
      j=i
      goto 20
      endif
10    continue
20    continue
      if (j.gt.40) then
      write(*,*) 'Error, plot filename too long .......',name
      write(*,*) ' the length of name must be less or equal 8 character'
      stop 111
      endif
c     read(*,'(a3)') ch3
      name(j:j)='.'
      name(j+1:j+4)=ch3
c     write(*,*) 'plot filename = ',name
      return
      end
 
      subroutine getext(ii,ch3)
      implicit real*8 (a-h,o-z)
      character ch3*3
      it = ii
      ch3 = '   '
      k = 0
      if (ii.ge.100) then
      n = it/100
      k = k+1
      call getchar(n,k,ch3)
      it = it - n*100
      endif
      if (ii.ge.10) then
      n = it/10
      k = k+1
      call getchar(n,k,ch3)
      it = it - n*10
      endif
      n = it
      k = k+1
      call getchar(n,k,ch3)
      return
      end
 
      subroutine getchar(n,k,ch3)
      implicit real*8 (a-h,o-z)
      character ch3*3
      if (n.eq.0) ch3(k:k) = '0'
      if (n.eq.1) ch3(k:k) = '1'
      if (n.eq.2) ch3(k:k) = '2'
      if (n.eq.3) ch3(k:k) = '3'
      if (n.eq.4) ch3(k:k) = '4'
      if (n.eq.5) ch3(k:k) = '5'
      if (n.eq.6) ch3(k:k) = '6'
      if (n.eq.7) ch3(k:k) = '7'
      if (n.eq.8) ch3(k:k) = '8'
      if (n.eq.9) ch3(k:k) = '9'
      return
      end

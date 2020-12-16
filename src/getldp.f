! getldp.f 
! Generate LDPs for main-chain probability map
! Copyright (C) 2020 Jiahua He

! This program is modified from MAINMAST.f (http://kiharalab.org/mainmast/) under GNU General Public License version 3

! Copyright (C) 2017 Genki Terashi, Daisuke Kihara 
!                    and Purdue University

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
*Technology

      module params
      parameter(msize=300)
      parameter(mlink=10)

      real gw,threshold,rshift,filter,rmerge
      end module params

      module commonmap
      use params,only:msize
      implicit integer (i-z)
      character*100 fname_mc,fname_ca
      real gstep,gbase(3)
      real amap(msize,msize,msize),cmap(msize,msize,msize)
      integer ngrid(3)

      integer sn,sn0
      dimension posi(3,msize**3)
      dimension fin(msize**3),fori(msize**3)
      dimension member(msize**2*10)
      end module commonmap

      program main
      use params
      use commonmap
      implicit integer (i-z)
      dimension ftmp(3)
      dimension after(3,msize**2*10),before(3,msize**2*10)
      character*100 argv(0:30)
      logical existance
      
      write(*,*) "# GETLDP.F"
      write(*,*) "# Generate LDPs for main-chain probability map"

      narg=iargc()
      do k=0,narg
        call getarg(k,argv(k))
      end do

      if(narg.lt.2) then
        write(*,*) " Usage: "//trim(argv(0))//
     *             " MC.situs CA.situs [OPTIONS]"
        write(*,*) " OPTIONS:"
        write(*,*) " -gw     [f]: bandwidth of the gaussian filter"
        write(*,*) "              def=2.0"
        write(*,*) " -thresh [f]: threshold of reading density values"
        write(*,*) "              def=0.0"
        write(*,*) " -rshift [f]: maximum shift distance"
        write(*,*) "              def=10.0"
        write(*,*) " -filter [f]: drop percentage of low density data"
        write(*,*) "              def=0.0"
        write(*,*) " -rmerge [f]: merging radius after mean-shifting"
        write(*,*) "              def=0.5"
        stop
      end if

      fname_mc=argv(1)
      inquire(file=fname_mc,exist=existance)
      if(.not.existance) then
        write(*,*) "# MC map file "//trim(fname_mc)//" not exist!"
        stop
      end if

      fname_ca=argv(2)
      inquire(file=fname_ca,exist=existance)
      if(.not.existance) then
        write(*,*) "# CA map file "//trim(fname_ca)//" not exist!"
        stop
      end if

      mode=0
      gw=2.0
      threshold=0.0
      rshift=10.0
      filter=0.0
      rmerge=0.5

      k=3
      do while(k.le.narg)
        if(argv(k).eq.'-gw') then
          if(k+1.le.narg) read(argv(k+1),*) gw
          k=k+2
        else if(argv(k).eq.'-thresh') then
          if(k+1.le.narg) read(argv(k+1),*) threshold
          k=k+2
        else if(argv(k).eq.'-rshift') then
          if(k+1.le.narg) read(argv(k+1),*) rshift
          k=k+2
        else if(argv(k).eq.'-filter') then
          if(k+1.le.narg) read(argv(k+1),*) filter
          k=k+2
        else if(argv(k).eq.'-rmerge') then
          if(k+1.le.narg) read(argv(k+1),*) rmerge
          k=k+2
        else 
          write(*,*) "# Unknown parameter input"
          stop
        end if
      end do

      write(*,*) "# Gaussian width=",gw
      write(*,*) "# Threshold of reading density value=",threshold
      write(*,*) "# Maximum shift distance=",rshift
      write(*,*) "# Drop percentage of low density data=",filter
      write(*,*) "# Merging radius after mean-shifting=",rmerge

      write(*,*) "# Reading MC map from file "//trim(fname_mc)
      call read_mc

      write(*,*) "# Reading CA map from file "//trim(fname_ca)
      call read_ca

      sn=0
      do k=1,ngrid(3)
        do j=1,ngrid(2)
          do i=1,ngrid(1)
            if(amap(i,j,k).gt.threshold) then
              sn=sn+1
              posi(1:3,sn)=(/i,j,k/)
              fin(sn)=amap(i,j,k)
            end if
          end do
        end do
      end do

      write(*,*) "# Initial SN=",sn

      sn0=sn

      fin(1:sn)=0.0
      do i=1,sn
        call meanshift(float(posi(1:3,i)),after(1:3,i),fin(i))
      end do

      before(1:3,1:sn)=after(1:3,1:sn)
      fori(1:sn)=fin(1:sn)

      write(*,*) "# Before merging SN0=",sn
      call  mergepoints(before(1:3,1:sn),after(1:3,1:sn))
      write(*,*) "# After Merging SN=",sn

      call ca_ldp_kmean(after(1:3,1:sn))

      call gen_ldp(after(1:3,1:sn),posi(1:3,1:sn0))

      stop 
      end program main

       
      subroutine read_mc
      use commonmap
      implicit integer (i-z)
      dimension buff(10)

      open(100,file=fname_mc,status='old')

      read(100,*) gstep,(gbase(k),k=1,3),(ngrid(k),k=1,3)

      write(*,*) "# Grid step=",gstep
      write(*,*) "# Grid origin=",gbase(:)
      write(*,*) "# Number of grids=",ngrid(:)

      if(maxval(ngrid(:)).gt.msize) then
        write(*,*) "#ERROR !! map size > max size !"
        stop
      end if

      maxlines=ngrid(1)*ngrid(2)*ngrid(3)/10
      write(*,*) '# Number of lines in .situs file=',maxlines

      read(100,*) ! spacing

      i=1
      j=1
      k=1
      do l=1,maxlines
        read(100,*) (buff(m),m=1,10)
        do m=1,10
          if(buff(m).gt.0) then !ignore negatve value
            amap(i,j,k)=buff(m)
          end if
          i=i+1
          if(i.gt.ngrid(1)) then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2)) then
            j=1
            k=k+1
          end if
        end do
      end do

      remain=ngrid(1)*ngrid(2)*ngrid(3)-maxlines*10
      if(remain.gt.0) then
        read(100,*) (buff(m),m=1,remain)
        do m=1,remain
          if(buff(m).gt.0) then
            amap(i,j,k)=buff(m)
          end if
          if(buff(m).gt.threshold) then
            sn=sn+1
            posi(1:3,sn)=(/i,j,k/)
            fin(sn)=buff(m)
          end if

          i=i+1
          if(i.gt.ngrid(1)) then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2)) then
            j=1
            k=k+1
          end if
        end do
      end if

      close(100)
     
      return  
      end subroutine read_mc

      subroutine read_ca
      use commonmap
      implicit integer (i-z)
      dimension buff(10)
      real gstep_tmp,gbase_tmp(3)
      integer ngrid_tmp(3)

      open(200,file=fname_ca,status='old')

      read(200,*) gstep_tmp,(gbase_tmp(k),k=1,3),(ngrid_tmp(k),k=1,3)

      if(gstep_tmp   .ne.gstep   .or.
     *   gbase_tmp(1).ne.gbase(1).or.
     *   gbase_tmp(2).ne.gbase(2).or.
     *   gbase_tmp(3).ne.gbase(3).or.
     *   ngrid_tmp(1).ne.ngrid(1).or.
     *   ngrid_tmp(2).ne.ngrid(2).or.
     *   ngrid_tmp(3).ne.ngrid(3)) then
        write(*,*) "# MC map and CA map not match!"
        stop
      end if
        
      maxlines=ngrid(1)*ngrid(2)*ngrid(3)/10

      read(200,*) ! spacing

      i=1
      j=1
      k=1
      do l=1,maxlines
        read(200,*) (buff(m),m=1,10)
        do m=1,10
          if(buff(m).gt.0) then !ignore negatve value
            cmap(i,j,k)=buff(m)
          end if
          i=i+1
          if(i.gt.ngrid(1)) then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2)) then
            j=1
            k=k+1
          end if
        end do
      end do

      remain=ngrid(1)*ngrid(2)*ngrid(3)-maxlines*10
      if(remain.gt.0) then
        read(200,*) (buff(m),m=1,remain)
        do m=1,remain
          if(buff(m).gt.0) then
            cmap(i,j,k)=buff(m)
          end if
          i=i+1
          if(i.gt.ngrid(1)) then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2)) then
            j=1
            k=k+1
          end if
        end do
      end if

      close(200)
     
      return  
      end subroutine read_ca

      subroutine ca_ldp_kmean(af)
      use params,only:msize
      use commonmap
      implicit integer(i-z)
      dimension af(3,msize**2*10),ftmp(3)

      do i=1,sn
        a=af(1,i)
        b=af(2,i)
        c=af(3,i)
        x0=floor(a)
        y0=floor(b)
        z0=floor(c)
        x1=x0+1
        y1=y0+1
        z1=z0+1
        dx=a-x0
        dy=b-y0
        dz=c-z0
        ca_prob=   dx    *dy    *dz    *cmap(x1,y1,z1)
     *            +(1-dx)*dy    *dz    *cmap(x0,y1,z1)
     *            +dx    *(1-dy)*dz    *cmap(x1,y0,z1)
     *            +dx    *dy    *(1-dz)*cmap(x1,y1,z0)
     *            +dx    *(1-dy)*(1-dz)*cmap(x1,y0,z0)
     *            +(1-dx)*dy    *(1-dz)*cmap(x0,y1,z0) 
     *            +(1-dx)*(1-dy)*dz    *cmap(x0,y0,z1) 
     *            +(1-dx)*(1-dy)*(1-dz)*cmap(x0,y0,z0)

        ftmp(1:3)=(af(1:3,i)-1.00)*gstep+gbase
        natm=natm+1
        write(*,'("ATOM  ",i5,"  CA  ALA A",
     *            i4,"    ",3f8.3,f6.2)'),
     *  natm,natm,ftmp,ca_prob/100.0
      end do
        
      return
      end subroutine 

      subroutine gen_ldp(af,bf)
      use params
      use commonmap
      implicit integer(i-z)
      dimension af(3,msize*msize*10)
      integer bf(3,msize*msize*10)
      dimension ftmp(3),itmp(3),fvec(3)
      logical flag
      integer cid(sn),nc(sn),contact(sn,sn)
      dimension edge(sn*sn),ide(2,sn*sn),order(sn*sn),
     *          edge_dens(sn*sn)
      logical activeg(msize**2*10)

      ! check LDP contact
      contact(:,:)=0
      do i=1,sn0
        p1=member(i)
        if(p1.eq.0) cycle
        do j=i+1,sn0
          p2=member(j)
          if(p2.eq.0) cycle
          if(p1.eq.p2)cycle
          itmp=abs(bf(1:3,i)-bf(1:3,j))
          if(maxval(itmp(1:3)).eq.1) then
            contact(p1,p2)=contact(p1,p2)+1
            contact(p2,p1)=contact(p2,p1)+1
          end if
        end do
      end do

      ! calculate contacting edges
      nedge=0
      do i=1,sn
        cid(i)=i
        do j=i+1,sn
          if(contact(i,j).eq.0) cycle
          ftmp=af(:,i)-af(:,j)
          dist=sqrt(dot_product(ftmp,ftmp))
          nedge=nedge+1
          edge(nedge)=dist
          ide(1,nedge)=i
          ide(2,nedge)=j
        end do
      end do
      write(*,*),"# Number of edges=",nedge

      call quick_sort(edge,order,nedge)

      do i=1,sn
        cid(i)=i
      end do
      
      ! build MST
      do i=1,nedge
        me=order(i) ! min edge
        v1=ide(1,me)
        v2=ide(2,me)
        if(cid(v1).ne.cid(v2)) then
          tmpid=cid(v2)
          do j=1,sn
            if(cid(j).eq.tmpid) cid(j)=cid(v1)
          end do
        end if
      end do

      ! remove island data
      nc(:)=0
      do i=1,sn
        nc(cid(i))=nc(cid(i))+1
      end do
      maxc=maxval(nc)
      do i=1,sn
        if(nc(i).eq.maxc) then
          usechain=i
          exit
        end if
      end do

      ! calculate edges again
      nedge=0
      do i=1,sn
        if(cid(i).ne.usechain) cycle
        do j=i+1,sn
          if(contact(i,j).eq.0) cycle
          if(cid(j).ne.usechain) cycle
          ftmp=af(:,i)-af(:,j)
          dist=sqrt(dot_product(ftmp,ftmp))
          nedge=nedge+1
          edge(nedge)=dist !distance
          ide(1,nedge)=i
          ide(2,nedge)=j
        end do
      end do

      ! calculate density weighted edge
      do i=1,nedge
        v1=ide(1,i)
        v2=ide(2,i)
        fvec=af(:,v1)-af(:,v2) ! vector v2->v1
        density_min=999999999.00
        do j=1,9
          ftmp=af(:,v2)+fvec*float(j)*0.1
          call meanshift_pos(ftmp,density)
          if(density.lt.density_min) density_min=density
        end do
        edge_dens(i)=density_min*edge(i)
      end do

      write(*,*) "# Maximum number of nodes used=",maxc
      write(*,*) "# Using root node id=",usechain
      write(*,*) "# Sum of edge density=",sum(edge_dens(1:nedge))
      write(*,*) "# Avg of edge density=",sum(edge_dens(1:nedge))/nedge

      do i=1,nedge
        write(*,'("BOND ",2i6,2f12.5)'),
     *          ide(1:2,i),edge(i),edge_dens(i)
      end do

      return
      end subroutine gen_ldp

      subroutine mergepoints(before,after)
      use params,only:msize,rmerge,filter
      use commonmap
      implicit integer(i-z)

      dimension before(3,msize**2*10),after(3,msize**2*10)
      logical stock(msize**2*10)
      dimension ftmp(3),member_tmp(msize**2*10)

      stock(1:sn0)=.true.

      dist=rmerge/gstep
      dist2=dist**2

      fmax=maxval(fori(1:sn0))
      fmin=minval(fori(1:sn0))
      frange=fmax-fmin

      do i=1,sn0
        member(i)=i
      end do

      do i=1,sn0-1
        if((fori(i)-fmin)/frange.lt.filter) stock(i)=.false.
        if(.not.stock(i)) cycle
        do j=i+1,sn0
          if(.not.stock(j)) cycle
          if((fori(j)-fmin)/frange.lt.filter) stock(j)=.false.
          if(.not.stock(j)) cycle
          ! check dist
          ftmp=before(1:3,i)-before(1:3,j)
          check=dot_product(ftmp,ftmp)
          if(check.lt.dist2) then
            if(fori(i).gt.fori(j)) then !keep i
              stock(j)=.false.
              member(j)=i
            else !keep j
              stock(i)=.false.
              member(i)=j
              exit
            end if
          end if
        end do
      end do
        
      do i=1,sn0
        now=member(i)
        do j=1,sn0
          if(now.eq.member(now)) exit
          now=member(now)
        end do
        member(i)=now
      end do

      ! copy
      ncnt=1
      member_tmp(:)=0
      do i=1,sn0
        if(stock(i)) then
          after(:,ncnt)=before(:,i)
          fin(ncnt)=fori(i)
          member_tmp(i)=ncnt
          ncnt=ncnt+1
        end if
      end do

      do i=1,sn0
        member(i)=member_tmp(member(i))
      end do

      sn=ncnt-1
        
      return
      end subroutine mergepoints

      subroutine meanshift(before,after,fout)
      use params,only:gw,rshift
      use commonmap
      implicit integer(i-z)
      dimension after(3),before(3),ftmp(3),tmpl(3),tmpu(3)

      nshift=0
      fsigma=(gw/gstep)*0.5
      fsigma=fsigma**2
      fsigma=1.000/fsigma
      fmaxd=(gw/gstep)*2.0
      dist2=(rshift/gstep)**2

      ftmp=before
      do while (nshift.le.10000)
        tmpl=nint(ftmp-fmaxd)
        tmpu=nint(ftmp+fmaxd)
        do ii=1,3
          if(tmpl(ii).lt.1) tmpl(ii)=1
          if(tmpu(ii).gt.ngrid(ii)) tmpu(ii)=ngrid(ii)
        end do

        ftotal=0.0
        after(:)=0.0
        do k=tmpl(3),tmpu(3)
          fz=ftmp(3)-float(k)
          fz2=fz**2
          do j=tmpl(2),tmpu(2)
            fy=ftmp(2)-float(j)
            fy2=fy**2
            do i=tmpl(1),tmpu(1)
              fx=ftmp(1)-float(i)
              fx2=fx**2

              fr2=fz2+fy2+fx2
              gaussian=exp(-1.50*fr2*fsigma)*amap(i,j,k)
              after=after+gaussian*(/i,j,k/)
              ftotal=ftotal+gaussian
            end do
          end do
        end do

        if(ftotal.eq.0.000) exit
        after=after/ftotal
        ftmp=ftmp-after
        if(dot_product(ftmp,ftmp).lt.0.001) exit
        ftmp=after-before
        if(dot_product(ftmp,ftmp).gt.dist2) exit

        ftmp=after
        nshift=nshift+1
      end do

      fout=ftotal

      return
      end subroutine meanshift

      subroutine meanshift_pos(before,fout)
      use params,only:gw
      use commonmap
      implicit integer(i-z)
      dimension before(3),ftmp(3),tmpl(3),tmpu(3)

      ! gaussian kernel gw=window size
      fsigma=(gw/gstep)*0.5
      fsigma=fsigma**2
      fsigma=1.0/fsigma
      fmaxd=(gw/gstep)*2.0

      ftmp=before

      tmpl=nint(ftmp-fmaxd)
      tmpu=nint(ftmp+fmaxd)
      do i=1,3
        if(tmpl(i).lt.1) tmpl(i)=1
        if(tmpu(i).gt.ngrid(i)) tmpu(i)=ngrid(i)
      end do

      ftotal=0
      do k=tmpl(3),tmpu(3)
        fz=ftmp(3)-float(k)
        fz2=fz**2
        do j=tmpl(2),tmpu(2)
          fy=ftmp(2)-float(j)
          fy2=fy**2
          do i=tmpl(1),tmpu(1)
            fx=ftmp(1)-float(i)
            fx2=fx**2

            fr2=fz2+fy2+fx2
            gaussian=exp(-1.50*fr2*fsigma)*amap(i,j,k)
            ftotal=ftotal+gaussian
          end do
        end do
      end do

      fout=ftotal

      return
      end subroutine meanshift_pos

      recursive subroutine quick_sort(list, order,n)
      use params,only:msize
      implicit none
      real  list(msize*msize*msize)
      integer order(msize*msize*msize)

      ! local variable
      integer :: i
      integer :: n

      do i = 1, n
        order(i) = i
      end do

      call quick_sort_1(1, n)

      contains

        recursive subroutine quick_sort_1(left_end, right_end)
        
        integer, intent(in) :: left_end, right_end

        !     local variables
        integer             :: i, j, itemp
        real                :: reference, temp
        integer, parameter  :: max_simple_sort_size = 6

        if (right_end < left_end + max_simple_sort_size) then
          ! use interchange sort for small lists
          call interchange_sort(left_end, right_end)

        else
          ! use partition ("quick") sort
          reference = list((left_end + right_end)/2)
          i = left_end - 1
          j = right_end + 1

          do
            ! scan list from left end until element >= reference is
            ! found
            do
              i = i + 1
              if (list(i) >= reference) exit
            end do
            ! scan list from right end until element <= reference is
            ! found
            do
              j = j - 1
              if (list(j) <= reference) exit
            end do

            if (i < j) then
              ! swap two out-of-order elements
              temp = list(i)
              list(i) = list(j)
              list(j) = temp

              itemp = order(i)
              order(i) = order(j)
              order(j) = itemp
            else if (i == j) then
              i = i + 1
              exit
            else
              exit
            end if
          end do

          if (left_end < j) call quick_sort_1(left_end, j)
          if (i < right_end) call quick_sort_1(i, right_end)
        end if

        end subroutine quick_sort_1


        subroutine interchange_sort(left_end, right_end)

        integer, intent(in) :: left_end, right_end

        !     local variables
        integer             :: i, j, itemp
        real                :: temp

        do i = left_end, right_end - 1
          do j = i+1, right_end
            if (list(i) > list(j)) then
              temp = list(i)
              list(i) = list(j)
              list(j) = temp
              itemp = order(i)
              order(i) = order(j)
              order(j) = itemp
            end if
          end do
        end do

        end subroutine interchange_sort

      end subroutine quick_sort

! getvox.f 
! Generate voxels on LDPs
! Copyright (C) 2020 Jiahua He

! Codes for reading SITUS files (line 119-164) were taken from MAINMAST.f (http://kiharalab.org/mainmast/) under GNU General Public License version 3

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

      MODULE COMMONPATH
      integer,parameter,private::MAXNLDP=100000
      integer nldp
      real ldps(3,MAXNLDP)
      protected::nldp,ldps
      CONTAINS
      SUBROUTINE READ_PATH(fpath)

      implicit none
      character*100 fpath,line

      nldp=0

      open(100,file=fpath,status='old')
      do
        read(100,'(A80)',end=101) line
        if(line(1:4).ne."ATOM") cycle
        
        nldp=nldp+1
        read(line(31:38),*),ldps(1,nldp)
        read(line(39:46),*),ldps(2,nldp)
        read(line(47:54),*),ldps(3,nldp)
        
      end do
101   close(100)
      END SUBROUTINE READ_PATH
      END MODULE COMMONPATH

      MODULE COMMONMAP
      integer,parameter::msize=300
      real gstep,origin(3)
      integer ngrid(3)
      real map(msize,msize,msize)
      real fin(msize**3)
      integer posi(3,msize**3)
      real posr(3,msize**3)
      integer member(msize*msize*10)
      integer mapping(msize*msize*10)
      END MODULE COMMONMAP

      PROGRAM MAIN
      USE COMMONMAP
      USE COMMONPATH
      
      implicit none
      integer,parameter::vsize=4
      integer max_line,remain
      integer i,j,k,l,m,n
      real tmpr(3),dist2,mdist2,tmpr1
      integer xlower,ylower,zlower
      integer xupper,yupper,zupper
      integer x0,y0,z0
      integer x1,y1,z1
      integer dxl,dyl,dzl
      integer dxu,dyu,dzu
      character*100,fmap,fpath
      integer narg
      character*100,allocatable :: argv(:)
      real buffer(10)
      real,allocatable::voxel(:,:,:)
      real maxv,minv
      character*100 outfmt
      integer nout
      character*100 noutc
      logical tf

      narg=iargc()
      allocate( argv(0:narg) )
      
      do k=0,narg
        call getarg(k,argv(k))
      end do


      if(narg.ne.2)then
        write(*,*) "# GETVOX.F"
        write(*,*) "# Generate voxels on LDPs"
        write(*,*) " Usage: "//trim(argv(0))//
     *             " map.situs LDPs.pdb"
        stop
      end if

      fmap=argv(1)
      inquire(file=fmap,exist=tf)
      if(.not.tf) then 
        write(*,*) "input map file not exist!" 
        stop
      end if
      
      fpath=argv(2)
      inquire(file=fpath,exist=tf)
      if(.not.tf) then 
        write(*,*) "input pdb file not exist!" 
        stop
      end if
      call READ_PATH(fpath)

      open(1, file=fmap,status='old')
      read(1,*) gstep,(origin(k),k=1,3),(ngrid(k),k=1,3)
      if(maxval(ngrid).gt.msize)then
        stop '# Size of map is greater than default max size'
       end if
      max_line=ngrid(1)*ngrid(2)*ngrid(3)/10
      read(1,*)
      i=1;j=1;k=1
      do l=1,max_line
        read(1,*) (buffer(m),m=1,10)
        do m=1,10
          if(buffer(m).gt.0.0) then
            map(i,j,k)=buffer(m)
          end if
          i=i+1
          if(i.gt.ngrid(1))then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2))then
            j=1
            k=k+1
          end if
        end do
      end do

      remain=ngrid(1)*ngrid(2)*ngrid(3)-max_line*10

      if(remain.gt.0)then
        read(1,*),(buffer(m),m=1,remain)
        do m=1,remain
          if(buffer(m).gt.0.0) then
            map(i,j,k)=buffer(m)
          end if
          i=i+1
          if(i.gt.ngrid(1))then
            i=1
            j=j+1
          end if
          if(j.gt.ngrid(2))then
            j=1
            k=k+1
          end if
        end do
      end if
      close(1)

      nout=(vsize*2+2)**3
      write(noutc,'(i5)') nout
      outfmt="(A3,"//trim(adjustl(noutc))//"f6.3)"

      allocate(voxel(2*vsize+2,2*vsize+2,2*vsize+2))
      do i=1,nldp
        tmpr=(ldps(1:3,i)-origin)/gstep + 1
        voxel=0.0
        x0=floor(tmpr(1))
        y0=floor(tmpr(2))
        z0=floor(tmpr(3))

        x1=ceiling(tmpr(1))
        y1=ceiling(tmpr(2))
        z1=ceiling(tmpr(3))

        xlower=max(1,x0-vsize)
        ylower=max(1,y0-vsize)
        zlower=max(1,z0-vsize)
        xupper=min(ngrid(1),x1+vsize)
        yupper=min(ngrid(2),y1+vsize)
        zupper=min(ngrid(3),z1+vsize)
        
        dxl=xlower-(x0-vsize)+1
        dyl=ylower-(y0-vsize)+1
        dzl=zlower-(z0-vsize)+1

        dxu=(x1+vsize)-xupper
        dyu=(y1+vsize)-yupper
        dzu=(z1+vsize)-zupper

        voxel(dxl:(2*vsize+2-dxu)
     *       ,dyl:(2*vsize+2-dyu)
     *       ,dzl:(2*vsize+2-dzu))
     *= map(   xlower:xupper
     *        ,ylower:yupper
     *        ,zlower:zupper)
        maxv=maxval(voxel)
        minv=minval(voxel)
        voxel=(voxel-minv)/(maxv-minv)

        write (*,outfmt),'UNK',voxel

      end do

      deallocate(voxel)
      deallocate(argv)
      stop
      END PROGRAM MAIN
C******************************************************************************C
C******************************************************************************C

! align.f 
! Align sequence to main-chain paths with deep learning predicted information
! Copyright (C) 2020 Jiahua He

! This program is modified from ThreadCA.f (http://kiharalab.org/mainmast/)
! under GNU General Public License version 3

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

      module params
      parameter(mldp=4000,mres=1000,mpath=360,mmod=5760)
      character*3 aa3(21)
      character*1 aa1(21)
      real aascore(0:5,0:5)
      real ssscore(0:4,0:4)

      real wca,wss,waa,whelix,wsheet,wcoil,rmsd
      integer nout
      character*20 prefix
      logical imd

      data aa1 / 'G','A','S','C','V','T','I','P',
     *           'L','D','N','E','Q','M',
     *           'K','R',
     *           'H','F','Y','W',
     *           'X'/

      data aa3 /'GLY','ALA','SER','CYS','VAL','THR','ILE','PRO',
     *          'LEU','ASP','ASN','GLU','GLN','MET',
     *          'LYS','ARG',
     *          'HIS','PHE','TYR','TRP',
     *          'UNK'/

      data aascore / 0.7,-0.4,-0.8,-1.4,-0.1,-1e4,
     *              -0.4, 0.6,-0.5,-1.2,-0.1,-1e4,
     *              -0.8,-0.5, 1.1,-1.2,-0.1,-1e4,
     *              -1.4,-1.2,-1.2, 1.2,-0.1,-1e4,
     *              -0.1,-0.1,-0.1,-0.1,-0.1,-1e4,
     *              -1e4,-1e4,-1e4,-1e4,-1e4,-1e4/

      data ssscore / 1.5,-1.0,-0.5,-0.1,-1e4,
     *              -1.0, 1.0,-0.5,-0.1,-1e4,
     *              -0.5,-0.5, 0.5,-0.1,-1e4,
     *              -0.1,-0.1,-0.1,-0.1,-1e4,
     *              -1e4,-1e4,-1e4,-1e4,-1e4/
      end module params

      module commonspd3
      use params,only:mres
      character*100 fname_spd3
      integer iseqaa(mres),iseqaa0(mres),iseqss(mres)
      integer nres
      end module commonspd3

      module commonali
      use params,only:mldp,mmod,mres
      integer nmod
      integer alignment(mldp,mmod)
      real score(mmod),score_tmp(mmod),models(3,mres,mmod)
      dimension modelsaa(mres,mmod),modelsss(mres,mmod)
      end module commonali

      module commonpath
      use params,only:mldp,mpath
      implicit integer (i-z)
      character*100 fname_path

      integer npath
      dimension nldp(mpath)
      dimension fpathxyz(3,mldp,mpath),fpredca(mldp,mpath)
      dimension ipredaa(mldp,mpath),ipredss(mldp,mpath)
      
      end module commonpath

      program main
      use params
      use commonpath
      use commonspd3
      use commonali
      implicit integer (i-z)
      character*100 argv(0:30)
      dimension ftag(3,mres),fref(3,mres),frefx(3,mres)
      dimension listout(mmod)
      dimension fcom(3),ftrans(3),frot(3,3)
      dimension fxyz(3,mldp),fca(mldp)
      dimension iaa(mldp),iss(mldp)
      dimension mloc(1),order(mmod)
      character*5 str5
      logical existance
      
      write(*,*) "# ALIGN.F"
      write(*,*) "# Align sequence to main-chain paths"
      narg=iargc()
      do k=0,narg
        call getarg(k,argv(k))
      end do

      if(narg.lt.2) then
        write(*,*) " Usage: "//trim(argv(0))//
     *             " [path.pdb] [seq.spd3] (option)"
        write(*,*) " OPTIONS:"
        write(*,*) "---parameters for alignment----"
        write(*,*) " -wca    [f]: Weight of CA probability score"
        write(*,*) "              def = 1.6"
        write(*,*) " -waa    [f]: Weight of AA matching score"
        write(*,*) "              def = 1.0"
        write(*,*) " -wss    [f]: Weight of SS matching score"
        write(*,*) "              def = 0.5"
        write(*,*) " -whelix [f]: Weight of CA-CA distance for helix"
        write(*,*) "              def = 1.0"
        write(*,*) " -wsheet [f]: Weight of CA-CA distance for sheet"
        write(*,*) "              def = 0.7"
        write(*,*) " -wcoil  [f]: Weight of CA-CA distance for coil"
        write(*,*) "              def = 0.8"
        write(*,*) "---parameters for output----"
        write(*,*) " -rmsd   [f]: RMSD cut-off for clustering"
        write(*,*) "              def = 5.0"
        write(*,*) " -nout   [i]: Number of output top scored models"
        write(*,*) "              def = 10"
        write(*,*) " -prefix [s]: Prefix for output models"
        write(*,*) "              def = 'model'"
        write(*,*) " -imd       : Flag to write models
     * into individual .pdb files"
        stop
      end if

      wca=1.6
      whelix=1.0
      wsheet=0.7
      wcoil=0.8
      waa=1.0
      wss=0.5
      rmsd=5.0
      nout=10
      prefix="model"
      imd=.false.

      fname_path=argv(1)
      inquire(file=fname_path,exist=existance)
      if(.not.existance) then
        write(*,*) "# LDP path file "//trim(fname_path)//" not exist!"
        stop
      end if

      fname_spd3=argv(2)
      inquire(file=fname_spd3,exist=existance)
      if(.not.existance) then
        write(*,*) "# SPD3 file"//trim(fname_spd3)//" not exist!"
      end if

      k=3
      do while(k.le.narg)
        if(argv(k).eq.'-wca') then
          if(k+1.le.narg) read(argv(k+1),*) wca
          k=k+2
        else if(argv(k).eq.'-whelix') then
          if(k+1.le.narg) read(argv(k+1),*) whelix
          k=k+2
        else if(argv(k).eq.'-wsheet') then
          if(k+1.le.narg) read(argv(k+1),*) wsheet
          k=k+2
        else if(argv(k).eq.'-wcoil') then
          if(k+1.le.narg) read(argv(k+1),*) wcoil
          k=k+2
        else if(argv(k).eq.'-waa') then
          if(k+1.le.narg) read(argv(k+1),*) waa
          k=k+2
        else if(argv(k).eq.'-wss') then
          if(k+1.le.narg) read(argv(k+1),*) wss
          k=k+2
        else if(argv(k).eq.'-rmsd') then
          if(k+1.le.narg) read(argv(k+1),*) rmsd
          k=k+2
        else if(argv(k).eq.'-nout') then
          if(k+1.le.narg) read(argv(k+1),*) nout
          k=k+2
        else if(argv(k).eq.'-prefix') then
          if(k+1.le.narg) read(argv(k+1),*) prefix
          k=k+2
        else if(argv(k).eq.'-imd') then
          imd=.true.
          k=k+1
        else 
          write(*,*) "# Unknown parameter input"
          stop
        end if
      end do

      write(*,*) "# Parameter setting:"
      write(*,*) "# wca =",wca
      write(*,*) "# whelix =",whelix
      write(*,*) "# wsheet =",wsheet
      write(*,*) "# wcoil =",wcoil
      write(*,*) "# waa =",waa
      write(*,*) "# wss =",wss
      write(*,*) "# rmsd =",rmsd
      write(*,*) "# nout =",nout
      write(*,*) "# prefix = ",trim(prefix)
      
      write(*,*) "# Reading LDP paths from file "//trim(fname_path)
      call read_path

      write(*,*) "# Reading SS prediction from file "//trim(fname_spd3)
      call read_spd3

      write(*,*) "# Start alignment!"

      nmod=0
      do  i=1,npath
        sn=nldp(i)
        if (sn.gt.1) then
          do reverse=1,2
            do j=1,sn
              if (reverse.eq.1) then
                jt=j
              else 
                jt=sn-j+1
              end if
              fxyz(1:3,j)=fpathxyz(1:3,jt,i)
              fca(j)=fpredca(jt,i)
              iaa(j)=ipredaa(jt,i)
              iss(j)=ipredss(jt,i)
            end do

            call get_ali(sn,fxyz(:,1:sn),fca(1:sn),
     *                      iaa(1:sn),iss(1:sn))
          end do
        end if
      end do

      score_tmp(1:nmod)=score(1:nmod)
      call quick_sort(-score(1:nmod),order(1:nmod),nmod)
      score(1:nmod)=score_tmp(1:nmod)

      if(imd) then
        write(*,*) "# write models into individual .pdb files"

        nout0=0
        do i=1,nmod
          k=order(i)
          
          if(nout0.gt.0) then
            ftag(1:3,1:nres)=models(1:3,1:nres,k)
            do j=1,nout0
              fref(1:3,1:nres)=models(1:3,1:nres,listout(j))
              fcom(:)=0.0
              frmsd=0.0
              do l=1,nres
                fcom(:)=fcom(:)+fref(1:3,l)/nres
              end do
              call orient_rot(nres,ftag,fref,fcom,frot,ftrans)
              call transform(nres,fref,fcom,frot,ftrans,frefx)
              frmsd=0.0
              do l=1,nres
                dx=frefx(1,l)-ftag(1,l)
                dy=frefx(2,l)-ftag(2,l)
                dz=frefx(3,l)-ftag(3,l)
                frmsd=frmsd+dx**2+dy**2+dz**2
              end do
              frmsd=sqrt(frmsd/nres)
              
              if(frmsd.lt.rmsd) then
                goto 10
              end if

            end do
          end if 

          nout0=nout0+1
          listout(nout0)=k
         
          write(str5,'(i5)') nout0
          open(300,file=trim(prefix)//'.'//
     *                  trim(adjustl(str5))//".pdb",status="unknown")

          write(300,'("MODEL",i5)') nout0
          write(300,'("# Score",x,f8.3)') score(k)
          do j=1,nres
            write(300,'(A4,i7,A6,A3,A2,I4,A4,3f8.3,2i6)')
     *                "ATOM",j,"  CA  ",aa3(iseqaa0(j))," A",
     *                j,'    ',
     *                models(1:3,j,k),modelsaa(j,k),modelsss(j,k)
          end do
          write(300,'("ENDMDL")')

          close(300)
          
          if(nout0.ge.nout) exit

10      end do

      else
        write(*,*) "# write models into one .pdb file"
        open(300,file=trim(prefix)//".pdb",status="unknown")

        nout0=0
        do i=1,nmod
          k=order(i)
          
          if(nout0.gt.0) then
            ftag(1:3,1:nres)=models(1:3,1:nres,k)
            do j=1,nout0
              fref(1:3,1:nres)=models(1:3,1:nres,listout(j))
              fcom(:)=0.0
              frmsd=0.0
              do l=1,nres
                fcom(:)=fcom(:)+fref(1:3,l)/nres
              end do
              call orient_rot(nres,ftag,fref,fcom,frot,ftrans)
              call transform(nres,fref,fcom,frot,ftrans,frefx)
              frmsd=0.0
              do l=1,nres
                dx=frefx(1,l)-ftag(1,l)
                dy=frefx(2,l)-ftag(2,l)
                dz=frefx(3,l)-ftag(3,l)
                frmsd=frmsd+dx**2+dy**2+dz**2
              end do
              frmsd=sqrt(frmsd/nres)
              
              if(frmsd.lt.rmsd) then
                goto 20
              end if

            end do
          end if 

          nout0=nout0+1
          listout(nout0)=k
         

          write(300,'("MODEL",i5)') nout0
          write(300,'("# Score",x,f8.3)') score(k)
          do j=1,nres
            write(300,'(A4,i7,A6,A3,A2,I4,A4,3f8.3,2i6)')
     *                "ATOM",j,"  CA  ",aa3(iseqaa0(j))," A",
     *                j,'    ',
     *                models(1:3,j,k),modelsaa(j,k),modelsss(j,k)
          end do
          write(300,'("ENDMDL")')

          if(nout0.ge.nout) exit

20      end do
        close(300)
      end if

      stop 
      end program main

      subroutine read_path
      use commonpath
      implicit integer (i-z)
      character*100 line

      open(100,file=fname_path,status='old') 
      npath=0
      nldp(:)=0
      sn=0
      do
        read(100,'(a)',end=101) line
        if(line(1:5).eq."MODEL") then 
          npath=npath+1
        else if(line(1:6).eq."ENDMDL") then
          write(*,*) "# Path ",npath," consists of ",sn," LDPs."
          nldp(npath)=sn
          sn=0
        else if(line(1:4).eq."ATOM") then
          sn=sn+1
          read(line(31:72),'(3f8.3,f6.2,2i6)')
     *         (fpathxyz(k,sn,npath),k=1,3),fpredca(sn,npath),
     *         ipredaa(sn,npath),ipredss(sn,npath)
        end if
      end do

101   close(100)

      write(*,*) "# Number of paths: ",npath
      
      return
      end subroutine read_path

      subroutine read_spd3
      use params,only:aa1
      use commonspd3
      implicit integer (i-z)
      character*100 line
      character aatype,sstype

      iseqaa0(:)=21
      iseqaa(:)=4
      iseqss(:)=3
      nres=0
      write(*,*) "# Reading predicted SS from .SPD3 file "
     *//trim(fname_spd3)
      open(200,file=fname_spd3,status='old')
      do 
        read(200,'(a)',end=201) line
        if(line(1:1).eq.'#') cycle
        read(line,*),nres,aatype,sstype

        j=0
        do i=1,21
          if(aatype.eq.aa1(i)) then 
            j=i
            exit
          end if
        end do

        if (j.eq.0) then 
          write(*,*) "# Skip unknown AA type: ",aatype
        end if

        iseqaa0(nres)=j

        if (1<=j.and.j<=8) then 
          iseqaa(nres)=0
        else if(9<=j.and.j<=14) then
          iseqaa(nres)=1
        else if(15<=j.and.j<=16) then
          iseqaa(nres)=2
        else if(17<=j.and.j<=20) then
          iseqaa(nres)=3
        else if(j.eq.21) then
          iseqaa(nres)=4
        end if
        
        if (sstype.eq.'H') then
          iseqss(nres)=0
        else if (sstype.eq.'E') then
          iseqss(nres)=1
        else if (sstype.eq.'C') then
          iseqss(nres)=2
        else 
          write(*,*) "# Skip unknown SS type: ",sstype
        end if
      end do
201   close(200)
      
      write(*,*) "# Number of residues: ",nres
      return 
      end subroutine read_spd3

      subroutine
     *get_ali(sn,fxyz,fca,iaa,iss)
      use commonspd3,only:nres,iseqaa,iseqss
      use commonali
      use params
      implicit integer (i-z)

      dimension iaa(sn),iss(sn),fxyz(3,sn),fca(sn)

      dimension ftmp(3)
      integer direc(0:sn,0:nres),prev(0:sn,0:nres)
      real dmat(0:sn,0:sn),pmat(0:sn,0:nres),smat(0:sn,0:nres)
      real diag,up,left

      gap=-10000.0
        
      dmat(:,:)=0.0
      do i=1,sn
        do j=i+1,sn
          ftmp=fxyz(:,i)-fxyz(:,j)
          dmat(i,j)=sqrt(dot_product(ftmp,ftmp))
          dmat(j,i)=dmat(i,j)
        end do
      end do

      do i=1,sn
        do j=i+2,sn
          dij=dmat(i,j)
          do k=i+1,j-1
            if(dmat(i,k).gt.dij) dij=dmat(i,k)
            if(dmat(j,k).gt.dij) dij=dmat(j,k)
          end do
          dmat(i,j)=dij
          dmat(j,i)=dmat(i,j)
        end do
      end do
      
      pmat(:,:)=0.0
      do m=1,sn
        do n=1,nres
          pmat(m,n)=aascore(iaa(m),iseqaa(n))*waa
     *             +ssscore(iss(m),iseqss(n))*wss
     *             +fca(m)*wca
        end do
      end do


c      do reverse=0,1
        do dcaca=3.1,3.8,0.1
      
          direc(:,:)=0 ! 0 null, 1 up, 2 left, 3 diag
          prev(:,:)=0 ! previous aligned seq id, 0 null
          smat(:,:)=0.0 ! S-W scoring matrix
      
          do m=1,sn
            smat(m,0)=0.0
            direc(m,0)=2
            prev(m,0)=0
          end do

          do n=1,nres
            smat(0,n)=smat(0,n-1)+gap
            direc(0,n)=1
            prev(0,n)=0
          end do        

          do m=1,sn
            do n=1,nres
              up=smat(m,n-1)+gap
              left=smat(m-1,n)
        
              refnow=m-1
              tagnow=n-1
              refnow=prev(refnow,tagnow)
              if(direc(refnow,tagnow).eq.0) then 
                diag=smat(m-1,n-1)+pmat(m,n)
              else
                dist=dmat(m,refnow)
                if(iseqss(n).eq.0) then
                  diag=smat(m-1,n-1)+pmat(m,n)-abs(dist-dcaca)*whelix
                else if(iseqss(n).eq.1) then
                  diag=smat(m-1,n-1)+pmat(m,n)-abs(dist-dcaca)*wsheet
                else 
                  diag=smat(m-1,n-1)+pmat(m,n)-abs(dist-dcaca)*wcoil
                end if
              end if

              if(diag.ge.left.and.diag.ge.up) then
                smat(m,n)=diag
                direc(m,n)=3
                prev(m,n)=m
              else if (up.ge.left) then
                smat(m,n)=up
                direc(m,n)=1
                prev(m,n)=prev(m,n-1)
              else 
                smat(m,n)=left
                direc(m,n)=2
                prev(m,n)=prev(m-1,n)
              end if
            end do
          end do
       
          nmod=nmod+1
          alignment(:,nmod)=0
          refnow=sn
          tagnow=nres
          do
            if(direc(refnow,tagnow).eq.0) exit
            if(direc(refnow,tagnow).eq.3) then
              alignment(tagnow,nmod)=refnow
              refnow=refnow-1
              tagnow=tagnow-1
            else if(direc(refnow,tagnow).eq.1) then
              tagnow=tagnow-1
            else if(direc(refnow,tagnow).eq.2) then
              refnow=refnow-1
            else
              write(*,*) "# Unknown S-W alignment direction"
              stop
            end if
          end do
 
          score(nmod)=smat(sn,nres)
          do j=1,nres
            k=alignment(j,nmod)
            models(1:3,j,nmod)=fxyz(1:3,k)
            modelsaa(j,nmod)=iaa(k)
            modelsss(j,nmod)=iss(k)
          end do
          write(*,'(" # Model",i5,x,f8.3)') nmod,score(nmod)

        end do ! end loop dcaca
c      end do ! end loop reverse

      return
      end subroutine get_ali
        
      recursive subroutine quick_sort(list,order,n)
      use params,only:mmod
      implicit none
      real  list(mmod)
      integer order(mmod)

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

      subroutine orient_rot (n,x0,x1,com_x1,rot,trans)

**    This subroutine is to generate the best rotational matrix using
**    Kabsch et al's algorithm (1978). The subroutine calls for another subroutine,
**    'cjcb', which uses the Jacobi algorithm to calculate the eigenvectors of
**    the correlation matrix.

      parameter(eps=0.00001)

      DOUBLE PRECISION R,RR,V,A,B,RA
      integer sflag

      dimension V(3,3)

      dimension x0(3,n),x1(3,n)
      dimension com0(3),com1(3),com_x1(3),trans(3)
      dimension R(3,3),RR(3,3),rot(3,3)
      dimension a(3,3),b(3,3),Ra(3,3),ip(3)

      com0=0.0
      com1=0.0

      do k=1,3
        do i=1,n
          com0(k)=com0(k)+x0(k,i)
          com1(k)=com1(k)+x1(k,i)
        end do
      end do

      an=1.0/n
      com0=com0*an
      com1=com1*an

      do i=1,3
        do j=1,3
          R(i,j)=0.0
          do k=1,n
            R(i,j)=R(i,j)+( x0(i,k)-com0(i) )*( x1(j,k)-com1(j) )
          end do
        end do
      end do

      do i=1,3
        do j=i,3
          RR(i,j)=0.0
          do k=1,3
            RR(i,j) = RR(i,j) + R(k,i)*R(k,j)
          end do
          RR(j,i)=RR(i,j)
        end do
      end do

      call cjcb(rr,3,eps,v,sflag)

      if(sflag.eq.0)then

        rot=0.0
        do i=1,3
          rot(i,i)=1
        end do

        call get_translation(com0,com1,com_x1,rot,trans,0)

        return

      end if

      ip=0
      do k1=1,3
        xmax=-100000.0
        do k2=1,3
         if(xmax < RR(k2,k2).and.ip(k2)==0)then
           xmax=RR(k2,k2)
           kmax=k2
         end if
        end do
        ip(kmax)=1
        a(k1,:)=v(:,kmax)
      end do

      temp=sqrt(a(1,1)*a(1,1)+a(1,2)*a(1,2)+a(1,3)*a(1,3))
      a(1,:)=a(1,:)/temp

      temp=sqrt(a(2,1)*a(2,1)+a(2,2)*a(2,2)+a(2,3)*a(2,3))
      a(2,:)=a(2,:)/temp

      a(3,1)= a(1,2)*a(2,3)-a(1,3)*a(2,2)
      a(3,2)= a(1,3)*a(2,1)-a(1,1)*a(2,3)
      a(3,3)= a(1,1)*a(2,2)-a(1,2)*a(2,1)

      b(1,1)=R(1,1)*a(1,1)+R(1,2)*a(1,2)+R(1,3)*a(1,3)
      b(1,2)=R(2,1)*a(1,1)+R(2,2)*a(1,2)+R(2,3)*a(1,3)
      b(1,3)=R(3,1)*a(1,1)+R(3,2)*a(1,2)+R(3,3)*a(1,3)

      temp=sqrt(b(1,1)*b(1,1)+b(1,2)*b(1,2)+b(1,3)*b(1,3))
      b(1,:)=b(1,:)/temp

      b(2,1)=R(1,1)*a(2,1)+R(1,2)*a(2,2)+R(1,3)*a(2,3)
      b(2,2)=R(2,1)*a(2,1)+R(2,2)*a(2,2)+R(2,3)*a(2,3)
      b(2,3)=R(3,1)*a(2,1)+R(3,2)*a(2,2)+R(3,3)*a(2,3)

      temp=sqrt(b(2,1)*b(2,1)+b(2,2)*b(2,2)+b(2,3)*b(2,3))
      b(2,:)=b(2,:)/temp

      b(3,1)= b(1,2)*b(2,3)-b(1,3)*b(2,2)
      b(3,2)= b(1,3)*b(2,1)-b(1,1)*b(2,3)
      b(3,3)= b(1,1)*b(2,2)-b(1,2)*b(2,1)

      do i=1,3
        do j=1,3
          rot(i,j)=0.0
          do k=1,3
            rot(i,j)=rot(i,j)+b(k,i)*a(k,j)
          end do
        end do
      end do

      call get_translation(com0,com1,com_x1,rot,trans,1)

      return
      end subroutine orient_rot


      subroutine get_translation(com0,com1,com_x1,rot,trans,kflag)

**  This subroutine is to calculate the translational vector between
**  the two match sets according to their centers of mass

      dimension com0(3),com1(3),com_x1(3),trans(3)
      dimension rot(3,3)

      if(kflag.eq.1)then

        do i=1,3
          trans(i) = rot(i,1) * (com_x1(1) - com1(1))
     *             + rot(i,2) * (com_x1(2) - com1(2))
     *             + rot(i,3) * (com_x1(3) - com1(3))
     *             + com0(i) - com_x1(i)
        end do

      else

        trans = com0 - com_x1

      end if

      return
      end subroutine get_translation

      subroutine transform (natm,old,coml,rot,trans,new)

**    This subroutine is to transform a set of old coordinates to
**    new coordinates by the given transform matrix.

      real old(3,natm), new(3,natm)
      dimension coml(3),rot(3,3),trans(3)

      do  i=1,natm
        do k=1,3
          new(k,i) = rot(k,1) * (old(1,i) - coml(1))
     *             + rot(k,2) * (old(2,i) - coml(2))
     *             + rot(k,3) * (old(3,i) - coml(3))
     *             + coml(k) + trans(k)
        end do
      end do

      return
      end subroutine transform

      subroutine cjcb(a,n,eps,v,l)

**    Use the Jacobi algorithm to calculate the eigenvalues and eigenvectors
**    for a real symmetric maxtrix.

      dimension a(n,n),v(n,n)
      double precision a,v,fm,cn,cn2,sn,sn2,omega,x,y
      integer p,q
      l=1

      v=0.0
      do i=1,n
        v(i,i)=1.0
      end do

30    fm=0.0

      do i=2,n
      do j=1,i-1
        if (abs(a(i,j)).gt.fm) then
          fm=abs(a(i,j))
          p=i
          q=j
        end if
      end do
      end do

      if (fm.lt.eps) then
        l=1
        return
      end if

      if (l.gt.100) then
        l=0
        return
      end if

      l=l+1
      x=-a(p,q)
      y=(a(q,q)-a(p,p))/2.0
      omega=x/sqrt(x*x+y*y)
      if (y.lt.0.0) omega=-omega
      sn=1.0+sqrt(1.0-omega*omega)
      sn=omega/sqrt(2.0*sn)
      sn2=sn*sn
      cn=sqrt(1.0-sn2)
      cn2=cn*cn
      fm=a(p,p)
      a(p,p)=fm*cn2+a(q,q)*sn2+a(p,q)*omega
      a(q,q)=fm*sn2+a(q,q)*cn2-a(p,q)*omega
      a(p,q)=0.0
      a(q,p)=0.0

      do j=1,n
        if ((j.ne.p).and.(j.ne.q)) then
          fm=a(p,j)
          a(p,j)=fm*cn+a(q,j)*sn
          a(q,j)=-fm*sn+a(q,j)*cn
        end if
      end do

      do i=1,n
        if ((i.ne.p).and.(i.ne.q)) then
          fm=a(i,p)
          a(i,p)=fm*cn+a(i,q)*sn
          a(i,q)=-fm*sn+a(i,q)*cn
        end if
      end do

      do i=1,n
        fm=v(i,p)
        v(i,p)=fm*cn+v(i,q)*sn
        v(i,q)=-fm*sn+v(i,q)*cn
      end do

      goto 30

      return
      end subroutine cjcb

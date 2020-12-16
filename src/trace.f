! trace.f 
! Trace the main-chain path with MAINMAST algorithm
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

      module params
      parameter(msize=300)
      parameter(mlink=10)
      character*5 ss(5)
      character*1 aa1(21)
      character*3 aa3(21)

      integer nrd,nnb,mtabu,ntraj
      real dkeep,constraint,rlocal
      end module params

      module commonldp
      use params,only:msize,mlink
      implicit integer (i-z)
      character*100 fname_ldp

      integer sn,nedge
      dimension fxyz(3,msize**2),fca(msize**2)
      dimension taa(msize**2),tss(msize**2)
      dimension edge(msize**2*mlink),ide(2,msize**2*mlink),
     *          edge_dens(msize**2*mlink)
      end module commonldp

      module commongraph
      use params,only:msize,mlink
      integer nlink(msize**2),
     *        link(msize**2,mlink),
     *        idelink(msize**2,mlink)
      end module commongraph

      program main
      use params
      use commonldp
      use commongraph
      implicit integer (i-z)
      character*100 argv(0:30)
      dimension path_best(100,msize**2),nv_best(100)
      dimension fpathxyz(3,msize**2)
      dimension ipredaa(msize**2),ipredss(msize**2),fpredca(msize**2)
      logical existance
      
      write(*,*) "# TRACE.F"
      write(*,*) "# Trace the main-chain path with MAINMAST"
      narg=iargc()
      do k=0,narg
        call getarg(k,argv(k))
      end do

      if(narg.lt.1) then
        write(*,*) " Usage: "//trim(argv(0))//
     *             " [LDP.pdb] (option)"
        write(*,*) " OPTIONS:"
        write(*,*) " -dkeep  [f]: keep edge shorter than"
        write(*,*) "              def=0.5"
        write(*,*) " -ntraj  [i]: number of trajectories"
        write(*,*) "              def=10"
        write(*,*) " -nrd    [i]: number of iterations"
        write(*,*) "              def=1000"
        write(*,*) " -nnb    [i]: number of neighborss"
        write(*,*) "              def=30"
        write(*,*) " -mtabu  [i]: maximum size of tabu-list"
        write(*,*) "              def=100"
        write(*,*) " -const  [f]: constraint of total length of edge"
        write(*,*) "              def=1.01"
        write(*,*) " -rlocal [f]: radius of locat mst"
        write(*,*) "              def=10.0"
        stop
      end if

      fname_ldp=argv(1)
      inquire(file=fname_ldp,exist=existance)
      if(.not.existance) then
        write(*,*) "# LDP.pdb file "//trim(fname_ldp)//" not exist!"
        stop
      end if

      dkeep=0.5
      ntraj=10
      nrd=1000
      nnb=30
      mtabu=100
      constraint=1.01
      rlocal=10.0

      wcaca=0.4
      wca=1.4
      waa=1.0
      wss=0.5

      k=2
      do while(k.le.narg)
        if(argv(k).eq.'-dkeep') then
          if(k+1.le.narg) read(argv(k+1),*) dkeep
          k=k+2
        else if(argv(k).eq.'-ntraj') then
          if(k+1.le.narg) read(argv(k+1),*) ntraj
          k=k+2
        else if(argv(k).eq.'-nrd') then
          if(k+1.le.narg) read(argv(k+1),*) nrd
          k=k+2
        else if(argv(k).eq.'-nnb') then
          if(k+1.le.narg) read(argv(k+1),*) nnb
          k=k+2
        else if(argv(k).eq.'-mtabu') then
          if(k+1.le.narg) read(argv(k+1),*) mtabu
          k=k+2
        else if(argv(k).eq.'-const') then
          if(k+1.le.narg) read(argv(k+1),*) constraint
          k=k+2
        else if(argv(k).eq.'-rlocal') then
          if(k+1.le.narg) read(argv(k+1),*) rlocal
          k=k+2
        else 
          write(*,*) "# Unknown parameter input"
          stop
        end if
      end do

      write(*,*) "# Keep edge shorter than=",dkeep
      write(*,*) "# Number of search trajectories=",ntraj
      write(*,*) "# Number of search rounds=",nrd
      write(*,*) "# Number of searcg neighbors=",nnb
      write(*,*) "# Maximum size of tabu-list=",mtabu
      write(*,*) "# Constraint of total edge length=",constraint
      write(*,*) "# Radius of local MST=",rlocal

      write(*,*) "# Reading LDPs from file "//trim(fname_ldp)
      call read_ldp

      call kruskal(path_best,nv_best)

      do i=1,ntraj
        if(nv_best(i).gt.1) then

          write(*,'("MODEL",i4)') i
          do j=1,nv_best(i)
            k=path_best(i,j)
            write(*,'("ATOM  ",i5,"  CA  ALA A",
     *                 i4,"    ",3f8.3,f6.2,2i6)'),
     *                 j,j,fxyz(1:3,k),fca(k),taa(k),tss(k)
          end do
          write(*,'("ENDMDL")')

        end if
      end do
       
      stop 
      end program main

      subroutine read_ldp
      use commonldp
      implicit integer (i-z)
      character*100 line

      sn=0
      nedge=0
      open(100,file=fname_ldp,status='old') 
      do
        read(100,'(a)',end=101) line
        if(line(1:4).eq."ATOM") then
          sn=sn+1
          read(line(31:72),'(3f8.3,f6.2,2i6)')
     *         (fxyz(k,sn),k=1,3),fca(sn),taa(sn),tss(sn)
        else if(line(1:4).eq."BOND") then
          nedge=nedge+1
          read(line(6:41),'(2i6,2f12.5)')
     *         (ide(k,nedge),k=1,2),edge(nedge),edge_dens(nedge)
        end if
      end do

101   close(100)
      
      write(*,*) "# Number of LDPs=",sn
      write(*,*) "# Number of edges=",nedge
      return
      end subroutine read_ldp

      subroutine kruskal(path_best,nv_best)
      use params
      use commonldp
      use commongraph
      implicit integer(i-z)
      dimension nlink_tmp(msize**2),
     *          link_tmp(msize**2,mlink),
     *          idelink_tmp(msize**2,mlink)
      dimension ftmp(3)
      logical flag
      integer cid(sn)
      dimension order(sn*sn),tree(sn*sn),edge_tmp(sn*sn),
     *  tree_traj(sn*sn),tree_best(sn*sn),tree_rd(sn*sn),
     *  tree_mst(sn*sn)
      logical activeg(msize**2),activeg_traj(msize**2)
      real score_best,score_rd,score_traj,score,score_traj_best,
     *     traj_score_best
      real length_rd,length_traj,length,length_mst,length_allow
      dimension ope(2),tabulist(2,msize),list_edge(msize**2)
      logical is_local(msize**2),usedv(sn),mstv(sn)
      dimension pathv(sn),path_best(100,msize**2),nv_best(100)
      dimension path_tmp(msize**2)
      dimension path_traj(msize**2),path_rd(msize**2)

      edge_tmp(1:nedge)=edge(1:nedge)
      call quick_sort(edge(1:nedge),order(1:nedge),nedge)
      edge(1:nedge)=edge_tmp(1:nedge)


      write(*,*) "# Sum of edge density=",sum(edge_dens(1:nedge))
      write(*,*) "# Avg of edge density=",sum(edge_dens(1:nedge))/nedge

      ! build MST again

      mstv(:)=.false.
       
      do i=1,sn
        cid(i)=i
      end do

      ntree=0
      activeg=.false.
      do i=1,nedge
        me=order(i)
        v1=ide(1,me)
        v2=ide(2,me)
        mstv(v1)=.true.
        mstv(v2)=.true.
        if(cid(v1).ne.cid(v2)) then
          ntree=ntree+1
          tree(ntree)=me
          activeg(me)=.true. ! using edge in initial mst
          tmpid=cid(v2)
          do j=1,sn
            if(cid(j).eq.tmpid) cid(j)=cid(v1)
          end do
        end if
      end do

      ! calculate local MST
      is_local(1:nedge)=.false.
      dist2=rlocal**2
      do pos=1,sn

        if(.not.mstv(pos)) cycle
        do i=1,sn
          cid(i)=i
        end do

        do i=1,nedge
          me=order(i) ! min edge
          v1=ide(1,me)
          v2=ide(2,me)
          if(cid(v1).eq.cid(v2)) cycle
          ftmp=fxyz(1:3,pos)-fxyz(1:3,v1)
          if(dot_product(ftmp,ftmp).gt.dist2) cycle
          ftmp=fxyz(1:3,pos)-fxyz(1:3,v2)
          if(dot_product(ftmp,ftmp).gt.dist2) cycle
          is_local(me)=.true.

          tmpid=cid(v2)
          do j=1,sn
            if(cid(j).eq.tmpid) cid(j)=cid(v1)
          end do

        end do
      end do

      nledge=0
      do i=1,nedge
        if(activeg(i)) is_local(i)=.true. ! add initial MST edge
        if(is_local(i)) nledge=nledge+1
      end do

      write(*,*) "# Number of local edges=",nledge

      ! calculated longest path of MST
      length=0.0
      do i=1,ntree
        length=length+edge(tree(i))
      end do
      length_mst=length
      call setup_tree(sn,ide,tree,ntree)
      call quality_tree(sn,edge_dens,
     *    ide,tree,ntree,score,usedv,pathv,nv)
      score_best=score

      write(*,*),"# MST",score,ntree,length_mst
      !--------------------------------!
      !!!         TABU search        !!!
      !--------------------------------!
      length_allow=constraint*length_mst
      tree_mst(1:ntree)=tree(1:ntree)
      nv_best(1:ntraj)=0

      do traj=1,ntraj ! loop trajectories
        traj_score_best=0.0
        ntabu=0
        tree_traj(1:ntree)=tree_mst(1:ntree)
        activeg_traj(1:nedge)=activeg(1:nedge)
        length_traj=length_mst

        maxtry=nledge
        maxtry2=1

        do rd=2,nrd ! loop rounds

          call setup_tree(sn,ide,tree_traj(1:ntree),ntree)
          tree_rd(1:ntree)=tree_traj(1:ntree)

          nlink_tmp(1:sn)=nlink(1:sn)
          do j=1,sn
            link_tmp(j,1:nlink(j))=link(j,1:nlink(j))
            idelink_tmp(j,1:nlink(j))=idelink(j,1:nlink(j))
          end do

          ntry2=0

          score_rd=0

          do nb=1,nnb ! loop neighbors

            nlink(1:sn)=nlink_tmp(1:sn)
            do j=1,sn
              link(j,1:nlink(j))=link_tmp(j,1:nlink(j))
              idelink(j,1:nlink(j))=idelink_tmp(j,1:nlink(j))
            end do
c            call check_tree(sn)

            ntry=0

3001        call random_number(a)

            ntry=ntry+1
            if(ntry.gt.maxtry) then
              write(*,*),"# Too many tries, clear the tabu list"
              ntabu=0
              ntry=0
              ntry2=ntry2+1
              !avoid too many loops
              if(ntry2.gt.maxtry2) then
                ntry2=0
                write(*,*),"# Too many tries, end current trajectory"
                goto 3003
              end if
            end if
        
            orderid=int(a*ntree)+1
            p11=tree_traj(orderid)
            ! keep short edge in MST
            if(activeg(p11).and.edge(p11).lt.dkeep) goto 3001

            call rmpath_tree(sn,ide,p11)
            ! split the tree into two sub-trees from the removed edge
            call split_chain(cid,sn,ide,ntree,p11)

            ! set the list of possible edges to add
            nee=0
            do j=1,nedge
              if(.not.is_local(j)) cycle
              ! not used by the tree
              if(activeg_traj(j)) cycle

              v1=ide(1,j);v2=ide(2,j)
              if(cid(v1).eq.cid(v2)) cycle
              ! not in the same sub-tree
              if(cid(v1)*cid(v2).eq.0) cycle
              ! restrain the total length of tree
              if(length_traj+edge(j)-edge(p11).
     *        gt.length_allow) cycle

              nee=nee+1
              list_edge(nee)=j
            end do

            ! if no edge can be added, choose another edge to remove
            if(nee.eq.0) then 
              call addpath_tree(sn,ide,p11)
              goto 3001
            end if

            flag=.false.

            do j=1,nee
              call random_number(a)
              p12=list_edge(int(a*nee)+1)
              length=length_traj+edge(p12)-edge(p11)

             ! tabu check
              do k=1,ntabu
                if((tabulist(1,k).eq.p11.and.tabulist(2,k).eq.p12).or.
     *             (tabulist(1,k).eq.p12.and.tabulist(2,k).eq.p11))
     *          goto 3002
              end do
              tree(1:ntree)=tree_traj(1:ntree)
              tree(orderid)=p12
              flag=.true.
              exit

3002        end do

            if(.not.flag) then
              call addpath_tree(sn,ide,p11)
              goto 3001
            end if

            call addpath_tree(sn,ide,p12)

            call quality_tree(sn,edge_dens,ide,tree,ntree,
     *             score,usedv,pathv,nv)
        
            ! vertices of removed MST edge should not in the longest path
            if(usedv(ide(1,p11)).and.usedv(ide(2,p11))
     *         .and.activeg(p11)) then
              call addpath_tree(sn,ide,p11)
              call rmpath_tree(sn,ide,p12)
              goto 3001
            end if

            ! keep the trival change that reduce total length
            if(score.eq.score_traj.and.length.gt.length_traj) then
              call addpath_tree(sn,ide,p11)
              call rmpath_tree(sn,ide,p12)
              goto 3001
            end if
               
            if(score_rd.le.score) then 

              ! keep the change with same score but shorter total length
c              if(length.ge.length_rd.and.score.eq.score_rd) then
              if(score.eq.score_rd.and.length.ge.length_rd) then
                call addpath_tree(sn,ide,p11)
                call rmpath_tree(sn,ide,p12)
                goto 3001
              end if

              ! check kendall tau
              if(traj.gt.1) then
                path_tmp(1:nv_traj)=pathv(1:nv_traj)
     *                             -path_traj(1:nv_traj)
                if(maxval(path_tmp(1:nv_traj)).ne.0.and.
     *             minval(path_tmp(1:nv_traj)).ne.0) then
                  do p=1,traj-1
                    fken=fkendall(pathv(1:nv)   ,nv        ,
     *                            path_best(p,:),nv_best(p),sn)
                    if(abs(fken).gt.0.95) then
                      call addpath_tree(sn,ide,p11)
                      call rmpath_tree(sn,ide,p12)
                      goto 3001
                    end if
                  end do
                end if
              end if
              ! end check kendall tau

              ! update round best state
              score_rd=score
              length_rd=length
              tree_rd(1:ntree)=tree(1:ntree)
              path_rd(1:nv)=pathv(1:nv)
              nv_rd=nv
              ope(1:2)=(/p11,p12/)
              if(traj_score_best.lt.score) then
                traj_score_best=score
                path_best(traj,1:nv)=pathv(1:nv)
                nv_best(traj)=nv
              end if
            end if
          end do ! end neighbor

          write(*,*),traj,rd,score_rd,score_best,nv_best(traj)

          ! update traj state
          tree_traj(1:ntree)=tree_rd(1:ntree)
          path_traj(1:nv_rd)=path_rd(1:nv_rd)
          nv_traj=nv_rd
          score_traj=score_rd

          length_traj=0.0
          do i=1,ntree
            length_traj=length_traj+edge(tree_traj(i))
          end do

         ! update best state, update tabu list
          if(score_best.lt.score_rd) then ! reset tabu list
            score_best=score_rd
            tabulist(1,1)=ope(2)
            tabulist(2,1)=ope(1)
            ntabu=1
          else
            ntabu=ntabu+1
            if(ntabu.gt.mtabu) ntabu=mtabu
            tabulist(1:2,2:ntabu)=tabulist(1:2,1:ntabu-1)
            tabulist(1,1)=ope(2)
            tabulist(2,1)=ope(1)
          end if  
        
c          activeg_traj(1:nedge)=.false.
c          do i=1,ntree
c            activeg_traj(tree_traj(i))=.true.
c          end do

          activeg_traj(ope(1))=.false.
          activeg_traj(ope(2))=.true.

        end do ! end round
3003  end do ! end traj

      return
      end subroutine kruskal

      subroutine setup_tree(n,ide,tree,nt)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      dimension ide(2,n*n),tree(n*n)

      nlink=0
      link=0

      do i=1,nt
        it=tree(i)
        id1=ide(1,it)
        id2=ide(2,it)
        !input 
        nlink(id1)=nlink(id1)+1
        link(id1,nlink(id1))=id2
        idelink(id1,nlink(id1))=it

        nlink(id2)=nlink(id2)+1
        link(id2,nlink(id2))=id1
        idelink(id2,nlink(id2))=it
      end do

      return
      end subroutine setup_tree
        
      subroutine addpath_tree(n,ide,iadd)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      dimension ide(2,n*n)

      id1=ide(1,iadd)
      id2=ide(2,iadd)

      !input 
      nlink(id1)=nlink(id1)+1
      link(id1,nlink(id1))=id2
      idelink(id1,nlink(id1))=iadd

      nlink(id2)=nlink(id2)+1
      link(id2,nlink(id2))=id1
      idelink(id2,nlink(id2))=iadd

      return
      end subroutine addpath_tree

      subroutine rmpath_tree(n,ide,irm)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      dimension ide(2,n*n)

      id1=ide(1,irm)
      id2=ide(2,irm)

      !remove
      do i=1,nlink(id1)
        if(link(id1,i).eq.id2) then
          !shift
          link(id1,i:nlink(id1)-1)=link(id1,i+1:nlink(id1))
          idelink(id1,i:nlink(id1)-1)=idelink(id1,i+1:nlink(id1))
          nlink(id1)=nlink(id1)-1
          exit
        end if
      end do
      do i=1,nlink(id2)
        if(link(id2,i).eq.id1) then
          !shift
          link(id2,i:nlink(id2)-1)=link(id2,i+1:nlink(id2))
          idelink(id2,i:nlink(id2)-1)=idelink(id2,i+1:nlink(id2))
          nlink(id2)=nlink(id2)-1
          exit
        end if
      end do

      return
      end subroutine rmpath_tree

      subroutine split_chain(cid,n,ide,nt,igp)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      logical visited(msize*msize)
      dimension list_act(msize*msize) !active list
      logical flag
      integer cid(n)
      dimension ide(2,n*n)

      !Depth First Search
      !Split the tree from the removed edge into chain 1 and 2
      cid(1:n)=0
      do idchain=1,2
        root=ide(idchain,igp)
        cid(root)=idchain
        !push
        na=1
        list_act(na)=root
        visited(1:n)=.false.
        do p=1,nt
          v=list_act(na) !pop
          na=na-1
          if(.not.visited(v)) then
            visited(v)=.true.
            do i=1,nlink(v)
              w=link(v,i)
              if(visited(w)) cycle
              na=na+1
              list_act(na)=w
              cid(w)=idchain
            end do
          end if
          if(na.eq.0) exit
        end do
      end do

      return
      end subroutine split_chain


      subroutine check_tree(sn)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      dimension nlink_tmp(msize**2),
     *          link_tmp(msize**2,mlink)
      logical deleted(sn)

      nlink_tmp=nlink
      link_tmp=link
      deleted=.false.
              
30    do i=1,sn
        if(.not.deleted(i)) then
          if(nlink_tmp(i).le.1) then
            do j=1,nlink_tmp(i)
              k=link_tmp(i,j)
              if(.not.deleted(k)) then
                do l=1,nlink_tmp(k)
                  if(link_tmp(k,l).eq.i) then
                    link_tmp(k,l:nlink_tmp(k)-1)=
     *              link_tmp(k,l+1:nlink_tmp(k))
                    nlink_tmp(k)=nlink_tmp(k)-1
                    exit
                  end if
                end do
              end if
            end do
            deleted(i)=.true.
          end if
        end if
      end do
        
      do i=1,sn
        if(.not.(deleted(i)).and.nlink_tmp(i).le.1) then 
          goto 30
        end if
      end do

      do i=1,sn
        if(.not.(deleted(i))) then
          write (*,*) "# Found ring in tree!"
          stop
        end if
      end do

      return 
      end subroutine check_tree

      subroutine
     *quality_tree(n,edge_dens,ide,tree,ntree,score,
     *                              usedv,pathv,nv)
      use params,only:msize,mlink
      use commongraph
      implicit integer(i-z)
      logical visited(n),used(n),usedv(n)
      dimension ide(2,n*n),edge_dens(n*n),
     *          tree(n*n)
      dimension lact(n)
      dimension pathv(n),nextv(n)
      integer bestv
      real score,cost(n)

      ! DFS
      root=ide(1,tree(ntree))
      root0=-1
      do icount=1,100
        cost(1:n)=0
        !push
        na=1
        lact(na)=root
        visited(1:n)=.false.
        fmaxcost=0
        max_leaf=-1
        do p=1,ntree
          now=lact(na) !pop
          na=na-1
          if(.not.visited(now)) then
            visited(now)=.true.
            do i=1,nlink(now)
              next=link(now,i)
              if(visited(next)) cycle
              na=na+1
              lact(na)=next
              cost(next)=cost(now)+edge_dens(idelink(now,i))
              nextv(next)=now
              if(cost(next).ge.fmaxcost) then
                fmaxcost=cost(next)
                max_leaf=next
              end if
            end do
          end if

          if(na.eq.0) exit

        end do

        if(max_leaf.eq.root0) exit

        root0=root
        root=max_leaf
      end do

      score=fmaxcost**2 ! 1st longest path
      !trace back
      now=max_leaf

      used(1:n)=.false.
      nv=1
      pathv(nv)=now
      !1st path
      do i=1,ntree
        next=nextv(now)
        nv=nv+1
        pathv(nv)=next
        used(next)=.true.
        if(next.eq.root) exit
        now=next
      end do
      usedv(1:n)=used(1:n)
      usedv(root)=.false. !ignore terminal
        
      !find long branches 2nd 3rd,....100th,
      do icount=2,100
        best_cost=0
        bestv=-1
        do i=1,n
          if(used(i)) cycle
          if(nlink(i).ne.1) cycle
          !traca back
          now=i
          do j=1,n
            next=nextv(now)
            if(used(next)) exit
            now=next
          end do
          diff_cost=cost(i)-cost(next)
          if(best_cost.lt.diff_cost) then
            best_cost=diff_cost
            bestv=i
          end if
        end do
        if(bestv.eq.-1) exit
        score=score+best_cost**2
        
        !trace back
        now=bestv
        used(now)=.true.
        do j=1,n
          next=nextv(now)
          if(used(next)) exit
          used(next)=.true.
          now=next
        end do
      end do

      return
      end subroutine quality_tree

      recursive subroutine quick_sort(list, order,n)
      use params,only:msize
      implicit none
      real  list(msize**3)
      integer order(msize**3)

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

      function fkendall(p1,n1,p2,n2,n)
      use params,only:msize
      implicit integer(i-z)
      dimension p1(msize**2),p2(msize**2),indx(msize**2)
      dimension posi1(n,n) ,posi2(n,n)
      logical ex1(n),ex2(n)

      !init 
      ex1(1:n)=.false.
      ex2(1:n)=.false.
      indx(1:maxval(p2(1:n2)))=0
      do i=1,n2
        ex2(p2(i))=.true.
        indx(p2(i))=i!index of path2
      end do
        
      ntrue=0
      nfalse=0
      do i=1,n1,3
        o11=i!order
        i11=p1(i)!id11
        if(.not.ex2(i11)) continue
        o21=indx(i11)
        do j=i+1,n1,3
          o12=j
          i12=p1(j)
          if(.not.ex2(i12)) continue
          o22=indx(i12)
          !o11 < o12 so...
          if(o21.lt.o22) then
            ntrue=ntrue+1
          else
            nfalse=nfalse+1
          end if
        end do
      end do
        
      fkendall=float(2*ntrue)/float(ntrue+nfalse)-1.00
      
c      return 
      end function fkendall 

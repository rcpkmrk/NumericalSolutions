!------------------------------------------------------------------------------|
!..A 2-D FINITE VOLUME SOLVER FOR THE DIFFUSION EQUATION                       |
!------------------------------------------------------------------------------|
Module data
  integer        :: ncell, nnode, ntmax, ntio=5000
  real           :: dt
  real,parameter :: alpha=22.5e-6, delT_allow=1.E-5
  integer,allocatable,dimension(:,:) :: node, neigh
  real,allocatable,dimension(:)      :: Tcell,Tbc,area
  real,allocatable,dimension(:,:)    :: xy, Tgrad
End module

program FINITEVOLUME
   Use data
   call INIT            !..Read the input data and initialize the solution
   do nt=1,ntmax        !..Start the solution loop 
      call GRADIENT     !..Evaluate T gradients for all the cells
      delT_max = 0.
      do n = 1,ncell    !..Sweep all the cells and solve for T^n+1
         delT     = -dt/area(n) * CellFLUX(n)
         Tcell(n) = Tcell(n) + delT
         delT_max = MAX(ABS(delT),delT_max)
      enddo
      if( nt < 20 .or. MOD(nt,100) .eq. 0 ) print*,' nt, DelT_max :',nt,delT_max
      if(delT_max .lt. delT_allow) exit
      if( MOD(nt,ntio) .eq. 0 .and. nt .ne. ntmax ) call Qout(nt)
   enddo                        
   call Qout(nt-1)       !..Output the final solution
   stop  'DONE'
 end
!-----------------------------------------------------------------------
subroutine GRADIENT
  Use data
  Tgrad = 0.
  do n = 1,ncell
     do ns = 1,3   
        n2 = MOD(ns,3) + 1
        dx = xy(1,node(n2,n)) - xy(1,node(ns,n))
        dy = xy(2,node(n2,n)) - xy(2,node(ns,n))
        ne = neigh(ns,n)
        if(ne .gt. 0) then                  !..real neighbor
          Tface = 0.5*(Tcell(n) + Tcell(ne))
         elseif (ne .eq. -3) then            !..insulated wall
          Tface =  0.5*(Tcell(n) + Tcell(n))
        else
          Tface =  0.5*(Tcell(n) + Tbc(-ne))  !..other BC's
        endif
        Tgrad(1,n) = Tgrad(1,n) + Tface*dy / area(n)
        Tgrad(2,n) = Tgrad(2,n) - Tface*dx / area(n)
     enddo 
  enddo

return 
end
!-----------------------------------------------------------------------
function CELLFLUX(n)
  Use data

  Cellflux = 0.
  do ns = 1,3                    !..Loop over the edges
     n2 = MOD(ns,3) + 1
     dx = xy(1,node(n2,n)) - xy(1,node(ns,n))
     dy = xy(2,node(n2,n)) - xy(2,node(ns,n))
     ne = neigh(ns,n)
     if( ne .gt. 0 ) then        !..Real neighbor... 
         f = 0.5*( Tgrad(1,n) + Tgrad(1,ne) )
         g = 0.5*( Tgrad(2,n) + Tgrad(2,ne) )
     elseif( ne .lt. 0 ) then
         f =  0.5 * Tgrad(1,n)                              !..BCs
         g =  0.5 * Tgrad(2,n)
     endif
     edgeflux = - alpha*(f*dy - g*dx)
     Cellflux = Cellflux + edgeflux     !..Add the edge flux
  enddo
return
end
!-----------------------------------------------------------------------
 subroutine INIT
  Use data
  character :: fn*32
  logical   :: ok

!..Read the grid data
   write(*,'(/(a))',advance='no')'  Enter the grid filename [grid.dat]: '
   read(*,'(a)') fn
   if( fn .eq. ' ') fn = 'grid.dat'
   inquire(FILE=fn,EXIST=ok)
   if( .not. ok ) then
       print*, '  ', fn, ' does not exist! \n\n'
       stop
   endif
   open(1,file=fn,form='formatted')
   read(1,*) ncell,nnode
   allocate( node(3,ncell),neigh(3,ncell),xy(2,nnode),area(ncell), &
             Tcell(ncell),Tgrad(2,ncell),Tbc(4) )
   read(1,*) (no,(xy(i,n),i=1,2),n=1,nnode)
   read(1,*) (no,(node(i,n),i=1,3),(neigh(i,n),i=1,3),n=1,ncell)
   close(1)
   print*, ' # of cells :',ncell
   print*, ' # of nodes :',nnode

   write(*,'(/(a))',advance='no')'  Enter the time step and max step: '
   read(*,*) dt,ntmax

!..Compute cell areas
   do n = 1,ncell
      n1 = node(1,n)
      n2 = node(2,n)
      n3 = node(3,n)
      area(n) = 0.5*( (xy(1,n2)-xy(1,n1))*(xy(2,n3)-xy(2,n1))   &
                     -(xy(2,n2)-xy(2,n1))*(xy(1,n3)-xy(1,n1)) )
   enddo
!..Initialize the solution and BC, output the initial condition
   Tbc(1) = 0.
   Tbc(2) = 10.
   Tbc(3) = 0.   !..note that it is insulated, the value is not used!
   Tbc(4) = 50.
   Tcell(:) = 10.
   call QOUT(0)

return
end
!-------------------------------------------------------------------
subroutine QOUT(nstep)   !..Output the solution/grid in PLOT3D/TECPLOT format
  Use data
  real, dimension(nnode)    :: Tnode,Gx,Gy,Outflux
  character :: fname*32, str*12, ext*6
   call QNODES(Tnode,Gx,Gy,Outflux)
!..Set the output file name
   write(str,'(f9.6)') FLOAT(nstep)/1e6
   read(str,'(3x,a6)') ext
   fname = 'q-'//ext//'.plt'
   open(1,file=fname, form='formatted')
   write(1,100) nnode,ncell
   write(1,101) (xy(1,n),xy(2,n),Tnode(n),Gx(n),Gy(n),Outflux(n),n=1,nnode)
   write(1,102) (node(1,n),node(2,n),node(3,n),n=1,ncell)
   close(1)
 100 format (' VARIABLES= "X", "Y", "T", "Gx", "Gy", "Outflux" '/,     &
             ' ZONE N=', I6,' E=', I6,' F=FEPOINT ',' ET=triangle'  )
 101 format (6(1x,e12.5))
 102 format (3(1x,i6))
return
end
!-------------------------------------------------------------------
subroutine QNODES(Tnode,Gx,Gy,Outflux)      !..Evaluate averaged node values
   Use data
   integer, dimension(nnode) :: npass
   real, dimension(nnode)    :: Tnode
   real, dimension(nnode)    :: Gx
   real, dimension(nnode)    :: Gy
   real, dimension(nnode)    :: Outflux
   npass = 0
   Tnode = 0.
   Gx = 0.
   Gy = 0.
   Outflux = 0.
   do n=1,ncell                !..Find the contribution of cells to the nodes 
   do nf=1,3
      nn = node(nf,n)
      Tnode(nn)=Tnode(nn)+Tcell(n)
      Gx(nn)=Gx(nn)+Tgrad(1,n)
      Gy(nn)=Gy(nn)+Tgrad(2,n)
      Outflux(nn) = sqrt(Gx(nn)**2+Gy(nn)**2)
      npass(nn)=npass(nn)+1
   enddo

   enddo
   Tnode(:)=Tnode(:)/npass(:)  !..Node averaged value
   Gx(:)=Gx(:)/npass(:)
   Gy(:)=Gy(:)/npass(:)
return 
end
!------------------------------END----------------------------------

!------------------------------------------------------------------------------|
!..A 2-D FINITE VOLUME SOLVER FOR THE DIFFUSION EQUATION                       |
!------------------------------------------------------------------------------|
Module data
  integer  :: ncell, nnode, ntmax, ntout=100
  real,parameter :: alpha_deg=5.0, delq_allow=1.E-5, pi=4*ATAN(1.)
  integer,allocatable,dimension(:,:) :: node, neigh
  real,allocatable,dimension(:,:) :: xy, qgrad
  real,allocatable,dimension(:) :: qcell,area,qbc, qnode,unode,vnode,cpnode
  real :: dt,uinf,vinf
End module

program FINITEVOLUME
   Use data
!..Read the input data and initialize the solution
   call INIT
!..Start the solution loop 
   do nt=1,ntmax
      call GRADIENT      !..Evaluate q gradients for all cells
      delq_max = 0.
      do n = 1,ncell     !..Sweep all the cells and solve for q^n+1
         delq     = -dt/area(n) * CellFLUX(n)
         qcell(n) = qcell(n) + delq
         delq_max = MAX(ABS(delq),delq_max)
      enddo
!..Output the intermediate solutions 
      if( nt < 10 .or. MOD(nt,100) .eq. 0 ) print*, ' nt, Delq_max :',nt,delq_max
      if( delq_max .lt. delq_allow) exit
      if( MOD(nt,ntout) .eq. 0 .or. nt .eq. ntmax ) call QOUT(nt)
   enddo
   call QOUT(nt-1)
   stop  'DONE'
 end
!------------------------------------------------------------------------
 subroutine INIT()
  Use data
  character :: fn*16
  logical   :: ok
!..Read the grid data
   write(*,'(/(a))',advance='no')'  Enter the grid filename [grid.dat]: '
   read(*,'(a)') fn
   if( fn .eq. ' ') fn = 'grid.dat'
   inquire(FILE=fn,EXIST=ok)
   if( .not. ok ) then
       print*, '  ', fn, ' does not exist! '
       stop
   endif 
   open(1,file=fn,form='formatted')
   read(1,*) ncell,nnode
   allocate( node(3,ncell),neigh(3,ncell),xy(2,nnode),area(ncell), &
             qcell(ncell),qnode(nnode),unode(nnode),vnode(nnode),cpnode(nnode), &
             qbc(4),qgrad(2,ncell) )
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

   alpha = alpha_deg/180.*pi
   uinf  = COS(alpha)
   vinf  = SIN(alpha)
!..Initialize the solution with 0
!   qcell = 0.
!..Initialize the solution with free stream
  do n =1,ncell
     xc = (xy(1,node(1,n))+xy(1,node(2,n))+xy(1,node(3,n)))/3.
     yc = (xy(2,node(1,n))+xy(2,node(2,n))+xy(2,node(3,n)))/3.
    qcell(n) = xc*uinf + yc*vinf
  enddo

   call GRADIENT
   call Qout(0)

   return 
 end
!-------------------------------------------------------------------
subroutine  GRADIENT()
  Use data
  qgrad = 0.
  do n = 1,ncell
     do ns = 1,3   
        n2 = MOD(ns,3) + 1
        dx = xy(1,node(n2,n)) - xy(1,node(ns,n))
        dy = xy(2,node(n2,n)) - xy(2,node(ns,n))
        ne = neigh(ns,n)
        if(ne .gt. 0) then                     !..real neighbor
           qface = 0.5*(qcell(n) + qcell(ne))
        elseif (ne .eq. -1) then           !..Farfield BC
           qgrad(1,n) = -uinf
           qgrad(2,n) = -vinf
           exit
        else                                   !..wall BC
           qface = qcell(n)
        endif
        qgrad(1,n) = qgrad(1,n) + qface*dy /area(n)
        qgrad(2,n) = qgrad(2,n) - qface*dx /area(n)
     enddo 
  enddo
  return 
 end
!------------------------------------------------------------------------
function CellFLUX(n)
  Use data

  Cellflux = 0.
  do ns = 1,3               !..Add the surface fluxes
     n2 = MOD(ns,3) + 1
     dx = xy(1,node(n2,n)) - xy(1,node(ns,n))
     dy = xy(2,node(n2,n)) - xy(2,node(ns,n))
!..Apply the proper BCs when computing the face fluxes
     ne = neigh(ns,n)
     if( ne .gt. 0 ) then            !..real neighbor... 
         f = -0.5*(qgrad(1,n) + qgrad(1,ne))
         g = -0.5*(qgrad(2,n) + qgrad(2,ne))
     elseif (ne .eq. -1) then        !..Farfield BC
         f = -uinf
         g = -vinf
     else                            !..wall BC
         f = 0.
         g = 0.
     endif
     Cellflux = Cellflux + (f*dy - g*dx)
  enddo 
return
end
!-------------------------------------------------------------------
subroutine  QOUT(nstep)     !..Output the solution/grid in TECPLOT format
  Use data
  character :: fname*32, str*12, ext*6
   call QNODES()
!..Set the output file name
   write(str,'(f9.6)') float(nstep)/1e6
   read(str,'(3x,a6)') ext
   fname = 'q-'//ext//'.plt'
   open(1,file=fname, form='formatted')
   write(1,100) nnode,ncell
   write(1,101) (xy(1,n),xy(2,n),qnode(n),unode(n),vnode(n),cpnode(n),n=1,nnode)
   write(1,102) (node(1,n),node(2,n),node(3,n),n=1,ncell)
   close(1)
  100 format (' VARIABLES= "X", "Y", "Phi", "u", "v", "Cp"'/,     &
              ' ZONE N=', I6,' E=', I6,' F=FEPOINT ',' ET=triangle'  )
  101 format (6(1x,e12.5))
  102 format (3(1x,i6))
return 
end
!-------------------------------------------------------------------
subroutine  QNODES()     !..Evaluate averaged q_node values
  Use data
  integer, dimension(:) :: npass(nnode) 
   qnode = 0.
   unode = 0.
   vnode = 0.
   npass = 0
!..Find the contribution of cells to the node Q
   do n=1,ncell
   do nf=1,3
      nn = node(nf,n)
      qnode(nn)=qnode(nn)+qcell(n)
      unode(nn)=unode(nn)+qgrad(1,n)
      vnode(nn)=vnode(nn)+qgrad(2,n)
      npass(nn)=npass(nn)+1
   enddo
   enddo
!..Average the total node Q with # of contributing cells 
   qnode(:)=qnode(:)/npass(:)
   unode(:)=unode(:)/npass(:)
   vnode(:)=vnode(:)/npass(:)
   cpnode(:)=1.-(unode(:)**2+vnode(:)**2)
return 
end
!------------------------------END----------------------------------

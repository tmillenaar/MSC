Program advDiv

Implicit none


!Declare variables
!------------------------------------------------------
integer :: i				!Used as index for general loops
integer :: j				!Used as index for the basin number, 1 = main basin, 2 = secund. basin
integer :: n				!Used to keep track of output time
! integer :: info
integer :: nmax				!Max number of months
integer :: ibcl				!Used to switch between forms of boundary condition; top
integer :: ibcr				!Used to switch between forms of boundary condition; bottom
integer :: t				!time in s

integer, dimension(2) :: imax		!
! integer, dimension(2) :: imld		!Depth of unstable column (from top down)
integer, dimension(2) :: imldbottom	!"Depth" of unstable column from (from bottom up)
integer, dimension(2) :: iconv		!

real :: yr2sec
real :: pi
real :: dx
real :: x
real :: dt
real :: tprev
real :: tmax
real :: tout
real :: dtout
real :: f
real :: minty
real :: bcr
real :: bcr1
real :: Dmax
real :: Wlimit
real :: flowdepth
real :: progress_output
real :: d_progress_output
real :: ramp
real :: heatFlux
real :: Cp
real :: Tbot
real :: Sbot
real :: Vol1_width
real :: Vol2_width
real :: Vol1_length
real :: Vol2_length
real :: Vol1
real :: Vol2
real :: density
real :: dtmax
real :: tau
real :: p
real :: pb
real :: pc
real :: Db
real :: Dc
real :: strait_depth
real :: strait_width
real :: U_lat
real :: Salt_diff_amplifier
real :: evapamp1
real :: evapamp2
real :: Temp_diff_amplifier
real :: Smean
real :: Tmean
real :: mean_evap_rate

real, dimension (2) :: L
real, dimension (2) :: tconv
real, dimension (2) :: evaporation
real, dimension (2) :: Rho_mean_flowdepth
real, dimension (2) :: Vol

real*16 :: bcl
real*16 :: bcl1

real*16, dimension (:,:), allocatable :: Temp
real*16, dimension (:,:), allocatable :: Sal
real*16, dimension (:,:), allocatable :: TempPrev
real*16, dimension (:,:), allocatable :: SalPrev
real*16, dimension (:,:), allocatable :: TempOut
real*16, dimension (:,:), allocatable :: SalOut
real*16, dimension (:,:), allocatable :: TempSalPrevOut
real*16, dimension (:,:), allocatable :: w_Temp
real*16, dimension (:,:), allocatable :: Temp1
real*16, dimension (:,:), allocatable :: Sal1
real*16, dimension (:,:), allocatable :: w
real*16, dimension (:,:), allocatable :: SS_S
real*16, dimension (:,:), allocatable :: SS_T
real*16, dimension (:,:), allocatable :: rho
real*16, dimension (:,:), allocatable :: D
real*16, dimension (:,:), allocatable :: D0
real*16, dimension (:,:), allocatable :: Output
real*16, dimension (:,:), allocatable :: PrevOutput
real*16, dimension (:,:), allocatable :: D_T
real*16, dimension (:,:), allocatable :: w_salt
real*16, dimension (:,:), allocatable :: imld
character :: ext, ext2*2, ext3*3, ext4*4, ext5*5, ext6*6!, outputnr*6
character*80 :: fna
! character*80 :: fnb, fnc, fnd, fne, fnf
character*80 :: comment

logical :: spinUp
logical :: noConv
logical :: check_conv
logical :: no_conv_in_col

common /bcs/ ibcl, ibcr, bcr, bcr1

!Setting the Variables
!------------------------------------------------------
!yr2sec= 365.25*24.e0*3600.e0
yr2sec= 360.00*24.e0*3600.e0
pi= 4.*atan(1.e0)
! Cp= 4.e3			! specific heat, used in heat flux

!Set update interval of progress bar
d_progress_output = (1./36.)*yr2sec
progress_output = d_progress_output


Dmax = max(maxval(D),maxval(D_T)) 
dt = 0.9*(dx*dx)/(2.e0*Dmax)
dtout= (1./12.)*yr2sec		! interval at which to write output

open (unit=15, file="setup_data.txt", action='read' )
read (15, *)  comment, L(1), comment, L(2), comment, tmax, comment, Cp, comment, dx, comment, &
      pb, comment, pc, comment, tau, comment, Tbot, comment, Sbot, comment, flowdepth, comment, &
      strait_depth, comment, strait_width, comment, U_lat,comment,Vol1_width,comment,Vol2_width, comment,&
      Vol1_length,comment,Vol2_length,comment,Salt_diff_amplifier, comment, Temp_diff_amplifier, comment, &
	  evapamp1, comment, evapamp2, comment, mean_evap_rate
tmax = tmax*yr2sec
! ! w0= -0./yr2sec		! advection velocity [m/s] (upward velocity = -ve)

! to keep an eye on equation of state
! write (*,'(a,f7.2)') 'density [kg/m3] at T= 20 C and S= 277: ', density(20.,277.)
!
! it is efficient to define these explicitely
!
Db= 10**pb
Dc= 10**pc

Vol(1) = Vol1_width*Vol1_length*dx		!Vol=width*length*dx NOTE: put dx here and not flowdepth ...
Vol(2) = Vol2_width*Vol2_length*dx			!...since larger flowdepth does not mean more volume	

! ! write (*,'(a,f8.2)') 'FTCS: max time step set by diffusion [s] = ', (dx*dx)/(2.*Dc)
! ! if (w0 .ne. 0.) then
! !    write (*,'(a,f8.2)') 'FTCS: max time step set by advection [s] = ', (2.*Db)/(w0*w0)
! !    write (*,'(a,f8.2)') 'CN: max time step set by advection [s] =   ', dx/w0
! !    write (*,'(a,f8.2)') 'background Peclet number main basin [-] =             ', (abs(w0)*L(1))/Db
! !    write (*,'(a,f8.2)') 'background Peclet number secundary basin [-] =             ', (abs(w0)*L(2))/Db
! ! endif

!Allocate allocatables
!------------------------------------------------------
! find required number of nodes and allocate array space
imax(1)= nint( (L(1)-dx)/dx )	!Number of cells in main basin column
imax(2)= nint( (L(2)-dx)/dx )	!Number of cells in secundary basin column
nmax= nint( tmax/dtout )
write(*,*) 'imax(1)= ', imax(1), '    imax(2)= ', imax(2)

allocate ( Temp(0:imax(1)+1,2), Sal(0:imax(1)+1,2), imld(0:imax(1)+1,2) )
allocate ( TempPrev(0:imax(1)+1,2), SalPrev(0:imax(1)+1,2) )
allocate ( Temp1(0:imax(1)+1,2), Sal1(0:imax(1)+1,2) )
allocate ( TempOut(0:imax(1)+1,2), SalOut(0:imax(1)+1,2) )
allocate ( Output(0:imax(1)+1,2), PrevOutput(0:imax(1)+1,2), w_Temp(0:imax(1)+1,2) )
allocate ( w(0:imax(1)+1,2), SS_S(0:imax(1)+1,2), SS_T(0:imax(1)+1,2), w_salt(0:imax(1)+1,2) )
allocate ( rho(0:imax(1)+1,2), D(0:imax(1)+1,2), D0(0:imax(1),2), D_T(0:imax(1)+1,2) )


! initialise 
!------------------------------------------------------
do i=0,imax(1)+1
  do j=1,2
    Temp(i,j)= Tbot
    Sal(i,j)= Sbot
    D(i,j)= Db*Salt_diff_amplifier
    D0(i,j)= Db
    D_T(i,j)= Db*Temp_diff_amplifier  !! 2 orders of magnitude quicker then salt diffusion [schmitt1994double] 
  end do
end do


! initialise other arrays
!------------------------------------------------------
do i= 0,imax(1)+1
  do j= 1,2
    if (j == 1) then
      Temp1(i,j)= 0.e0
      Sal1(i,j)= 0.e0
      SalPrev(i,j)= 0.e0
      TempPrev(i,j)= 0.e0
      PrevOutput(i,j)= 0.e0
      w(i,j)=0.e0
      w_Temp(i,j)=0.e0
      w_salt(i,j) = 0.e0
      SS_S(i,j)= 0.e0
      SS_T(i,j)= 0.e0
      rho(i,j)= density( Temp1(i,j), Sal1(i,j) )
    else if (i < imax(2)+1) then
      Temp1(i,j)= 0.e0
      Sal1(i,j)= 0.e0
      SalPrev(i,j)= 0.e0
      TempPrev(i,j)= 0.e0
      PrevOutput(i,j)= 0.e0
      SS_S(i,j)= 0.e0
      SS_T(i,j)= 0.e0
      w(i,j)=0.e0
      w_Temp(i,j)=0.e0
      w_salt(i,j) = 0.e0
      rho(i,j)= density( Temp1(i,j), Sal1(i,j) )
    endif
  end do
end do


! write initial profile
!
open (20,file= 'Main_basin_output_0.xy')
open (30,file= 'Secundary_basin_output_0.xy')
do i= 0,imax(1)+1
   write (20,*) real(i)*dx, Temp(i,1), Sal(i,1), density( Temp(i,1), Sal(i,1) )
   if (i < imax(2)+1) then
      write (30,*) real(i)*dx, Temp(i,2), Sal(i,2), density( Temp(i,2), Sal(i,2) )
   endif
end do
close (20)
close (30)


open (11, file= 'main_mean_values.xy')
open (12, file= 'secundary_mean_values.xy')
open (21,file= 'main_mixing.xy')
open (22,file= 'secundary_mixing.xy')
open (31,file= 'main_forcing.xy')
open (32,file= 'secundary_forcing.xy')

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop
!
t= 0.e0
tout= dtout
tprev= 0.e0
noConv= .true.
no_conv_in_col = .true.
evaporation= 0.e0
iconv= 0
tconv= 0.e0
spinUp= .true.
minty = 0
check_conv=.false.

do while (t.lt.tmax) 		!Main loop
    ramp= min( t/(4*yr2sec), 1.e0 )
    n= nint(tout/dtout)
    if (t > 10.e0*yr2sec) then
	if (spinUp) then
	    spinUp= .false.
	    write(*,*) ''
	    write(*,*) 'End of ramp.   Start of convection, advection and source-sink'
	    write(*,*) ''
	    write(*,*) 'Progress convection calculations (runspeed is dependent on diff. vel. i.e. quick = no convection):'
	    call progress(t,tmax,10.*yr2sec,1,n,nint(tmax/dtout)) ! reset the progress bar.
	endif
	if (t .ge. progress_output) then
	    call progress(t,tmax,10.*yr2sec,0,n,nint(tmax/dtout)) ! generate the progress bar.
	    progress_output = progress_output + d_progress_output
	endif
	noConv= .false.
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!  Checking for, calculating and setting: Source-Sink and Advection  !!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!!!!!Determine average density in flowdepth
	do i=0,nint(flowdepth)
	    Rho_mean_flowdepth(1) = Rho_mean_flowdepth(1) + density( Tbot*Temp1(imax(1)-i,1), Sbot*Sal1(imax(1)-i,1) )
	    Rho_mean_flowdepth(2) = Rho_mean_flowdepth(2) + density( Tbot*Temp1(imax(2)-i,2), Sbot*Sal1(imax(2)-i,2) )
	end do
	Rho_mean_flowdepth(1) = Rho_mean_flowdepth(1)/(flowdepth+1)
	Rho_mean_flowdepth(2) = Rho_mean_flowdepth(2)/(flowdepth+1)
	if (Rho_mean_flowdepth(2) > Rho_mean_flowdepth(1) .and. t>11*yr2sec ) then !Initate flow since "density in bottom of basin 2" > "density in bottom basin 1"
	    Rho_mean_flowdepth(1)=0. 
	    Rho_mean_flowdepth(2)=0.
	    ! do i=0,nint(flowdepth)-1
! ! 		!!Assuming Vin = Vout en Utop=Ubot 
! ! 		!!U_lat only determines the transport amount of Vol of water/sec, there is no travel distance or &
! ! 		!! 	travel time of the lateral water flow in this model (for now at least)
		! SS_S(imax(1)-i,1)= max(0.0, (1./(real(flowdepth)))*(Sal(imax(2)-i,2) * strait_depth * strait_width * U_lat)/Vol(1) )	! @ maxdepth of basin 1: SS= + salinity at maxdepth basin 2
		! SS_S(imax(2)-i,2) = min( 0.0, (1./(real(flowdepth)))*(-1.)*(Sal(imax(2)-i,2) * strait_depth * strait_width * U_lat)/Vol(2) )	! @ maxdepth of basin 2: SS= - salinity at maxdepth basin 2
		! SS_S(i,1)= min( 0.0, (1./(real(flowdepth)))*(-1.)*(Sal(i,1) * strait_depth * strait_width * U_lat)/Vol(1))!8.**(-6) )		! @ surface of basin 1: SS= - salinity at surface basin 1
		! SS_S(i,2)= max( 0.0 , (1./(real(flowdepth)))*(Sal(i,1) * strait_depth * strait_width * U_lat)/Vol(2) )			! @ surface of basin 2: SS= + salinity at surface basin 1
	
		! SS_T(imax(1)-i,1)= max(0.0, (1./(flowdepth))*(Temp(imax(2)-i,2) * strait_depth * strait_width * U_lat)/Vol(1) )	! @ maxdepth of basin 1: SS= + salinity at maxdepth basin 2
		! SS_T(imax(2)-i,2) = min( 0.0, (1./(flowdepth))*(-1.)*(Temp(imax(2)-i,2) * strait_depth * strait_width * U_lat)/Vol(2) )	! @ maxdepth of basin 2: SS= - salinity at maxdepth basin 2
		! SS_T(i,1)= min( 0.0, (1./(flowdepth))*(-1.)*(Temp(i,1) * strait_depth * strait_width * U_lat)/Vol(1))!8.**(-6) )		! @ surface of basin 1: SS= - salinity at surface basin 1
		! SS_T(i,2)= max( 0.0 , (1./(flowdepth))*(Temp(i,1) * strait_depth * strait_width * U_lat)/Vol(2) )			! @ surface of basin 2: SS= + salinity at surface basin 1
	    ! end do
	    do i=1,imax(1)
! ! ! 	!set upward advection of main column to that of the lateral velocity (water up/m^2  = water in/surface area vol1)
		w(i,1)=-dx*(U_lat*strait_depth*strait_width)/(Vol(1))
		w_salt(i,1)=w(i,1)*Sal(i,1)
		w_Temp(i,1)=w(i,1)*Temp(i,1)
		if (i<imax(2)+1) then
		    w(i,2)=dx*(U_lat*strait_depth*strait_width)/(Vol(2))		!set downward advection of secundary column ''
		    w_salt(i,2)=w(i,2)*Sal(i,2)
		    w_Temp(i,2)=w(i,2)*Temp(i,2)
		end if
	    end do
	else !In which case the density trigger is not activated
	    do i=0,imax(1)	! Set advection and SS to 0 since there is no flow
		w(i,1)=0.	
		w_salt(i,1)=0.
		w_Temp(i,1)=0.
		SS_S(i,1)= 0.
		SS_T(i,1)= 0.
		if (i<imax(2)+1) then
		    w(i,2)=0.	
		    w_salt(i,2)=0.
		    w_Temp(i,2)=0.
		    SS_S(i,2) = 0.
		    SS_T(i,2) = 0.
		end if
	    end do
	endif
	
	else !In which case: t < 10 yr
	    if (t .ge. progress_output) then
	    call progress(t,10.*yr2sec,minty,0,n,nint(tmax/dtout)) ! generate the progress bar.
	    progress_output = progress_output + d_progress_output
	endif
	noConv= .true.
    endif


!  !!!Obtain size of next timestep
    if ( max(maxval(w_salt),maxval(w_Temp)) /= 0 ) then  
      Dmax = max(maxval(D),maxval(D_T))
      Wlimit = (2.*Dmax)/(max( (maxval(w_salt)*maxval(w_salt)), (maxval(w_Temp)*maxval(w_Temp)) ))
      dt = min(0.9*(dx*dx)/(2.e0*Dmax), 0.9*Wlimit )	!timestep set by diffusion and advection
      if (.NOT. check_conv) check_conv=.true.
    else
      check_conv=.false.
      Dmax = max(maxval(D),maxval(D_T)) 
      dt = 0.9*(dx*dx)/(2.e0*Dmax)		!timestep set by diffusion only
    endif
    if ( dt < 0.0005 ) then !To avoid extreme runtimes. It generally works fine in practice but may bring artifacts. &
			      ! If this happens, the results are clearly wrong or 'NAN'
	dt = 0.0005
	print*, 'Error, dt too small, dt set to 0.0005s for this iteration'
    endif
     
		do j= 1,2			!Dual basin loop, where j(1)= main basin and j(2)= secundary (shallower) basin

		! (evaporation in m/year)
		  if (j==1) then
		    evaporation(j)= ramp*0.3*sin((t/yr2sec)*2.*pi)+mean_evap_rate
		  else if (j==2) then
		    evaporation(j)= ramp*0.3*sin((t/yr2sec)*2.*pi)+mean_evap_rate
		  end if

		! set boundary conditions
		! ibcl = type of bc at "l"eft = at i= 0
		! ibcl= 0 == specified T or S (value = bcl and bcl1 at next time step)
		! ibcl= 1 == specified flux (value = idem)
		! likewise for ibcr, bcr, and bcr1
		!
		
		!!! TEMP
		
		bcl= ramp*(40.e0/( density(Temp(0,j),Sal(0,j))*Cp ))*sin( (t/yr2sec)*2.e0*pi )
		heatFlux= bcl*density(Temp(0,j),Sal(0,j))*Cp
		ibcl= 1
		ibcr= 1
		bcl1= bcl
		bcr1= bcr
		bcr= 0.e0



		call FTCS (imax(j), dx, dt, Temp(:,j), Temp1(:,j), SS_T(:,j), D_T(:,j), w_Temp(:,j),bcl1,Db)

		!!!SALINITY

		bcl= (evaporation(j)/yr2sec)*Sal(0,j)
		bcl1= bcl

		call FTCS (imax(j), dx, dt, Sal(:,j), Sal1(:,j), SS_S(:,j), D(:,j), w_salt(:,j), bcl1,Db)

		! calculate density profile
		!
		do i= 0,imax(j)+1
		     rho(i,j)= density( Temp1(i,j), Sal1(i,j) )
		end do

	! account for convection
		!
		
	if (spinUp) then
	  iconv(j)= 0
	  do i= 0,imax(j)+1
	    imld(i,j)= 0
	    p= max( pb, pc + (pb-pc)*(tconv(j)/tau) )
	    if (noConv) p= pb
	    D(i,j)= Salt_diff_amplifier*10**p
	    D_T(i,j)= Temp_diff_amplifier*10**p
	  end do
	else
	i=0
	  do while(i<imax(j))  
! 	    threshold = rho(i+1,j)
! 	    threshold = threshold + 10.
	    if (rho(i,j) > (rho(i+1,j) +0.0005)) then
	      if (iconv(j).eq.0) then
			tconv(j)= dt
	      else
			tconv(j)= tconv(j) + dt
	      endif
	      p= min( pc, pb + (pc-pb)*(tconv(j)/tau) )
	      if (noConv) p= pb
	      D(i,j)= 10.**p
	      ! D(i-1,j)= 10.**p
	      ! D(i+1,j)= 10.**p
	      D_T(i,j)= 10.**p
	      ! D_T(i-1,j)= 10.**p
	      ! D_T(i+1,j)= 10.**p
	      imld(i,j)=i
	      ! imld(i+1,j)=i+1
	      ! imld(i-1,j)=i-1
	      no_conv_in_col = .false.
	      i=i+1
	    else
	      imld(i,j)=0
	      
	      D(i,j)= Db*Salt_diff_amplifier
	      D_T(i,j)= Db
	      
	      i=i+1
	      
! 	      p= max( pb, pc + (pb-pc)*(tconv(j)/tau) )
! 	      if (noConv) p= pb
! 	      D(i,j)= Salt_diff_amplifier*10**p
! 	      D_T(i,j)= 10**p 
	    end if
	  end do
	  if (no_conv_in_col) then
	    iconv(j)=0
	    	D(i,j)= Db*Salt_diff_amplifier
		D_T(i,j)= Db
	  else
	    iconv(j) = iconv(j) + 1
	    do i= 0,imax(j)
		p= min( pc, pb + (pc-pb)*(tconv(j)/tau) )
		if (noConv) p= pb
		    D(i,j)= 10.**p
		    D_T(i,j)= 10.**p
	    end do
! 	    if (rho(imax(j),j) > (rho(imax(j)-1,j) -0.001)) then
! 	      i=0
! 	      do while ((rho(imax(j),j) > (rho(imax(j)-i-1,j)-0.001)) .and. i <= imax(j)-10)
! 		i=i+1
! 	      end do
! 	      imldbottom(j) = i
! 	      i=0
! ! 	      print*, imax(j)-imldbottom(j)
! 	      if (i>2) then
! 		do i= (imax(j)-imldbottom(j)),imax(j)
! 		    D(i,j)= Db*Salt_diff_amplifier
! 		    D_T(i,j)= Db
! 		end do
! 	      endif
! 	    endif
	   end if
	  no_conv_in_col = .true.
	end if
	
	
! 		if (rho(0,j) > rho(1,j) .and. .NOT. noConv) then
! 
! 		      if (iconv(j).eq.0) then
! 			 tconv(j)= dt
! 		      else
! 			 tconv(j)= tconv(j) + dt
! 		      endif
! 		      iconv(j)= 1
! 		! assess depth extent of unstable column
! 		!
! 		      i= 1
! 		      do while (rho(0,j) > rho(i+1,j) .and. i <= imax(j))
! 			 i= i+1
! 		      end do
! 		      imld(j)= i
! 		      
! 		      do i= 0,imld(j)
! 			 p= min( pc, pb + (pc-pb)*(tconv(j)/tau) )
! 			 if (noConv) p= pb
! 			 D(i,j)= 10.**p
! 			 D_T(i,j)= 10.**p
! 		      end do
! 		      do i= imld(j)+1,imax(j)+1
! 			 D(i,j)= Db*Salt_diff_amplifier
! 			 D_T(i,j)= Db*Temp_diff_amplifier
! 		      end do
! 
! 		else !if (.NOT. spinUp) then
! 
! 		      if (iconv(j).eq.1) then
! 			 tconv(j)= dt
! 		      else
! 			 tconv(j)= tconv(j) + dt
! 		      endif
! 		      iconv(j)= 0
! 
! 		      imld(j)= 0
! 		      do i= 0,imax(j)+1
! 			 p= max( pb, pc + (pb-pc)*(tconv(j)/tau) )
! 			 if (noConv) p= pb
! 			 D(i,j)= Db*Salt_diff_amplifier !10**p*Salt_diff_amplifier
! 			 D_T(i,j)= Db*Temp_diff_amplifier !10**p*Temp_diff_amplifier
! 		      end do
! 
! 		endif
! 		imldbottom(j) = 0
! 		if (rho(imax(j),j) < (rho(imax(j)-1,j) - 0.0001) .and. .NOT. noConv) then
! 		      i=0
! 		      do while ((rho(imax(j),j) < (rho(imax(j)-i-1,j)-0.0001)) .and. i <= imax(j))
! 			i=i+1
! 			D(imax(j)-i,j)= 10**pc
! 			D_T(imax(j)-i,j)= 10**pc
! 		      end do
! 		      imldbottom(j) = i
! 		      i=0
! 	! 	      print*, imax(j)-imldbottom(j)
! 	! 	      if (imldbottom(j)<imax(j)-1) then
! 	! 		do i= (imax(j)-imldbottom(j)),imax(j)
! 	! 		    D(i,j)= 10**pc
! 	! 		    D_T(i,j)= 10**pc
! 	! 		end do
! 	! 	      endif
! 		endif
		   
!		   if (t.gt.(11.25*yr2sec)) write (30+j,*) ((t/yr2sec)-11)*1e3, imld, -log10(D(0,j)), D(0,j)
		   do i= 0,imax(j)+1
			SalPrev(i,j)= Sal(i,j)
			TempPrev(i,j)= Temp(i,j)
			Sal(i,j)= Sal1(i,j)
			Temp(i,j)= Temp1(i,j)
		   end do
	

		end do ! End dual basin loop
	   	tprev = t
		t= t + dt
		! write
		!
		if (t.ge.tout) then
			do j=1,2 ! Second go of dual basin loop

			      !write (20+j,*) t/dtout, imld(j)!, -log10(D(0))
! 			      write (20+j,*) t/dtout, imld(j), imax(j)+1-imldbottom(j)
			! find h at output time by interpolation
			!
			      f= (tout - tprev)/(t - tprev)
			      do i= 0,imax(j)+1
				    if (imld(i,j) /= 0) write (20+j,*) t/dtout, i
				    SalOut(i,j)= SalPrev(i,j) + f*( Sal(i,j) - SalPrev(i,j) )
				    TempOut(i,j)= TempPrev(i,j) + f*( Temp(i,j) - TempPrev(i,j) )
			      end do

			      n= nint(tout/dtout)

			      if (j==1) then
! 				      if (mod(n,10).eq.0) write (*,'(a,i3)') 'writing output time       n= ', n

				      if (n.lt.10) then
					 write (ext,'(i1)') n
					 fna= 'Main_basin_output_'//ext//'.xy'
				      else if (n.ge.10.and.n.lt.100) then
					 write (ext2,'(i2)') n
					 fna= 'Main_basin_output_'//ext2//'.xy'
				      else if (n.ge.100.and.n.lt.1000) then
					 write (ext3,'(i3)') n
					 fna= 'Main_basin_output_'//ext3//'.xy'
				      else if (n.ge.1000.and.n.lt.10000) then
					 write (ext4,'(i4)') n
					 fna= 'Main_basin_output_'//ext4//'.xy'
				      else if (n.ge.10000.and.n.lt.100000) then
					 write (ext5,'(i5)') n
					 fna= 'Main_basin_output_'//ext5//'.xy'
				      else if (n.ge.100000.and.n.lt.1000000) then
					 write (ext6,'(i6)') n
					 fna= 'Main_basin_output_'//ext6//'.xy'
				      else
					 write (*,*) 'error: file name cannot be constructed.'
					 stop
				      endif

				else if (j==2) then
!				      if (mod(n,10).eq.0) write (*,'(a,i3)') 'writing secundary basin output time  n= ', n

				      if (n.lt.10) then
					 write (ext,'(i1)') n
					 fna= 'Secundary_basin_output_'//ext//'.xy'
				      else if (n.ge.10.and.n.lt.100) then
					 write (ext2,'(i2)') n
					 fna= 'Secundary_basin_output_'//ext2//'.xy'
				      else if (n.ge.100.and.n.lt.1000) then
					 write (ext3,'(i3)') n
					 fna= 'Secundary_basin_output_'//ext3//'.xy'
				      else if (n.ge.1000.and.n.lt.10000) then
					 write (ext4,'(i4)') n
					 fna= 'Secundary_basin_output_'//ext4//'.xy'
				      else if (n.ge.10000.and.n.lt.100000) then
					 write (ext5,'(i5)') n
					 fna= 'Secundary_basin_output_'//ext5//'.xy'
				      else if (n.ge.100000.and.n.lt.1000000) then
					 write (ext6,'(i6)') n
					 fna= 'Secundary_basin_output_'//ext6//'.xy'
				      else
					 write (*,*) 'error: file name cannot be constructed.'
					 stop
				      endif
				else
					write(*,*) 'error, j out of bounds, write file cannot be constructed'
				end if
				!Writing log file for in and output via S term:
		!		write(44,*) t/yr2sec, '    in S bot:', s(imax),'    out S top:', S(0),'    in D top:', S_daughter(0), &
		!			'    out D bot:', S_daughter(imax_daughter)
		!		write(45,*) 's_in-s_out:',s(imax)-s(0),s(imax)-s(0)+evaporation,&
		!	'	D_in-D_out:',s_daughter(0)-s_daughter(imax_daughter),s_daughter(0)-s_daughter(imax_daughter)+evaporation


			      open (20,file= fna)
			      do i= 0,imax(j)+1
				 write (20,*) real(i)*dx, TempOut(i,j), SalOut(i,j), density( TempOut(i,j), SalOut(i,j))
			      end do
			      close (20)


			     Tmean= Temp(0,j)*dx/2.e0
			     Smean= Sal(0,j)*dx/2.e0
			!      HC= h(0,1)*Cp*density(h(0,1),h(0,2))*dx/2.e0
			     do i= 1,imax(j)
			        Tmean= Tmean + TempOut(i,j)*dx
			        Smean= Smean + SalOut(i,j)*dx
			        ! HC= HC + h(i,1)*Cp*density(h(i,1),h(i,2))*dx
			     end do
			     Tmean= Tmean + Temp(imax(j)+1,j)*dx/2.e0
			     Smean= Smean + Sal(imax(j)+1,j)*dx/2.e0
			     Tmean= Tmean/L(j)
			     Smean= Smean/L(j)
				 write(10+j,*) t, Tmean, Smean
			!      HC= HC + h(imax+1,1)*Cp*density(h(imax+1,1),h(imax+1,2))*dx

			!      Tmean2= 0.e0
			!      Smean2= 0.e0
			!      do i= 101,220
			!         Tmean2= Tmean2 + hout(i,1)
			!         Smean2= Smean2 + hout(i,2)
			!      end do
			!      Tmean2= Tmean2/120.
			!      Smean2= Smean2/120.

			   !   write (22,*) t/yr2sec, Tmean, Smean, h(200,2)
			      write (30+j,*) t/yr2sec, heatFlux, evaporation(j)

			! now set next time at which to write output
			!
			end do ! End second go of dual basin loop
		tout= tout + dtout
		endif
end do  !End main loop

close(11)
close(12)
close(20)
close(30)
close(21)
close(22)
close(31)
close(32)

print*, ''
print*, ''
print*, 'Flow through strait = ', strait_depth * strait_width * U_lat/1.e6, ' Sv'

end
!
! end of main
! --------------------------------------------------------------

! --------------------------------------------------------------

subroutine FTCS (imax, dx, dt, h, h1, S, D, w, bcl1,Db)

implicit none

common /bcs/ ibcl, ibcr, bcr, bcr1

integer :: i, j, n, ibcl, ibcr, imax
! integer :: info
real :: dx, dt, Db
real :: bcr, bcr1
real*16 :: bcl1, bcl
real*16:: h(0:imax+1)
real*16 :: S(imax)
real*16 :: w(imax)
real*16 :: h1(0:imax+1), D(0:imax+1)
real :: an, dn
real*16 :: Dph, Dmh

! advance solution for internal nodes
!
do i= 1,imax

   h1(i)= h(i) - (( w(i)*dt )/( 2.e0*dx )) * ( h(i+1) - h(i-1) ) + S(i)*dt

   Dph= ( D(i+1)+D(i) )/2.e0
   Dmh= ( D(i)+D(i-1) )/2.e0

   h1(i)= h1(i) + ( (Dph*dt)/(dx*dx) )*( h(i+1) - h(i) )
   h1(i)= h1(i) - ( (Dmh*dt)/(dx*dx) )*( h(i) - h(i-1) )

end do

! use bcs for external nodes
!
if (ibcl == 0) then

   h1(0)= bcl1

else if (ibcl == 1) then

   Dph= ( D(1)+D(0) )/2.e0
   h1(0)= h1(1) + (dx/Dph)*bcl1

endif

if (ibcr == 0) then

   h1(imax+1)= bcr1

else

! NOT TESTED
!
   Dmh= ( D(imax)+D(imax+1) )/2.e0
   h1(imax+1)= h1(imax) - (dx/Dmh)*bcr1

endif

return
end

! ............................................

real function density(T,S)

implicit none

real(16) :: T, S

! linear equation of state from Johnson et al. (2007) - CHECK 
!real, parameter :: r0= 1027.5e0, alpha= 2.e-4, beta= 8.e-4
!real, parameter :: T0= 5.e0, S0= 35.e0
!
! Wunsch (2015) writes r= r0*(1 - alpha*T + beta*S)
!  with r0= 1038 kg/m3, alpha= 0.5-3e-4 /C, beta= 0.8e-4
!  citing Thorpe (2005).
!real, parameter :: r0= 1038.e0, alpha= 20.e-5, beta= 0.8e-4
!
! Cessi et al. (2014)
!real, parameter :: r0= 1029.e0, alpha= 2.3e-4, beta= 7.5e-4
!
! Ivanov et al. (2002) and Anati (1997)
real, parameter :: r0= -15.654, alpha= 0.4309, beta= 0.936

!density= r0*(1.e0 - alpha*(T-T0) + beta*(S-S0))

! Cessi et al. (2014)
!density= r0*(1.e0 - alpha*T + beta*S)

! Ivanov et al. (2002) and Anati (1997)
density= 1.e3 + ( r0 - alpha*T + beta*S )

return
end

! ............................................
subroutine progress(t,maximum,minimum,p,n,outputnr)
  implicit none
  real ::maximum,minimum
  integer :: t
  integer:: p
  integer::k, n, outputnr
  character(len=75)::bar="???% |                                        | output time = ????/????    "
  write(unit=bar(1:3),fmt="(i3)") nint((100*(abs(t)-minimum)*(1/(maximum-minimum))))
  write(unit=bar(63:66),fmt="(i4)") n
  write(unit=bar(68:71),fmt="(i4)") outputnr
  do k=1, nint((40*(abs(t)-minimum)*(1/(maximum-minimum))))
    bar(6+k:6+k)="#"
  enddo
  ! print the progress bar.
  if (p==1) then
    bar="???% |                                        | output time = ????/????    "
  endif
  write(unit=6,fmt="(a1,a1,a75)",advance="no") '-',char(13), bar
  !print*, ((t-minimum)*(1/(maximum-minimum)))
  !(t-min)*(max/(max-min))
  return
end subroutine progress

! ............................................

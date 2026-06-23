! -----------------------------------------------------------------------------------------
! Description
! -----------------------------------------------------------------------------------------
! This program generates stochastic bridges for a two-dimensional stochastic process described by the following
! stochastic differential equations (SDEs):
!    dx/dt = Fx(x,y) + sqrt(D) Gx(x,y) xi_x(t),
!    dy/dt = Fy(x,y) + sqrt(D) Gy(x,y) xi_y(t),
! where xi_x(t) and xi_y(t) are independent Gaussian white noises. The equations are interpreted in the Itô sense. 
!
! A bridge is a trajectory conditioned to start at (x(0), y(0)) = (x0, y0) and end at (x(T), y(T)) = (xT, yT).
! The program generates M independent bridges and computes their corresponding transition times between 
! neighbourhoods around the initial and target states.
!
! -----------------------------------------------------------------------------------------
! Bridge parameters
! -----------------------------------------------------------------------------------------
! M      : number of bridges to generate
! dt     : integration time step
! tT     : bridge duration
! (xT,yT): target state
! (x0,y0): initial state
! C      : constant used in the acceptance probability
! R      : radius of the neighbourhoods around the initial and target states used to define transition times
!
! In this example, stochastic bridges are generated for the genetic toggle switch model. 
! They represent transitions from the undifferentiated state (x0,y0) to the differentiated state (xT,yT).
!
! -----------------------------------------------------------------------------------------
! SDEs parameters
! -----------------------------------------------------------------------------------------
! D     : noise intensity
! Fx, Fy: drift functions
! Gx, Gy: noise functions
! Additional model parameters may be required depending on the chosen system.
!
! In this example, Fx, Fy, Gx and Gy correspond to the genetic toggle switch model.
!
! -----------------------------------------------------------------------------------------
! Quasi-stationary distribution
! -----------------------------------------------------------------------------------------
! The bridge-sampling algorithm requires the quasi-stationary probability distribution PQS(x,y). 
! Here it is assumed to be available from a previous numerical computation and stored in a text file.
!
! [xmin,xmax] × [ymin,ymax]: spatial domain
! Nx, Ny                   : number of mesh points in x and y directions
! filename_PQS             : file containing PQS(x,y) on an Ny x Nx grid
!
! In this example, we use the file "PQS_geneticswitch.txt". 
! It corresponds to the toggle switch model under demographic noise with intensity D = 0.005.
! The distribution describes the metastable state (x0,y0) (differentiated state) from which we want to sample
! escape paths.
!
! -----------------------------------------------------------------------------------------
! Output
! -----------------------------------------------------------------------------------------
! bridges.txt       : sampled stochastic bridges
! transitiontime.txt: corresponding transition times 
!----------------------------------------------------------------------------------------------------------
program Bridges2D
	implicit double precision(a-h,o-z)
	!!! Bridge parameters
	integer, parameter :: M = 10
	double precision, parameter :: dt = 0.0005, tT = 30, xT = 0.0039216, yT = 1.99608, x0 = 1, y0 = 1, C = 0.5, R = 0.20
	!!! SDEs parameters 
	double precision, parameter :: D = 0.005
	double precision, parameter :: a = 1.0, b = 1.0, k = 1.0, S = 0.5, n = 4, Sn = S**n ! Parameters of the toggle switch model
	!!! Quasi-stationary distribution parameters
	double precision, parameter :: xmin = 0, xmax = 3, ymin = 0, ymax = 3
	integer, parameter :: Nx = 601, Ny = 601
	character(len=40) :: filename_PQS = "PQS_geneticswitch.txt" 
	!!! General definitions
	double precision, allocatable :: x_bridge(:), y_bridge(:), t_bridge(:), PQS(:,:), xvec(:), yvec(:)
	integer :: flagstart, flagend
	
	call cpu_time(t1)
	call dran_ini(12345) ! Initialization of the random number generator (program dranxor.f90, available in this repository)
	
	!-------------------------------------------------------------------------------------
	! Derived quantities and array allocation
	!-------------------------------------------------------------------------------------
	sqrtdt = sqrt(dt)
	sqrtD = sqrt(D)
	pi = 4.0d0 * datan(1.0d0)
	Nsteps = nint(tT / dt) + 1 ! Total number of temporal steps (including t = tT and t = 0)
    allocate(x_bridge(Nsteps), y_bridge(Nsteps), t_bridge(Nsteps))
	
	!-------------------------------------------------------------------------------------
	! Open output files
	!-------------------------------------------------------------------------------------
	open(unit = 1, file = "bridges.txt", status = "replace", action = "write")
	open(unit = 2, file = "transitiontimes.txt", status = "replace", action = "write")

	!-------------------------------------------------------------------------------------
	! Load the quasi-stationary distribution
	!-------------------------------------------------------------------------------------
 	call import_PQS(xmin, xmax, ymin, ymax, Nx, Ny, filename_PQS, PQS, xvec, yvec)

	!-------------------------------------------------------------------------------------
	! Generation of M independent stochastic bridges
	!-------------------------------------------------------------------------------------
	do ijk = 1, M
		print *, "######", ijk

		! Auxiliary quantities to compute the transition time
		flagstart = 0
		flagend = 0
		tstart = 0
		tend = 0
		
		!----------------------------------------------
		! Current state (x(t+dt), y(t+dt))
		!---------------------------------------------
		t = tT
		xtdt = xT
		ytdt = yT
		Pxtdtytdt = zxy_bilinear_interp(xtdt, ytdt, Nx, Ny, xvec, yvec, PQS)
		
		x_bridge(1) = xT 
		y_bridge(1) = yT 
		t_bridge(1) = tT 
		
		counter = 1 ! Number of performed time steps (including t = tT)

		do while (counter < Nsteps)
			!----------------------------------------------
			! Proposal for the state (x(t), y(t))
			!---------------------------------------------
			x = sample_gaussian(xtdt - Fx(xtdt, ytdt) * dt, Gx(xtdt, ytdt) * sqrtdt) ! Proposal for x(t)
			y = sample_gaussian(ytdt - Fy(xtdt, ytdt) * dt, Gy(xtdt, ytdt) * sqrtdt) ! Proposal for y(t)
			
			if (x.ge.xmin .and. x.le.xmax .and. y.ge.ymin .and. y.le.ymax) then ! Only evaluate proposals lying inside the domain
				!----------------------------------------------
				! Evaluation of the proposal
				!---------------------------------------------- 
				Pxy = zxy_bilinear_interp(x, y, Nx, Ny, xvec, yvec, PQS)
				Pratio = Pxy / Pxtdtytdt
				Gxratio = gaussian(xtdt, x + Fx(x, y) * dt, Gx(x, y) * sqrtdt) / gaussian(x, xtdt - Fx(xtdt, ytdt) * dt,  Gx(xtdt, ytdt) * sqrtdt)
				Gyratio = gaussian(ytdt, y + Fy(x, y) * dt, Gy(x, y) * sqrtdt) / gaussian(y, ytdt - Fy(xtdt, ytdt) * dt,  Gy(xtdt, ytdt) * sqrtdt)
				H = C * Gxratio * Gyratio * Pratio ! Acceptance probability
				if (H.ge.1) print *, "Warning: H >= 1. Try reducing C. H =", H
				if (dran_u().lt.H) then 
					!----------------------------------------------
					! Proposal accepted: (x(t),y(t)) = (x, y)
					!----------------------------------------------
					counter = counter + 1
					xtdt = x
					ytdt = y
					Pxtdtytdt = Pxy
					t = t - dt
					x_bridge(counter) = x
					y_bridge(counter) = y
					t_bridge(counter) = t	

					!----------------------------------------------
					! Check whether the trajectory has left the neighbourhood of the target state or entered the neighbourhood of the initial state
					!----------------------------------------------
					if (((x-xT)**2+(y-yT)**2.gt.R**2) .and. flagend.eq.0) then ! First exit from the neighbourhood of the target state
						tend = t
						flagend = 1
					end if

					if (((x-x0)**2+(y-y0)**2.lt.R**2) .and. flagstart.eq.0 .and. flagend.eq.1) then ! First entry into the neighbourhood of the initial state
						tstart = t
						flagstart = 1
					end if				
				end if	
			end if
		end do
		
		!----------------------------------------------
		! Export bridge data
		!----------------------------------------------
		if (ijk.eq.1) then  ! Bridge
			write(1, '( *(F12.6, 1X) )') t_bridge ! t_bridge is common to all bridges, so we write this vector only once
			write(1, '( *(F12.6, 1X) )') x_bridge
			write(1, '( *(F12.6, 1X) )') y_bridge
		else 
			write(1, '( *(F12.6, 1X) )') x_bridge
			write(1, '( *(F12.6, 1X) )') y_bridge
		end if
		
		if (flagstart.eq.1) then ! Transition time
			transitiontime = tend - tstart
			print *, "Transition time =", transitiontime
			write(2,*) transitiontime
		else
			print *, "Transition not completed"
		end if
		
	end do 
   
   	print *, "######"
   	call cpu_time(t2)
   	print *, "Simulation time =", t2 - t1
   	
contains

	double precision function Fx(x, y) ! Drift function in the x-direction
        double precision, intent(in) :: x, y
        Fx = a * x**n / (Sn + x**n) + b * Sn / (Sn + y**n) - k * x
	end function Fx
	
	double precision function Fy(x, y) ! Drift function in the y-direction
        double precision, intent(in) :: x, y
        Fy = a * y**n / (Sn + y**n) + b * Sn / (Sn + x**n) - k * y
	end function Fy

	double precision function Gx(x, y) ! Noise function in the x-direction
        double precision, intent(in) :: x, y
        Gx = sqrtD * sqrt(a * x**n / (Sn + x**n) + b * Sn / (Sn + y**n) + k * x)
	end function Gx
	
	double precision function Gy(x, y) ! Noise function in the y-direction
        double precision, intent(in) :: x, y
        Gy = sqrtD * sqrt(a * y**n / (Sn + y**n) + b * Sn / (Sn + x**n) + k * y)
	end function Gy

	double precision function gaussian(x, mu, sigma) 
        double precision, intent(in) :: x, mu, sigma
		gaussian = exp(-(x - mu)**2 / (2.0 * sigma**2)) / (sigma * sqrt(2.0 * pi))		
	end function gaussian
 	
	double precision function sample_gaussian(mu, sigma) 
	    ! Samples a random variable from a Gaussian distribution with mean mu and standard deviation sigma
		double precision, intent(in) :: mu, sigma
		sample_gaussian = sigma * dran_gbmw() + mu
	end function sample_gaussian
	
	double precision function zxy_bilinear_interp(x, y, Nx, Ny, xvec, yvec, zvec)
		! Bilinear interpolation of z(x,y)
		! zvec: function z evaluated at the 2D grid defined by xvec (Nx components) and yvec (Ny components).
	  	implicit double precision(a-h,o-z)
  		double precision, intent(in) :: x, y
  		integer, intent(in) :: Nx, Ny
  		double precision, dimension(Nx), intent(in) :: xvec
  		double precision, dimension(Ny), intent(in) :: yvec
  		double precision, dimension(Ny, Nx), intent(in) :: zvec

		if (x <= xvec(1)) then
			i = 1
		else if (x >= xvec(Nx)) then
			i = Nx - 1
		else
			do i = 1, Nx - 1
				if (xvec(i) <= x .and. x <= xvec(i + 1)) exit
			end do
		end if
	
		if (y <= yvec(1)) then
			j = 1
		else if (y >= yvec(Ny)) then
			j = Ny - 1
		else
			do j = 1, Ny - 1
				if (yvec(j) <= y .and. y <= yvec(j + 1)) exit
			end do
		end if
	
		! Values x1, x2, y1, y2 such that x1 <= x <= x2 and y1 <= y <= y2
		x1 = xvec(i) 
		x2 = xvec(i + 1)
		y1 = yvec(j)
		y2 = yvec(j + 1)
	
		! Funcion z evaluated at the four surrounding grid points
		z11 = zvec(j, i)
		z12 = zvec(j+1, i)
		z21 = zvec(j, i+1)
		z22 = zvec(j+1, i+1)
		
		dx = x2 - x1
		dy = y2 - y1
		
		zxy_bilinear_interp = 1 / (dx * dy) * (z11 * (x2 - x) * (y2 - y) + z21 * (x - x1) * (y2 - y) + z12 * (x2 - x) * (y - y1) + z22 * (x - x1) * (y - y1))
  	end function zxy_bilinear_interp
     
	subroutine import_PQS(xmin, xmax, ymin, ymax, Nx, Ny, filename_PQS, PQS, xvec, yvec)
		implicit double precision(a-h,o-z)
    	double precision, intent (in) :: xmin, xmax, ymin, ymax
    	integer, intent (in) :: Nx, Ny
    	double precision, allocatable, intent(out) :: PQS(:,:), xvec(:), yvec(:)
		character(len=*), intent(in) :: filename_PQS 

		allocate(PQS(Ny, Nx), xvec(Nx), yvec(Ny))
		
		dx = (xmax - xmin) / dble(Nx - 1) ! Spatial resolution in the x-direction 
    	dy = (ymax - ymin) / dble(Ny - 1) ! Spatial resolution in the y-direction 
    	
		open(unit = 10, file = filename_PQS, status = "old", action = "read")
    
    	do i = 1, Nx
        	xvec(i) = xmin + (i - 1) * dx
		end do
	
		do j = 1, Ny
			yvec(j) = ymin + (j - 1) * dy
		end do
			
		do j = 1, Ny
			read(10, *) PQS(j, :)
		end do
		
					
	end subroutine import_PQS


end program Bridges2D
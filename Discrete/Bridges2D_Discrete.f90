! -----------------------------------------------------------------------------------------
! Description
! -----------------------------------------------------------------------------------------
! This program generates stochastic bridges for a two-dimensional Markov jump process (nx, ny).
!
! A bridge is a trajectory conditioned to start at (nx(0), ny(0)) = (nx0, ny0) and end at (nx(T), ny(T)) = (nxT, nyT).
!
! The dynamics are simulated using the Gillespie algorithm with reversed transition rates, which are obtained as:
! W[(nx, ny)<-(nx', ny')] = W[(nx, ny)->(nx', ny')] * PQS(nx, ny) / PQS(nx', ny'),
! where PQS(nx, ny) is the quasi-stationary distribution. 
!
! The program generates M independent bridges and computes their corresponding transition times between 
! neighbourhoods around the initial and target states.
!
! -----------------------------------------------------------------------------------------
! Bridge parameters
! -----------------------------------------------------------------------------------------
! - M        : number of bridges to generate
! - tT       : bridge duration
! - (nxT,nyT): target state
! - (nx0,ny0): initial state
! - nxnylim  : size of the square neighbourhoods around the initial and target states used to define entry/exit
!   conditions for transition time measurements
!
! In this example, stochastic bridges are generated for the reaction-based genetic toggle switch model. 
! They represent transitions from the undifferentiated state (nx0,ny0) to the differentiated state (nxT,nyT).
!
! -----------------------------------------------------------------------------------------
! Model parameters
! -----------------------------------------------------------------------------------------
! - NN: system size (inverse noise intensity)
! - Reversed transition rates: they appear in lines 111-115 and 123-131, which must be adapted to the particular problem under study
! In this example, they correspond to the genetic toggle switch model. 
!
! -----------------------------------------------------------------------------------------
! Quasi-stationary distribution
! -----------------------------------------------------------------------------------------
! The bridge-sampling algorithm requires the quasi-stationary probability distribution PQS(nx,ny). 
! Here it is assumed to be available from a previous numerical computation and stored in a text file.
!
! - [0,nx_max] x [0,ny_max]: spatial domain
! - filename_PQS           : file containing PQS(nx,ny) on a uniform 2D grid
!
! In this example, we use the file "PQS_geneticswitch_Discrete.txt". 
! It corresponds to the toggle switch model with system size NN = 100.
! The distribution describes the metastable state (nx0,ny0) (differentiated state) from which we want to sample
! escape paths.
!
! -----------------------------------------------------------------------------------------
! Output
! -----------------------------------------------------------------------------------------
! bridges.txt       : sampled stochastic bridges
! transitiontime.txt: corresponding transition times 
!----------------------------------------------------------------------------------------------------------
program Bridges2D_Discrete
	implicit double precision(a-h,o-z)
	!!! Model parameters
	integer, parameter :: NN = 100
	double precision, parameter :: a = 1.0, b = 1.0, k = 1.0, S = 0.5, n = 4, Sn = S**n ! Parameters of the reaction-based toggle switch model
	!!! Bridge parameters
	integer, parameter :: M = 10, nxT = nint(0.0039216 * NN), nyT = nint(1.99608 * NN), nx0 = nint(1.0d0 * NN), ny0 = nint(1.0d0 * NN), nxnylim = nint(0.20 * NN)
	double precision, parameter :: tT = 30
	!!! Quasi-stationary distribution parameters
	integer, parameter :: nx_max = 1500, ny_max = 1500
	character(len=40) :: filename_PQS = "PQS_geneticswitch_Discrete.txt"
	!!! General definitions
	double precision :: PQS(0:nx_max, 0:ny_max)
	integer :: flagstart, flagend
	
	call cpu_time(t1)
	call dran_ini(12345) ! Initialization of the random number generator (program dranxor.f90, available in this repository)
	
	!-------------------------------------------------------------------------------------
	! Open output files
	!-------------------------------------------------------------------------------------
	open(unit = 1, file = "bridges.txt", status = "replace", action = "write")
	open(unit = 2, file = "transitiontimes.txt", status = "replace", action = "write")

 	!-------------------------------------------------------------------------------------
	! Load the quasi-stationary distribution
	!-------------------------------------------------------------------------------------
 	call import_PQS(PQS, nx_max, ny_max, filename_PQS)

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
		! Initialize trajectory at the final (target) state and evolve backward in time
		! This is the current state (nx(t+dt), ny(t+dt))
		!---------------------------------------------- 
		nx = nxT
		ny = nyT
		t = tT
		Pnxny = PQS(nx,ny)
				
		do while (t.gt.0)
			!----------------------------------------------
			! We sample the state (nx(t), ny(t)) using the Gillespie algorithm with reversed transition rates
			!---------------------------------------------- 
			
			! Model-specific reversed transition rates 
			W1 = k * (nx+1) * PQS(nx+1,ny) / Pnxny ! Rate W[(nx+1,ny)<-(nx,ny)]
			W2 = k * (ny+1) * PQS(nx,ny+1) / Pnxny ! Rate W[(nx,ny+1)<-(nx,ny)]
			W3 = NN * (a * (dble(nx-1)/NN)**n/(Sn+(dble(nx-1)/NN)**n) + b * Sn/(Sn+(dble(ny)/NN)**n)) * PQS(nx-1,ny) / Pnxny ! Rate W[(nx-1,ny)<-(nx,ny)]	
			W4 = NN * (a * (dble(ny-1)/NN)**n/(Sn+(dble(ny-1)/NN)**n) + b * Sn/(Sn+(dble(nx)/NN)**n)) * PQS(nx,ny-1) / Pnxny ! Rate W[(nx,ny-1)<-(nx,ny)]	
			W = W1 + W2 + W3 + W4 ! Total escape rate from current state

			! Time interval until the next jump
			dt = -log(dran_u()) / W ! dran_u(): generates a float random number in the interval [0,1)
			t = t - dt
			
			! State after the jump
			r = dran_u() * W
			if (r.lt.W1) then
				nx = nx + 1
			else if (r.lt.W1+W2) then
				ny = ny + 1
			else if (r.lt.W1+W2+W3) then
				nx = nx - 1
			else
				ny = ny - 1
			end if
			Pnxny = PQS(nx,ny)
						
			!----------------------------------------------
			! Check whether the trajectory has left the neighbourhood of the target state or entered the neighbourhood of the initial state
			!----------------------------------------------
			if ((nx.lt.nxT-nxnylim .or. nx.gt.nxT+nxnylim .or. ny.lt.nyT-nxnylim .or. ny.gt.nyT+nxnylim) .and. flagend.eq.0) then
				tend = t
				flagend = 1
			end if
				
			if ((nx.le.nx0+nxnylim .and. nx.ge.nx0-nxnylim .and. ny.le.ny0+nxnylim .and. ny.ge.ny0-nxnylim) .and. flagstart.eq.0 .and. flagend.eq.1) then
				tstart = t
				flagstart = 1
			end if
			
			!----------------------------------------------
			! Export bridge data
			! In discrete-state stochastic dynamics, jump times are random, so trajectories must be recorded event-by-event.
			!----------------------------------------------
			write(1,'(ES24.14, 3X, I6, 3X, I6)') t, nx, ny
			
		end do
		
		write(1,*) ! Blank line between trajectories 

		if (flagstart.eq.1) then 
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
     
	subroutine import_PQS(PQS, nx_max, ny_max, filename_PQS)
		implicit double precision(a-h,o-z)
    	integer, intent(in) :: nx_max, ny_max
    	double precision, intent(out) :: PQS(0:nx_max, 0:ny_max)
		character(len=*), intent(in) :: filename_PQS 
	
		open(unit = 10, file = filename_PQS, status = "old", action = "read")

    	do ny = 0, ny_max
    		read(10,*) (PQS(nx, ny), nx = 0, nx_max)
		end do
    	
    	close(10)
			
	end subroutine import_PQS

end program Bridges2D_Discrete
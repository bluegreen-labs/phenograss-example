module SANN
contains

subroutine SA(dataFile)
  use MODEL
  use FUNC
  implicit none
  
  ! SANN settings (temperature, cooling rate etc)
  real 						:: T, EPS, RT, FOPT,  FP, VMtmp, Ctmp
  integer 					:: NS, NT, IER, MAXEVL, IPRINT, NEPS

  ! SANN parameter ranges and associated files (dynamically sized)
  real, allocatable  		:: X(:), LB(:), UB(:), C(:), VM(:),&
								FSTAR(:),XOPT(:), XP(:)

  ! more parameter mess yay
  integer 					:: N, NACC,NOBDS, NFCNEV
  logical 					:: MAX, QUIT, QUIT_IT
  real 						:: F, P, PP, RATIO, RNUMBER,&
								acceptance_rate,rejection_rate
  integer 					:: NUP, NDOWN, NREJ, NNEW, LNOBDS,&
								H, I, J, M, ENOUT, PAROUT

  ! data
  real						:: dataFile(:,:,:)

  ! stuff needed to read the input data files
  character(len=200) 		:: line, filepath, sannPar, modRange
  integer 					:: io, count, nrsites
  integer 					:: nrPar, extPar
  
  ! random number generated related declarations
  integer 					:: rand_seed(18), seed_size
  integer, allocatable 		:: seed(:),  NACP(:)

  sannPar = "./parameters/SANN_parameters.txt"
  modRange= "./parameters/model_parameter_ranges.txt"

! ---------------------- IMPORT PARAMETERS -----------------------------

  ! read first argument, and determine the number of files in the 
  ! sites file
  count = command_argument_count()
  if (count < 1) then
	call exit(0)
  else 
	call get_command_argument(1, filepath)

	! read the number of sites in the file
	call LinesInFile(filepath,nrsites)
  end if

  ! import the SANN settings: Temperature / cooling rate etc
  open(1, file = trim(adjustl(sannPar))) 				
  i = 1
  sannpar_loop: do
  read(1,'(A)',iostat=io) line
	if (io < 0) then
			exit sannpar_loop 											
		else if ( index(line, "#") /= 0 ) then 						
		else if ( line == "" ) then 									
		else
			select case (i)
				case (1)
				 read(line,*) T 											
				case (2)
			         read(line,*) EPS 										
				case (3)
				 read(line,*) RT 										
				case (4)
				 read(line,*) NS 										
				case (5)
				 read(line,*) NT 									
				case (6)
				 read(line,*) MAXEVL 									
				case (7)
				 read(line,*) VMtmp 										
				case (8)
				 read(line,*) Ctmp 										
				case (9)
				 read(line,*) NEPS 										
				case (10)
				 read(line,*) IPRINT 									
			end select
			i = i + 1
	end if
  end do sannpar_loop
  close(1)

  ! get the number of parameters
  call LinesInFile(modRange,nrPar)
  
  ! after reading in parameters, determine how many there are
  N = nrPar !+ nrsites
    
  ! allocate memory to upper / lower boundaries and par vector
  allocate(LB(N),UB(N),X(N))
  
  ! but also allocate memory to the matrices holding the intermediate
  ! optimization results which should be of the same size!!
  allocate(XP(N),XOPT(N),C(N),VM(N),&
			FSTAR(NEPS),NACP(N)) 
  
  ! set step length for the different parameters
  do i = 1, N										
    VM(i) = VMtmp
	C(i) = Ctmp	
  end do

  ! set all vectors to 0
  LB = 0
  UB = 0
  X = 0
    
  ! read all parameters  
  open(1, file =trim(adjustl(modRange)))
  i = 1
  parameter_loop: do
  read(1,'(A)',iostat=io) line
	if (io < 0) then
			exit parameter_loop
		else if ( index(line, "#") /= 0 ) then 
		else if ( line == "" ) then
		else
		read(line,*) lb(i), ub(i)
		i = i + 1
	end if
  end do parameter_loop
  close(1)
  
  do i=1,nrPar
  	X(i) = (LB(i) + UB(i))/2
  end do

! ---------------------- IMPORT PARAMETERS - DONE ----------------------

! Print information on the input parameters as defined in the
! SANN_parameters file
	  write(*,1000) N, T, RT, EPS, NS, NT, NEPS, MAXEVL

1000 format(/,15x,'INITIATING SIMULATED ANNEALING'                   ,/,&
                5x ,'Number of parameters - no obs. operators       : ',I3/,&
                5x ,'Initial Temperature (T)                        : ',F8.4/,& 
                5x ,'Temperature reduction factor (RT)              : ',F8.4/,&
                5x ,'Error tolerance for termination (EPS)          : ',F8.4/,&
                5x ,'Number of Cycles (NS)                          : ',I10/,&
                5x ,'Itt. before Temp. reduction (NT)               : ',I10/,&
                5x ,'Nr. stable function values before exit (NEPS)  : ',I10/,&
                5x ,'Maximum function evaluations (MAXEVL)          : ',I10) 

  ! Generate random seed using machine state (CPU / clock)
  call random_seed() 													! initialize with system generated seed
  call random_seed(size=seed_size) 										! find out size of seed
  allocate(seed(seed_size))
  call random_seed(get=seed) 											! get system generated seed

  ! general flag to quit the routine, set processing
  QUIT 		= .false.
  QUIT_IT 	= .false.
  MAX 		= .false.

  ENOUT 	= 200
  PAROUT 	= 300

  ! set initial values
  NACC 		= 0
  NOBDS 	= 0
  NFCNEV 	= 0
  IER 		= 99

  do I = 1, N
    XOPT(I) = X(I)	! copy the start parameter values inot XOPT (optimal parameters)
    NACP(I) = 0   	! set number of accepted trials per parameter to 0	
  end do

  do I = 1, NEPS
    FSTAR(I) = 1.0D+20  ! set F* to a very high value (to be minimized)
  end do

  ! Give warning when temperature is negative, exit the SANN run
  ! by setting QUIT to TRUE
  if (T .LE. 0.0) then
     QUIT = .true.
     write(*,'(/,5x,''  THE INITIAL TEMPERATURE IS NOT POSITIVE. CHANGE VALUE!'',/)')
     IER = 3
  end if

  ! Give a warning when the parameters are out of bounds, set QUIT
  do I = 1, N
      if ((X(I) .gt. UB(I)) .or. (X(I) .lt. LB(I))) THEN
      	QUIT = .true.
	call PRT1
      end if
  end do
    
  XP = X
  do I = 1, 5
	!  Evaluate the function with input X and return value as F
	call cost(XP,dataFile,F)
        X = XP
  end do

  if(.not. MAX) F = -F
  NFCNEV = NFCNEV + 1
  FOPT = F
  FSTAR(1) = F
  if(IPRINT .eq. 1) call PRT2(MAX,X,F)

! write output to file
  open(ENOUT,file="summarizing_parameters.txt")

!-- MAIN LOOP
whileloop:  do

	NUP = 0
	NREJ = 0
	NNEW = 0
	NDOWN = 0
	LNOBDS = 0

	!-- should we not start the routine, as certain conditions are not met (see above)
	if (QUIT) exit whileloop

tempitterations: do M = 1, NT !-- loop which contains the statements before quenching temperature
updatesteps:	 do J = 1, NS !-- loop which contains statements before updateing step sizes 
parameters:	 do H = 1, N  !-- loop which goes through X parameters and updates them

				!-- we update the parameters to be evaluated
				!-- uniform random selection within the parameter range 
				do I = 1, N
				
	 				!-- for every parameter tweak slightly
					if (I .eq. H) then
						call random_number(RNUMBER)
						XP(I) = X(I) + (RNUMBER*2.-1.) * VM(I)
					else
						XP(I) = X(I)
					end if
					
					!-- if out of bound constrain
					if ((XP(I) .lt. LB(I)) .or. (XP(I) .gt. UB(I))) then
					call random_number(RNUMBER)
					XP(I) = LB(I) + (UB(I) - LB(I))*RNUMBER
                    			LNOBDS = LNOBDS + 1
                    			NOBDS = NOBDS + 1
					if(IPRINT .eq. 3) call PRT3(MAX,XP,X,FP,F)
					end if
					
				end do

				!-- evaluate the cost function with the new parameters
				!-- write results to a variable FP
				call cost(XP,dataFile,FP)

		                if (.not. MAX) FP = -FP
		                NFCNEV = NFCNEV + 1
		                if(IPRINT .eq. 3) call PRT4(MAX,XP,X,FP,F)

				!-- check if the maximum number of function evaluations
				!-- has been exceeded, if so flag as last itteration and exit		
				if (NFCNEV .ge. MAXEVL) then
					call PRT5
					if (.not. MAX) FOPT = -FOPT
					IER = 1
					QUIT_IT = .true.	!-- quit the loop when maximum eval exceeded, write last
					exit tempitterations	!-- parameters to file
                end if
				
				!-- if the new model error value is lower than the previous one
				!-- update parameter array to reflect this change 
				if(FP .ge. F) then
					if(IPRINT .eq. 3) then
                    	write(*,'(''  POINT ACCEPTED'')')
					end if
						X = XP						!-- update parameters
                  		F = FP 						!-- set new error minimum 
                  		NACC = NACC + 1 			!-- add a successful trial to the counter (global)
                  		NACP(H) = NACP(H) + 1 		!-- add a successful trial to the counter for this parameter
                  		NDOWN = NDOWN + 1			!-- mark move as down

					!-- if the new value is smaller than ALL previous ones
					!-- update parameters to reflect this change (set global optimum)
					if (FP .gt. FOPT) then 	!-- minimized value is lower than lowest on record

						if(IPRINT .eq. 3) then
							write(*,'(''  NEW OPTIMUM'')')
						end if
						
						XOPT = XP
						FOPT = FP
						NNEW = NNEW + 1

					end if

				!-- this is the METROPOLIS section
				!-- which allows random 'uphill' moves
				!-- in our 'downhill' optimization... 
				else	
					P = exp((FP - F)/T) 		!-- determine metropolis 'cut'
					call random_number(PP)		!-- random number
					if (PP .lt. P) THEN			!-- compare and decide
						if (IPRINT .eq. 3) call PRT6(MAX)
	
						X = XP						!-- update the parameter even if it is 'not optimal' 
						F = FP						!-- set new minimum / maximum
						NACC = NACC + 1				!-- add a successful trial to the counter (global)
						NACP(H) = NACP(H) + 1		!-- add a successful trial to the counter for this parameter
						NUP = NUP + 1				!-- mark move as going up, away from minimum
					else
						NREJ = NREJ + 1				!-- if the move is not accepted, do nothing, mark as rejected
						if (IPRINT .eq. 3) call PRT7(MAX)
					end if
				end if
				write(ENOUT,*) T,FOPT,X
					end do parameters
				end do updatesteps
				
				!-- dump everything to a line
				if (IPRINT .eq. 5) write(*,*) T,FOPT,F,X

			!--  Adjust VM so that approximately half of all evaluations are accepted.
         		do I = 1, N
            			RATIO = dfloat(NACP(I)) /dfloat(NS)
            			acceptance_rate = .5
            			rejection_rate = 1 - acceptance_rate
            			if (RATIO .gt. acceptance_rate) then
               				VM(I) = VM(I)*(1. + C(I) * (RATIO - acceptance_rate)/rejection_rate)
            			else if (RATIO .lt. rejection_rate) then
               				VM(I) = VM(I)/(1. + C(I) * (rejection_rate - RATIO)/rejection_rate)
            			end if
				if (VM(I) .gt. (UB(I)-LB(I))) then
					VM(I) = UB(I) - LB(I)
				end if
			end do

			if (IPRINT .eq. 2) then
				call PRT8(VM,XOPT,X)
			end if

			!-- reset acceptance values
			do I = 1, N
				NACP(I) = 0
			end do

		 end do tempitterations

	! provide feedback on this temperature itteration
	if (IPRINT .eq. 1) then
		call PRT9(max,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
	end if

	!-- Check termination criteria.
	QUIT = .false.
	FSTAR(1) = F
	if ((FOPT - FSTAR(1)) .LE. EPS) QUIT = .true.
	do I = 1, NEPS
		if (abs(F - FSTAR(I)) .gt. EPS) QUIT = .false.
	end do

	!-- Terminate SA if appropriate.									
	if (QUIT .or. QUIT_IT) then
		
		X = XOPT

		IER = 0
		if (.not. max) FOPT = -FOPT
		
		!-- give some feedback on termination criteria
		call PRT10(NFCNEV,NACC,NOBDS,T,F)

		!-- write everything to file !!!
		open(1,file="./parameters/optimized_model_parameters.txt")
	    write(*,*) "Optimal Parameters:"
		do i=1,size(X)
			write(*,*) X(i)
			write(1,*) X(i)
		end do
		close(1)

		exit whileloop !-- quit processing if criteria met
	end if

	!-- If termination criteria is not met, prepare for another loop.
	T = RT*T !-- adjust temperature

	do I = NEPS, 2, -1
		FSTAR(I) = FSTAR(I-1)
	end do

	F = FOPT
	X = XOPT

  end do whileloop
  close(ENOUT)
  
!-- END MAIN LOOP
end subroutine SA


!-- These are plotting subroutines for verbose output of results
subroutine PRT1
      write(*,'(5x,''  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS '',/,&
               5x,''  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY'',/,&
               5x,''  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  '',/,&
               5x,''  LB(I) .lt. X(I) .lt. UB(I), I = 1, N. ''/)')
end subroutine PRT1

! 
subroutine PRT2(max,X,F)
      real ::  X(:), F
      logical  max

      write(*,'(5x,''  '')')
      call PRTVEC(X,'INITIAL X')
      if (max) then
         write(*,'(5x,''  INITIAL F: '',/, G25.18,/)') F
      else
         write(*,'(5x,''  INITIAL F: '',/, G25.18,/)') -F
      end if
end subroutine PRT2


subroutine PRT3(max,XP,X,FP,F)
      real ::  XP(:), X(:), FP, F
      logical  max

      write(*,'(5x,''  '')')
      call PRTVEC(X,'CURRENT X')
      if (max) then
         write(*,'(5x,''  CURRENT F: '',G25.18)') F
      else
         write(*,'(5x,''  CURRENT F: '',G25.18)') -F
      end if
      call PRTVEC(XP,'TRIAL X')
      write(*,'(5x,''  POINT REJECTED SINCE OUT OF BOUNDS'')')
end subroutine PRT3


subroutine PRT4(max,XP,X,FP,F)
      real ::  XP(:), X(:), FP, F
      logical  max

      write(*,'(5x,''  '')')
      call PRTVEC(X,'CURRENT X')
      if (max) then
         write(*,'(5x,''  CURRENT F: '',G25.18)') F
         call PRTVEC(XP,'TRIAL X')
         write(*,'(5x,''  RESULTING F: '',G25.18)') FP
      else
         write(*,'(5x,''  CURRENT F: '',G25.18)') -F
         call PRTVEC(XP,'TRIAL X')
         write(*,'(5x,''  RESULTING F: '',G25.18)') -FP
      end if
end subroutine PRT4

subroutine PRT5
      write(*,'(5x,''  TOO MANY FUNCTION EVALUATIONS; CONSIDER '',/,&
                5x,''  INCREASING MAXEVL OR EPS, OR DECREASING '',/,&
                5x,''  NT OR RT. THESE RESULTS ARE LIKELY TO BE '',/,&
                5x,''  POOR.'',/)')
end subroutine PRT5

subroutine PRT6(max)
      logical  max

      if (max) then
         write(*,'(5x,''  THOUGH LOWER, POINT ACCEPTED'')')
      else
         write(*,'(5x,''  THOUGH HIGHER, POINT ACCEPTED'')')
      end if
end subroutine PRT6


subroutine PRT7(max)
      logical  max

      if (max) then
         write(*,'(''  LOWER POINT REJECTED'')')
      else
         write(*,'(''  HIGHER POINT REJECTED'')')
      end if
end subroutine PRT7

subroutine PRT8(VM,XOPT,X)
      real ::  VM(:), XOPT(:), X(:)

      write(*,'(/,5x,'' INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT'',/)')
      call PRTVEC(VM,'NEW STEP LENGTH (VM)')
      call PRTVEC(XOPT,'CURRENT OPTIMAL X')
      call PRTVEC(X,'CURRENT X')
      write(*,'(5x,'' '')')
end subroutine PRT8


subroutine PRT9(max,T,XOPT,VM,FOPT,NUP,NDOWN,NREJ,LNOBDS,NNEW)
      real 			::  XOPT(:), VM(:), T, FOPT
      integer  	:: NUP, NDOWN, NREJ, LNOBDS, NNEW, TOTMOV
      logical  	:: max

      TOTMOV = NUP + NDOWN + NREJ
         write(*,'(/,'' INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION'',/)')
	     write(*,'(''  CURRENT TEMPERATURE         :  '',G12.5)') T
      if (max) then
         write(*,'(''  MAX FUNCTION VALUE SO FAR   :  '',G25.18)') FOPT
         write(*,'(''  TOTAL MOVES                 :  '',I8)') TOTMOV
         write(*,'(''     UPHILL                   :  '',I8)') NUP
         write(*,'(''     ACCEPTED DOWNHILL        :  '',I8)') NDOWN
         write(*,'(''     REJECTED DOWNHILL        :  '',I8)') NREJ
         write(*,'(''  OUT OF BOUNDS TRIALS        :  '',I8)') LNOBDS
         write(*,'(''  NEW MAXIMA THIS TEMPERATURE :  '',I8)') NNEW
      else
         write(*,'(''  MIN FUNCTION VALUE SO FAR   :  '',G25.18)') -FOPT
         write(*,'(''  TOTAL MOVES                 :  '',I8)') TOTMOV
         write(*,'(''     DOWNHILL                 :  '',I8)')  NUP
         write(*,'(''     ACCEPTED UPHILL          :  '',I8)')  NDOWN
         write(*,'(''     REJECTED UPHILL          :  '',I8)')  NREJ
         write(*,'(''  TRIALS OUT OF BOUNDS        :  '',I8)')  LNOBDS
         write(*,'(''  NEW MINIMA THIS TEMPERATURE :  '',I8)')  NNEW
      end if
      call PRTVEC(XOPT,'CURRENT OPTIMAL PARAMETERS')
      call PRTVEC(VM,'STEP LENGTH (VM)')
      write(*,'('' '')')
end subroutine PRT9

subroutine PRT10(NFCNEV,NACC,NOBDS,T, F)
	  integer 		:: NFCNEV, NACC, NOBDS
	  real 			:: T, FSTAR
	  write(*,1001) NFCNEV, NACC, NOBDS, T, F
1001  format(/,10x ,'SANN ACHIEVED TERMINATION CRITERIA. IER = 0.',/,&
              5x ,' Number of functions evaluated   :      ',I10 ,/,&
              5x ,' Number of accepted functions    :      ',I10 ,/,&
              5x ,' Number of out of bound functions:      ',I10 ,/,&
              5x ,' Final Temperature               :      ',G20.13,/,&
              5x ,' Final Model Error               :      ',G20.13,/)
end subroutine PRT10

subroutine PRTVEC(VECTOR,NAME)
  !  This subroutine prints the double precision vector named VECTOR.
  !  Elements 1 thru NCOLS will be printed. NAME is a character variable
  !  that describes VECTOR. Note that if NAME is given in the call to
  !  PDA_PRTVEC, it must be enclosed in quotes. If there are more than 10
  !  elements in VECTOR, 10 elements will be printed on each line.

  real :: VECTOR(:)
  integer ::  NCOLS,I
  character *(*) NAME
  
  NCOLS = size(VECTOR)
  write(*,1001) NAME
  if (NCOLS .gt. 10) then
      LINES = INT(NCOLS/10.)

   do 100, I = 1, LINES
          LL = 10*(I - 1)
        write(*,1000) (VECTOR(J),J = 1+LL, 10+LL)
  100  continue
         write(*,1000) (VECTOR(J),J = 11+LL, NCOLS)
      else
         write(*,1000) (VECTOR(J),J = 1, NCOLS)
      end if

 1000 format( 10(G12.5,1X))
 1001 format(/,5X,A)

      return
end subroutine PRTVEC

end module SANN

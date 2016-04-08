module FUNC
! these are stand alone functions which only return values or data
! and don't rely on global parameters all input should be provided

contains
subroutine LinesInFile(dataFile,nrLines)
  implicit none

  ! assumed length (len=*) is required within subroutine
  character(len=*), intent(in)  :: dataFile
  integer, intent(out)          :: nrLines
  character(len=200)            :: line
  integer :: io

  ! open file
  open(1, file=trim(adjustl(dataFile)))

  ! initiate line number counter
  nrLines = 0

  ! read line by line and count valid lines
  read_loop: do
  read(1,'(A)',iostat=io) line
	if (io < 0) then
			exit read_loop
		else if ( index(line, "#") /= 0 ) then
		else if ( line == "" ) then
		else
			! count lines
			nrlines = nrLines + 1
	end if
  end do read_loop
  close(1)
end subroutine LinesInFile

! Calculates Potential Evapotranspiration on a daily timestep
subroutine PET(Tmax,Tmin,DOY,latitude,ET)
  implicit none
  real, intent(in) 	:: latitude
  real, intent(in) 	:: DOY(:)
  real, intent(in) 	:: Tmax(:), Tmin(:)
  real 					    :: R(size(DOY)), L(size(DOY))
  real, intent(out) 	:: ET(size(DOY))

  ! calculate daily radiation MJ m-2 day-1
  call Ra(DOY,latitude,R)

  ! conversion factor to translate radiation into
  ! evapotransipiration equivalents
  L = 0.0864 * (28.9 - 0.0028 * ((Tmax + Tmin)/2))

  ! calculate PET based upon Hargreaves equation
  ! but with (DVWK) daily Tmean dependent conversion factors L
  ET = 0.0023 * R/L * sqrt(Tmax - Tmin) * ( ((Tmax + Tmin)/2) + 17.8)

  ! set negative values to 0
  where (ET <= -0.0)
    ET = 0
  end where

end subroutine PET

! Calculates total extraterrestrial radiation, Ra on a dialy timestep
subroutine Ra(DOY,latitude,R)
  implicit none

  real, intent(in) :: DOY(:)
  real, intent(in) :: latitude
  real, intent(out) :: R(size(DOY)) ! calculated radiation vector
  real :: GSC												! Global Solar Constant
  real :: Omega	(size(DOY))					! latitude in radians
  real :: Delta	(size(DOY))					! solar declination in radians
  real :: Ws(size(DOY))
  real :: Dr(size(DOY))
  real :: pi

  ! asign pi
  pi=4.D0*datan(1.D0)

  ! convert latitude (degrees) to radians
  Omega = latitude * ( pi / 180 )

  ! solar constant [MJ m-2 min-1]
  GSC = 0.082

  ! inverse distance to the Earth-Sun
  Dr = 1 + 0.033 * cos( 2. * pi / 365 * DOY)

  ! solar declination [rad]
  Delta = 0.409 * sin ( 2. / 365 * pi * DOY - 1.39 )

  ! sunset hour angle [rad]
  Ws = acos(-1*tan(Omega)*tan(Delta))

  ! extraterrestrial radiation Ra [MJ m-2 day-1]
  R = (24*60)/pi * GSC * Dr * (Ws*sin(Omega)*sin(Delta) &
		+ cos(Omega)*cos(Delta)*sin(Ws))

end subroutine Ra

! Calculates total extraterrestrial radiation, Ra on a dialy timestep
subroutine DayLength(DOY,latitude,DL)
  implicit none

  real, intent(in) 		:: DOY(:)			     	! day of year
  real, intent(in) 		:: latitude				  ! latitude of a site
  real, intent(out) 	:: DL(size(DOY))		! calculated radiation vector
  real 					      :: Omega(size(DOY))	! latitude in radians
  real 					      :: Phi(size(DOY))		! solar declination in radians
  real 					      :: pi,p
  integer 				    :: i
  logical				      :: l

  ! set daylength coefficient
  p = 0.0

  ! asign pi
  pi=4.D0*datan(1.D0)

  ! Forsythe et al. 1995 eq. 1
  Omega = 0.2163108 + 2 * atan(0.9671396 * tan (0.00860 * (DOY - 186)))

  ! eq. 2
  Phi = asin(0.39795* cos(Omega))

  ! eq. 3 / returns daylength D
  DL = 24 - 24 / pi * acos((sin(p*pi/180)+sin(latitude*pi/180)*sin(Phi)) &
							/(cos(latitude*pi/180)*cos(Phi)))

  ! Normal output return NaN values, fill these in with either 0 or
  ! 24 for high / low latitudes
  do i=1,size(DL)
    l = DL(i-1) > 20
    if ( l  .and.  isnan(DL(i)) ) then
      DL(i) = 24
    end if
    if ( isnan(DL(i)) ) then
      DL(i-1) = 0
    end if
  end do
end subroutine DayLength

! number of space-separated items in a line
integer function ntokens(line)
	character,intent(in):: line*(*)
	integer i, n, toks
	i = 1;
	n = len_trim(line)
	toks = 0
	ntokens = 0
	do while(i <= n)
	   do while(line(i:i) == ' ')
		 i = i + 1
		 if (n < i) return
	   enddo
	   toks = toks + 1
	   ntokens = toks
	   do
		 i = i + 1
		 if (n < i) return
		 if (line(i:i) == ' ') exit
	   enddo
	enddo
end function ntokens

end module

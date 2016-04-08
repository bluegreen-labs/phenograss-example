module INC

  ! Everything defined here is shared between the subroutines
  ! and the main program. These variables are limited to what is
  ! strictly needed. All other variables are local variables within
  ! subroutines.
  ! TODO: properly implement dealllocate()

  ! counter increment variables
  integer 							:: i,j,r,k

  ! define variables for subroutine ImportDataFiles
  ! define 3-dim. array to hold all model data, pass this around
  real, allocatable     :: dataFile(:,:,:)
  integer               :: nrsites, nrfiles, nrDataPoints,&
	                         nrPar,argcount,nrParSet
  character(len=300)    :: fileindex,outdir
  real, allocatable     :: par(:,:)
  character(len=300)    :: modPar
  character(len=300)    :: modelrun

  contains
  subroutine importData
  use FUNC
  implicit none

  ! file reading variables, locally defined not needed globally
  integer               :: isite, ilat, io, iWCAP, iWP, count
  character(len=300)    :: filepath, tmpchar
  character(len=200)    :: line

  integer, allocatable  :: dataLength(:)
  real, allocatable     :: latitude(:),WCAP(:),WP(:)
  character(len=300),allocatable :: siteName(:), siteData(:)
  logical               :: mask(size(dataFile,1))
  integer,allocatable   :: doy(:)
  character(len=200)    :: empty, natrap

  ! read first argument -- file listing sites to process
  argcount = command_argument_count()
  if (argcount < 1) then
    call exit(0)
  else
    call get_command_argument(1, filepath)
    call LinesInFile(filepath,nrsites)
		modelrun = "modelrun"
		outdir="./output/"
	  modPar = "./parameters/optimized_model_parameters.txt"
  end if

  ! allocate the array for the sites and paths to the data
  allocate(siteName(nrsites),siteData(nrsites),dataLength(nrsites),&
			latitude(nrsites),WCAP(nrsites),WP(nrsites))

  ! read data location, sites file
  open(10, file = trim(adjustl(filepath)))
  I = 1
  sites_loop: do
  read(10,'(A)',iostat=io) line
	if (io < 0) then
			exit sites_loop
		else if ( index(line, "#") /= 0 ) then
		else if ( line == "" ) then
		else
		read(line,*) siteData(I)
		I = I + 1
	end if
  end do sites_loop
  close(10)

  ! For every site get the maximum length of the data file
  do i=1,size(siteName)
		call LinesInFile(siteData(i),dataLength(i))
  end do

  ! get the overall maximum length of the datafile
  ! this value is needed in case of uneven input files
  nrDataPoints = maxval(dataLength)

  ! allocate 3D array to contain the data with rows
  ! equal to nrDataPoints and columns equal to variables
  ! layers equal to # sites

  allocate(dataFile(nrDataPoints,18,nrsites))
  dataFile = 0 ! set everything to zero

  ! define some counters for the loop
  ! NOTE: I use string matching to get to the site
  isite=1
  ilat=1
  iWCAP=1
  iWP=1

  do j=1,nrsites

	  ! open site data file
	  open(1, file = trim(adjustl(siteData(j))))
	  i = 1
	  data_loop: do
	  read(1,'(A)',iostat=io) line
	  if (io < 0) then
			exit data_loop ! exit read loop if at end of file
		else if ( index(line,"# Site:") /= 0 ) then
			read(line,*) tmpchar, tmpchar, siteName(isite)
			isite = isite + 1
		else if ( index(line,"# Lat:") /= 0 ) then
			read(line,*) tmpchar, tmpchar, latitude(ilat)
			ilat = ilat + 1
		else if ( index(line,"# WCAP:") /= 0 ) then
			read(line,*) tmpchar, tmpchar, WCAP(iWCAP)
			iWCAP = iWCAP + 1
		else if ( index(line,"# WP:") /= 0 ) then
			read(line,*) tmpchar, tmpchar, WP(iWP)
			iWP = iWP + 1

		! check the first character and if == # skip line
		else if ( index(line, "#") /= 0 ) then
		else if ( line == "" ) then
		else
			read(line,*) empty, empty, empty, empty, natrap
			if ( trim(natrap) == "NA" ) then
				continue
			else
				! year(1), doy(2), validation(3), gcc(4), tmax(5), tmin(6), precip(8)
				! remaining columns remain blank for now -> contain model output
				! all data is passed around in this array format to make indexing easy
				! this is not optimal but works relatively well and fast in serial mode
				read(line,*) dataFile(i,1,j), dataFile(i,2,j),&
				dataFile(i,3,j),dataFile(i,4,j), dataFile(i,5,j),&
				dataFile(i,6,j), dataFile(i,8,j)
			end if
			  i = i + 1
	    end if
    end do data_loop
    close(1)

    ! calculate mean temperature
    dataFile(:,7,j) = (dataFile(:,5,j)+dataFile(:,6,j))/2

    ! takes tmax, tmin, doy (convert to int), latitude as input
    ! returns PET
	  call PET(dataFile(:,5,j),dataFile(:,6,j),dataFile(:,2,j),&
			   latitude(j),dataFile(:,9,j))

	  ! calculate daily radiation
	  call Ra(dataFile(:,2,j),latitude(j),dataFile(:,10,j))

	  ! mean annual precip
	  dataFile(:,14,j) = sum(dataFile(:,8,j)) / (maxval(dataFile(:,1,j))- minval(dataFile(:,1,j)) + 1)

	  ! pass along the AWC for each site
	  dataFile(:,15,j) = WCAP(j)

	  ! pass along the WP for each site
	  dataFile(:,16,j) = WP(j)

	  ! running mean temperature
	  do i=1,size(dataFile,1)
		  dataFile(i,17,j) = sum(dataFile((i-15):i,7,j)) / 15
	  end do
  end do
end subroutine importData

end module INC

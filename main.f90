! Defines the main program setup, which will call upon
! different subroutines for flexibility.

program main

  use, intrinsic :: iso_fortran_env ! subroutines to deal with i/o
  use FUNC 			    ! load helper functions
  use INC                           ! setup of the model
  use MODEL		            ! main model
  use SANN

  implicit none
  character(len=10000) 	:: customfilename,tmpchar
  character(len=10000)	:: line
  CHARACTER(LEN=300)  	:: Format
  call importData		    ! import data for all sites (if multiple)


  ! optimize the parameters if requested, otherwise read the parameter
  ! file and run the model using these parameters
  if (optim == "o") then
		! call the optimization routine
		call SA(dataFile)
  end if

  ! Read the optimized
  ! parameters for plotting
  call LinesInFile(modPar,nrPar)

  ! count the number of columns in a file
  open(unit=10,file=trim(adjustl(modPar)))
  !reading in data file
  read(10,'(A)') line	   ! reads first line to determine # of columns
  nrParSet = ntokens(line) ! calculate how many columns are in the file
  close(10)

  allocate(par(nrPar,nrParSet))
  open(1,file = trim(adjustl(modPar)))
  do i=1,nrPar
    read(1,*) par(i,:)
  end do
  close(1)

  ! run the model with the above optimized or loaded parameters
  ! if multiple parameter sets are given these are calculated
  ! and assigned to different output files (for cross validation)
  do r=1,nrparset
    do i=1,nrsites
      ! call the model and process a site for a given parameter set
      call fcover(par(:,r),dataFile,i)
    end do

    ! write the r parameter set integer to temporary character
    write (tmpchar, "(I2.2)") r
    customfilename = trim(adjustl(outdir))//"/"//trim(adjustl(modelrun))//"_"//trim(adjustl(tmpchar))//"_fCover.txt"
    open(1,file=customfilename)
    do k=1,nrDataPoints
      write(tmpchar,"(I4)") nrsites
      Format = "("//trim(adjustl(tmpchar))//"(1X,F7.6))" ! space delimited
      write(1,Format) (dataFile(k,11,j), j=1,nrsites)
    end do
    close(1)
   
    ! clear all previous output from the dataFile array
    dataFile(:,11,:) = 0

  end do

  ! write input data to file for plotting
  open(1,file="./output/precip.txt")
  do i=1,nrDataPoints
    write(1,*) (dataFile(i,8,j), j=1,nrsites)
  end do
  close(1)

  open(1,file="./output/MAP.txt")
  do i=1,nrDataPoints
    write(1,*) (dataFile(i,14,j), j=1,nrsites)
  end do
  close(1)

  open(1,file="./output/Gcc.txt")
  do i=1,nrDataPoints
    write(1,*) (dataFile(i,4,j), j=1,nrsites)
  end do
  close(1)

end program main   !-- end of main program

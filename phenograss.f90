module MODEL

contains
  !------------- FRACTIONAL VEGETATION COVER MODEL -----------------------
  subroutine fcover(par,dataFile,site)
    implicit none
    real					      :: dataFile(:,:,:) ! this is the data matrix
    real, intent(in) 		:: par(:)
    real					      :: precip(size(dataFile,1)),&
    evap(size(dataFile,1))
    real				 	      :: V(size(dataFile,1))
    real					      :: W(size(dataFile,1)),Dt(size(dataFile,1)),&
    Sd(size(dataFile,1)),Ra(size(dataFile,1)),g2(size(dataFile,1))
    real					      :: T(size(dataFile,1)),Tm(size(dataFile,1))
    real 					      :: b1,b2,b3,b4,b5,Rmax,a,b
    integer 				    :: L,d
    real 					      :: Wcap, Wstart,Dtl,Dtl1,Topt,Wp,dor
    integer				      :: site,nr_sites
    integer 				    :: i,values,nrPar
    real 					      :: Vmin, Vmax, Tmin
    real 					      :: intercept, slope, g, Tmax, Phmin,Phmax, m, c

    ! how many sites are we calculating
    nr_sites = size(dataFile,3)

    ! asign human readable drivers (extract from dataFile array)
    precip  = dataFile(:,8,site)											! precipitation
    evap    = dataFile(:,9,site)											! potential evapotranspiration
    T       = dataFile(:,6,site)											! mean temperature
    Ra      = dataFile(:,10,site)											! TOA radiation MJ m-2 s-1
    Tm      = dataFile(:,17,site)											! mean temperature vector
    Wcap    = dataFile(1,15,site)											! field capacity
    Wp      = dataFile(1,16,site)											! wilting point

    ! assign parameters
    b1      = par(1)
    b2 		  = par(2)
    b3 		  = par(3)
    b4 		  = par(4)
    L       = par(5)
    Phmin   = par(6)
    slope		= par(7)
    Topt		= par(8)
    Phmax		= par(9)

    ! Maximum temperature of the growth response curve
    ! corresponding to sensible literature values
    Tmin = 0
    Tmax 	= 45

    ! assign constant values to declared variables Vmin/Vmax
    Vmin 	= 0.001	! set to a small value but not 0
    Vmax 	= 1.		! set to 100% (GCC values are scaled between 0-1
    d     = 0		 	! decay flag

    ! assign values to both the W and V arrays
    ! these are initial conditions that influence the final
    ! results to a great extent!!
    W		     = 0.
    Wstart   = 0.
    V 	     = 0.001
    Sd 	     = 0.
    b1 	     = Wp
    m	       = 3600

    ! This loop takes the array of values
    ! above and processes the values sequentially,
    ! on a daily timestep.
    values = size(precip)-1
    do i=1,values ! -- main loop

      ! first if statement traps initial conditions
      ! where the values would be out of bound (looking
      ! at place such as -1)
      ! within the statement define Dtl and Dtl1
      ! state of the soil water at time day and day-1
      ! including a lag of L days
      if ( i - L - 1 < 0 ) then
        Dt(i) = max(0.0,W(i) - b1)
        Dtl = Wstart
        Dtl1 = Wstart
      else
        Dt(i) = max(0.0,W(i) - b1)
        Dtl = Dt(i-L)
        Dtl1 = Dt(i-L-1)
      end if

      ! if there is more precipitation (SWC)
      ! compared to the previous day then decay = 0
      ! otherwise the decay parameter is set to 1
      ! and senescence sets in
      if ( Dtl > Dtl1 ) then
        d = 0
      else
        d = 1
      end if

      ! Temperature response function
      ! and capture out of bound values
      ! for optimization
      g = ((Tmax - Tm(i))/(Tmax - Topt)) * (((Tm(i) - Tmin) / (Topt - Tmin)) ** (Topt/(Tmax - Topt)))

      ! Unfavourable conditions (too cold / hot) can return NA values
      ! which corrupt further calculations, trap these
      ! and set to growth 0 instead of NA (causing errors)
      if ( isnan(g) ) then
        g = 0
      end if

      ! set dormancy conditions based upon the
      ! available radiation
      dor = (Ra(i) - Phmin)/(Phmax - Phmin)
      if ( Ra(i) >= Phmax ) then
        dor = 1
      end if
      if ( Ra(i) <= Phmin ) then
        dor = 0
      end if

      ! SOIL WATER CONTENT SECTION:
      W(i+1) = W(i) + precip(i) - (1-V(i)) * ((( Dt(i)/(Wcap - b1)))**2) * evap(i) - g * b4 * V(i) * Dt(i)

      ! don't allow negative SWC values
      W(i+1) = max(0.0,min(Wcap,W(i+1)))

      ! VEGETATION GROWTH SECTION:
      V(i+1) = V(i) + g * dor * b2 * Dtl * (1 - V(i)/Vmax) - d * b3 * V(i) * (1-V(i))

      ! do not allow negative vegetation cover or values > 1
      V(i+1) = max(Vmin,min(Vmax,V(i+1)))

    end do ! -- end main loop

    ! put data back into the array
    ! packaging everything neatly together for
    ! writing to file
    dataFile(:,11,site) = V
    dataFile(:,7,site) = g2

  end subroutine fcover

!------------------ COST FUNCTION -----------------

subroutine cost(pars,dataFile,FIT)
  implicit none
  real					:: dataFile(:,:,:)
  real 					:: pars(:)
  real					:: siteError(size(dataFile,3))									
  real					:: GCCtmp(size(dataFile,1))
  real,intent(out) 		:: FIT
  integer 				:: masksize, i, nrsites, nrPar, N
  logical				:: mask(size(dataFile,1))
  real, allocatable 	:: Xi(:), Yi(:), years(:)
  real					:: a, b, scalingFactor
  real					:: Xhat, Yhat, R2, SSres, SStot, MAE, yrmin,&
							yrmax, nryr,RMSE
  
  N = size(pars)
  nrsites = size(dataFile,3)
  nrPar = N - nrsites
  
  do i=1,nrsites
	call fcover(pars,dataFile,i)

	! How many valid validation locations are there
	masksize = sum(dataFile(:,3,i))
		
	! allocate space for this subset of validation data
	allocate(Xi(masksize),Yi(masksize))
	
	! extract this data from the original dataset, create subset
	mask=dataFile(:,3,i) == 1 
	
	! dump all data in a temporary matrix
	GCCtmp = dataFile(:,4,i)	
	
	! after donohue et al.	else
	scalingFactor = (1 * dataFile(1,14,i)) / (dataFile(1,14,i) + 245.478027)
	
	! subset GCC as validation set, and fCover values (V)
	!obs = pack(GCCtmp,mask)
	Xi = pack(GCCtmp,mask)
	Yi = pack(dataFile(:,11,i),mask)
	Xi = scalingFactor * Xi
	years = pack(dataFile(:,1,i),mask)
	yrmin = minval(years)
	yrmax = maxval(years)
	
	if (yrmax == yrmin) then
		nryr = 1
	else
		nryr = yrmax - yrmin
	end if
	
	! calculate the means of observed and predicted
	Xhat = sum(Xi)/size(Xi)
	Yhat = sum(Yi)/size(Yi)
	
	! calculate sum of squares observed vs. mean / obs. vs. pred.
	SStot = sum((Xi - Xhat)**2)
	SSres = sum((Xi - Yi)**2)
	RMSE = sqrt(SSres/size(Xi))

	! calculate mean absolute error weighted by the mean of the series
	! giving less weight to very pronounced series
	MAE = (sum(abs(Xi - Yi)) / nryr ) / Xhat
	siteError(i) = MAE

	! clear formatting / memory allocation of the obs / pred values
	! this is necessary as the size of the validation data can change
	! between the sites
	deallocate(Xi,Yi)
    
  end do
  
  ! average error values across the sites should be minimized
  FIT = sum(siteError)/size(siteError)
  
end subroutine cost

end module MODEL

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

end module MODEL

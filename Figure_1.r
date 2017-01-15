#!/usr/bin/Rscript

# This is a rather ugly R scripts which generates a rather beautiful plot
# as shown in my Nature Climate Change paper: "Productivity of North American grasslands
# is increased under future climate scenarios despite rising aridity" (doi:10.1038/nclimate2942).
# You can execute the file as a script on linux or source the file on an OSX system
# to recreate Figure 1 from the study, using the Fortran code included in this example.
#
# Running the code will generate Figure 1 (Figure_1.pdf) in the output directory of the
# root of the project. All paths are relative.
#
# Before using this example and the especially the Fortran model, be aware of the
# LICENSE agreement. All modifications of this code should be accessible by me
# and the public at large as stipulated under the Affero General Public License V3.

# load libraries
require(zoo)

# compile and run the code
system('gfortran -ffree-form -ffree-line-length-200 -g func.f90 inc.f90 phenograss.f90 main.f90 -o phenograss', wait = T)
system('./phenograss ./parameters/sites.txt r',wait = T)

# do a system call to see which sites to plot
sites_data = system("sed -e '/^#.*$/d' -e '/^$/d' -e 's/\"//g' ./parameters/sites.txt",intern=TRUE)

# remove leading directory
sites = basename(sites_data)

# extract first chunk of after split == site name
sites = sapply(strsplit(sites,split="_"),"[[",1)

# read parameter values
model_par = as.matrix(read.table('./parameters/optimized_model_parameters.txt'))
nr_par = dim(model_par)[1]
nr_sites = length(sites)
slope = model_par[7,1]
intercept = model_par[8,1]

# read all data
gcc = read.table("./output/Gcc.txt")
gcc[gcc==-9999] = NA

# find areas which don't have Gcc values
# to cull the areas where there is no PhenoCam data
gcc_loc = apply(gcc,2,function(x)na.approx(x,maxgap=60,na.rm = F))

precip = read.table("./output/precip.txt")
precip[precip==0] = NA

# don't include precip data for areas without
# Gcc data
precip[is.na(gcc_loc)] = NA

fcover = read.table("./output/modelrun_01_fCover.txt")
# limit model output to the range of the Gcc data
# don't include spin-up
fcover[is.na(gcc_loc)] = NA 

# scaling factor to be calculated from MAP and the
# optimized slope parameter
MAP = read.table("./output/MAP.txt")
scaling_factor = (1 * MAP[1,]) / (MAP[1,] + slope)

# padding to accomodate he inset scatter plot
offset = 2095 

leg = c('a','b','c','d','e','f')
phenocam_site_letter = c('g','c','j','m','f','d')
site_order=c(1,2,3,4,5,6)
metadata = read.table('./data/metadata.csv',sep=',',header=TRUE)

pdf("./output/Figure_1.pdf",7,13)
par(xaxs='i',yaxs='i')
par(mar=c(0,5,1,2),oma=c(4,1,1,3),lwd=2,cex=1.2)
layout(matrix(c(1,2,3,4,5,6,7), 7,1, byrow=T))
k=0

for (i in site_order){
  k=k+1
  
  # read data from original input file
  input_data = read.table(sites_data[i],skip=12)
  
  # get modis data if available
  lat_long = scan(sites_data[i],what=character(),nmax=9)
  lat = as.numeric(lat_long[6])
  long = as.numeric(lat_long[9])
  
  # construct dates string
  input_dates = as.Date(paste(input_data[,1],input_data[,2],sep="-"),format="%Y-%j")
  
  # calculate spinput length
  spinup_length = input_dates[min(which(!is.na(gcc[,i]) == TRUE),na.rm=TRUE)]
  
  # length of the time series
  ts_length = length(input_dates)
  
  # plot container
  plot(input_dates,gcc[1:ts_length,i],
       bty='n',
       type='n',
       ylab='',
       yaxt='n',
       xaxt='n',
       ylim=c(0,1.3),
       xlim=c(as.numeric(input_dates[ts_length-offset]),as.numeric(input_dates[ts_length]))
  )
  
  # format dates axis
  locs <- tapply(X=input_dates, FUN=min, INDEX=format(input_dates, '%Y%m'))
  at = input_dates %in% locs
  at = at & format(input_dates, '%m') %in% c('01', '07') 
  
  # call a new plot
  par(new=TRUE)
  
  # plot the precipitation
  plot(input_dates,precip[1:ts_length,i],
       bty='n',
       xaxt='n',
       yaxt='n',
       xlab='',ylab='',
       type='h',
       lwd=2,
       col='deepskyblue',
       xlim=c(as.numeric(input_dates[ts_length-offset]),as.numeric(input_dates[ts_length])),
       ylim=c(0,150)
  )
  axis(4,cex.axis=1.2,tck=0.03,las=2,col='deepskyblue',col.axis='deepskyblue',at=c(0,70,140),lwd=1.5)
  
  # call a new plot to plot GCC/fcover
  par(new=TRUE)
  
  # plot gcc
  plot(input_dates,gcc[1:ts_length,i] * scaling_factor[1,i],
       bty='n',
       ylab='fCover',
       yaxt='n',
       xaxt='n',
       ylim=c(0,1),
       pch=20,
       lty=1,
       col='grey25',
       lwd=2,
       cex.lab=1.5,
       xlim=c(as.numeric(input_dates[ts_length-offset]),as.numeric(input_dates[ts_length]))
  )
  
  # add fcover
  lines(input_dates,gcc[1:ts_length,i] * scaling_factor[1,i],
        bty='n',
        yaxt='n',
        xaxt='n',
        ylim=c(0,1),
        lty=1,
        col='grey25',
        lwd=2.3,
        cex.lab=1.3,
        xlim=c(as.numeric(input_dates[ts_length-offset]),as.numeric(input_dates[ts_length])))
  
  # add fcover
  lines(input_dates,fcover[1:ts_length,i],
        ylim=c(0,1),
        lty=1,
        cex=0.2,
        lwd=2.3,
        col='red')
  
  legend('topright',legend=leg[k],bty='n',cex=2.5)
  
  # find location based upon name
  loc = which(metadata$site == sites[i])

  # site letter
  legend('top',paste(sites[i],"-",metadata$letter[loc]),
         bty='n',
         cex=1.5,
         bg='white')
  
  # add nice axis
  if(i == 6){
    axis(side=1, at=seq(as.Date("2010/1/1"), as.Date("2015/1/1"), by = "year"),labels=F,tck=0.03,cex.axis=1.2,lwd=1.5)
    axis(side=1, at=seq(as.Date("2010/6/30"), as.Date("2015/6/30"), by = "year"),labels=2010:2015,tck=F,cex.axis=1.2,lwd=1.5)
  } else {
    axis(side=1, at=seq(as.Date("2010/1/1"), as.Date("2015/1/1"), by = "year"),labels=F,tck=0.03,cex=1.5,lwd=1.5)
  }  

  axis(2,at=c(0,0.5,1),cex.axis=1.2,tck=0.03,las=2,lwd=1.5) 
  mtext("Precip. (mm)",4,3,cex=1,col='deepskyblue')
  
}

# legend
frame()
legend('center',horiz = T,legend=c('Observed fCover','Predicted fCover','Precipitation'),
       pch=c(20,NA,NA),lty=c(NA,1,1),lwd=2,col=c('grey25','red','deepskyblue'),bty='n',cex=1.5)

# inset specs
k = 0
top = seq(1/7,1,1/7) - 0.01
bottom = top - 0.09
top = rev(top[2:7])
bottom= rev(bottom[2:7])

for (j in site_order){
  df = na.omit(cbind(gcc[,j] * scaling_factor[1,j],fcover[,j]))
  
  k = k + 1
  par(fig=c(0.1,0.40,bottom[k],top[k]), new = TRUE)
  plot(df,
       xaxs='i',yaxs='i',
       xaxt='n',
       yaxt='n',
       bty='n',
       xlab='',
       ylab='',
       type='n',
       xlim=c(0,0.95),
       ylim=c(0,0.95)
  )
  
  R = cor(df)[1,2]
  
  rect(0,0,0.8,0.8,col = "white",border=0)
  points(df,
         pch=20,
         col='grey50')
  
  # labels
  lines(c(0,0.8),c(0,0.8),lty=2)
  legend(-0.15,0.9,legend=paste("R = ",round(R,2),sep=''),bty='n',cex=1.3)
  axis(1,at=c(0,0.4,0.8),tck=0.03,labels=c(0,'',0.8),cex.axis=1.2,lwd=1.5)
  axis(2,at=c(0,0.4,0.8),tck=0.03,labels=c(0,'',0.8),cex.axis=1.3,lwd=1.5)
  mtext('Observed',1,2,cex=0.9)
  mtext('Predicted',2,2,cex=0.9)
  
}
dev.off()
require(zoo)
graphics.off()

error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y, angle=90, code=3, length=length, ...)
}

# compile code
system('gfortran -ffree-form -ffree-line-length-200 -g func.f90 inc.f90 phenograss.f90 sann.f90 main.f90 -o phenograss', wait = T)

# run the code for the konza site, using site specific meteorological drivers
system('./a.out ./parameters/konza.txt r',wait = T)

# local weather
drivers = read.table("data/drivers/konza_grass.csv")
drivers_year = drivers[,1]
precip = drivers[,7]

# read konza yearly productivity measurements
df_bio = read.table('data/ancillary/konza_total_annual_biomass.csv',sep=',',header=T)
df_bio[is.na(df_bio)] = 0

# use florence soil: The Florence series consists of deep, well drained, moderately slowly permeable soils on uplands.
df_bio = df_bio[c(df_bio$SOILTYPE=="fl"),]

# extract the dates on which the measurements are made as well as the DOY and the year
dates_bio = as.Date(paste(df_bio$RECYEAR,df_bio$RECMONTH,df_bio$RECDAY,sep='-'),"%Y-%m-%d")
doy = as.numeric(format(dates_bio,"%j"))
year = as.numeric(format(dates_bio,"%Y"))

# make the variables human readable
grass = as.numeric(as.matrix(df_bio$LVGRASS))
forbs = as.numeric(as.matrix(df_bio$FORBS))
death_cur = as.numeric(as.matrix(df_bio$CUYRDEAD))
woody = as.numeric(as.matrix(df_bio$WOODY))

death_prev = as.numeric(as.matrix(df_bio$PRYRDEAD))
death_prev = as.numeric(c(death_prev[-1],NA))

# we will consider all biomass of perennial nature (death from this season or alive)
biomass = grass + forbs + death_cur
biomass = biomass * 10 # values are in g / 0.1 m2

# read in the modelled results
model= unlist(read.table('output/modelrun_01_fCover.txt'))

# read in the input data to extract the dates vector
dates =as.Date(paste(drivers[,1],drivers[,2],sep='-'),"%Y-%j")

bio_year = format(dates_bio,"%Y")
year = as.numeric(format(dates[order(dates)],"%Y"))
doy_drivers = as.numeric(format(dates[order(dates)],"%j"))
year_drivers = as.numeric(format(dates[order(dates)],"%Y"))
bio_doy = as.numeric( format(dates_bio,"%j"))

growth = model
growth_increase = c(diff(growth),NA)
growth_increase[which(growth_increase < 0) ] = NA

annual_model = by(growth,INDICES = list(year),function(x,...)sum(x[1:260],na.rm=T))

sums = by(biomass,INDICES=list(bio_year),FUN=mean,na.rm=T)
sds = by(biomass,INDICES=list(bio_year),FUN=sd,na.rm=T)

model_names = names(annual_model)
bio_names = names(sums)

pdf("output/Figure4.pdf",21,10)
par(mfrow=c(1,2),mar=c(5,5,3,1),tck=0.03,lwd=2.5,cex=2.5)
dd = cbind(sums,sds,annual_model[which(model_names %in% bio_names)])

plot(dd[,3],dd[,1],
     xlab='Annual intergral of fCover',
     col='olivedrab4',
     ylab=expression("Biomass (g m" ^-2 *" yr" ^-1*")"),
     pch=19,
     bty='n',
     xlim=c(60,120),
     ylim=c(125,700),
     yaxt='n',
     xaxt='n',
     xaxs='i',
     yaxs='i')
axis(2,at=c(125,250,500),tck=0.03,lwd=2)
axis(1,tck=0.03,lwd=2)

# omit first 10 years as spinup
fit = lm(dd[,1]~dd[,3])
print(cor(dd[,3],dd[,1]))
print(summary(fit))
abline(fit,lty=2,col='grey')
coef = round(fit$coefficients,2)

text(65,550,expression("R" ^2 * "= 0.54; p < 0.001"),pos=4)
text(65,490,paste("y = ",coef[1]," + ",coef[2],"x",sep=""),pos=4)
legend('topright',legend="a",bty='n',cex=2)

##
par(mar=c(5,1,3,5))
p <- barplot(t(dd[,1]),ylim=c(125,700),xpd=FALSE,col='wheat',
             xlab="Year",
             ylab="",
             border=NA,
             yaxt='n')
axis(2,at=c(125,250,500),tck=0.03,labels=FALSE,lwd=2)
#mtext(expression("Biomass (g m" ^-2 *" yr" ^-1*")"),2,3,adj=0.2)

error.bar(p,dd[,1],dd[,2],col='wheat3')

par(new=TRUE)
lines(p,-28.97+dd[,3]*4.25,type='b',
      xaxt='n',
      yaxt='n',
      col='olivedrab4',
      lwd=3,
      xlab='',
      ylab='',
      lty=2,
      cex=0.8,
      pch=19)
axis(4,at=c(125,250,500),labels=round(c(125,250,500)/4.25),tck=0.03,lwd=2)
mtext("Annual integral of fCover",4,3,cex=2.5)
legend('topright',legend="b",bty='n',cex=2)

graphics.off()

# Note: this relies on old libraries YMMV
# this is a simplified version of the figure
# in the final publication - for reference only

library(mapdata)
library(maptools)
library(maps)
library(raster)
library(RColorBrewer)

#---- set colours ----

red = colorRampPalette(brewer.pal(9,"YlOrRd"))(40)
blue = colorRampPalette(brewer.pal(9,"PuBu"))(20)
jet.col.arid = rev(c(rev(red),blue))

ylgn = colorRampPalette(brewer.pal(5,"YlGn"))(54)
ylorbr = colorRampPalette(brewer.pal(5,"YlOrBr"))(11)
ylgn = ylgn[c(-1:-4)]
ylorbr = ylorbr[-1]
jet.col.fcover = c(rev(ylorbr),ylgn)

#---- load data ----

fcover = raster("data/NCC_map_data/figure_2_fCover_change.tif")
aridity = raster("data/NCC_map_data/figure_2_aridity_index_change.tif")

#---- plot data ----

# margins etc
par(mfrow=c(2,1),oma=c(2,0,0,0),mar=c(4,4,4,4))

# plot
plot(fcover,box=F,axes=F,legend.only=F,zlim=c(-20,100),add=F,
     ylim=c(22,52.875),
     col=jet.col.fcover,
     horiz = T,
     lwd=1,
     cex.lab=1.5,
     legend.width=0.5, legend.shrink=0.9,
     axis.args=list(cex.axis=1.5),
     legend.args=list(
             text='Fractional cover change (%)',
             side=1,
             font=1,
             line=3,
             cex=1.5)
    )

# country outlines
map("worldHires","usa",resolution=0,add=T,lwd=1.1,xlim=c(-130,-70),ylim=c(22,52.875))
map("worldHires","canada",resolution=0,add=T,lwd=1.1,xlim=c(-130,-70),ylim=c(22,52.875))
map("worldHires","mexico",resolution=0,add=T,lwd=1.1,xlim=c(-130,-70),ylim=c(22,52.875))
legend('bottomright',legend='a',bty='n',cex=2.5)

plot(aridity,box=F,axes=F,legend.only=F,
     zlim=c(-0.1,0.2),
     add=F,
     ylim=c(22,52.875),
     col=jet.col.arid,
     horiz = T,
     lwd=1,
     cex.lab=1.5,
     legend.width=0.5, legend.shrink=0.9,
     axis.args=list(cex.axis=1.5),
     legend.args=list(text=expression(Delta*" Aridity Index"), side=1, font=1, line=3, cex=1.5))

# country outlines
map("worldHires","usa",resolution=0,add=T,lwd=1.3,xlim=c(-130,-70),ylim=c(25.125,52.875))
map("worldHires","canada",resolution=0,add=T,lwd=1.3,xlim=c(-130,-70),ylim=c(25.125,52.875))
map("worldHires","mexico",resolution=0,add=T,lwd=1.3,xlim=c(-130,-70),ylim=c(25.125,52.875))
legend('bottomright',legend='b',bty='n',cex=2.5)

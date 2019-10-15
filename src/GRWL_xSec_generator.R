# GRWL_xSec_generator.R

# by George Allen, Sept 2016
# built for AMHG Grant 

# runs through each GRWL centerline, calculates the orthogonal
# direction to the along stream direction at each vertex. 
# this runs assuming a greographic lat/lon projection.


# load packages (you must install these before you can run this script):
require(foreign)
require(geosphere)
require(rgdal)
require(shapefiles)


########################################################################
# input parameters:

# Define paths: 
wd = "E:/research/2019_10_13_QfromSpace/git/QfromSpace"

# path of directory containing GRWL shapefiles file(s): 
inDir = paste0(wd, "/in/GRWL/GRWL_vector_V01.02")
# where cross section shapfiles will be written:
outDir = paste0(wd, "/out/GRWL/")


n = 5 # N centerline vertices overwhich to calculate direction (must be odd numbers >1)
res = 30 # spatial resolution of dataset (m)
wt = c(5,5,3,1,1)/15 # weights for the weighted mean calculation
xLength = 6 # multiplier controling length to draw Xsection lines (e.g. 3x GRWL width)
xSpacing = 28 # n pixel spacing between each cross section (1 km = ~27.7 pixels)
minSegLen = 28 # n pixels (10 km = ~277 pixels)




########################################################################
# functions:

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,length(existingDF)+1)] = existingDF[seq(r,length(existingDF))]
  existingDF[r] = newrow
  return(existingDF)
}

spatialJumps <- function(tab, lag=1){
  # calculate spatial jumps over a specified nextdoor neighbor vertices:
  x = tab$lon
  y = tab$lat
  w = tab$width
  l = nrow(tab)
  
  p1x = x[-(l-lag+1)]
  p1y = y[-(l-lag+1)]
  p2x = x[-lag]
  p2y = y[-lag]
  
  p1 = cbind(p1x, p1y)
  p2 = cbind(p2x, p2y)
  
  sd = distGeo(p1, p2)
  sj = which(sd > res*2)+1
  
  return(sj)
}


########################################################################
# get list of dbf files to alter:
dbfPs = list.files(inDir, 'dbf', full.names=T)
dbfNs = list.files(inDir, 'dbf', full.names=F)
prjPs = sub('dbf', 'prj', dbfPs)

outPs = sub('.dbf', '', paste0(outDir, dbfNs))
pdfPs = paste0(outPs, ".pdf")
outPrjPs = paste0(outPs, ".prj")


########################################################################
# for each GRWL shapefile:
print(paste("N shapefies to process:", length(dbfPs)))

for (h in 1:length(dbfPs)){ # h = 190 #yellow  h = 528 # Lena # h= grep("SB20", dbfPs)
  
  ########################################################################
  # calculate cross sectional direction at each vertex:
  
  # read in GRWL shapefile dbf: 
  tab = foreign::read.dbf(dbfPs[h])
 
  # exclude non-river water bodies:
  tab = tab[tab$lakeFlag == 0,]
  
  # if N records in table is less than minimum segment length, skip to the next tile:
  if (nrow(tab) < minSegLen){
    print(paste(h, dbfNs[h], "skipped - too few records"))
    next
  }
  
  
  # get index of spatial jumps, calculated by comparing adjacent locations:
  sj = spatialJumps(tab, lag=1)
  sj = c(1, sj, nrow(tab))
 
  # remove any very short segments that throw a wrench in the process below:
  shortSegs = which(diff(sj) < minSegLen) # 277 = ~10 km
  
  if (length(shortSegs) > 0){
    rmInd = vector()
    shortSegTab = cbind(sj[shortSegs], sj[shortSegs+1])
    for (i in 1:nrow(shortSegTab)){
      rmInd = c(rmInd, shortSegTab[i,1]:shortSegTab[i,2])
    }
    tab = tab[-rmInd,]
  }
  
  # if N records in table is less than minimum segment length, skip to the next tile:
  if (nrow(tab) < minSegLen){
    print(paste(h, dbfNs[h], "skipped - too few records"))
    next
  }
  
  # define variables from updated table:
  x = tab$lon
  y = tab$lat
  w = tab$width
  l = nrow(tab)
  
  # chop start and end of vectors calculate bearing bewtween neighbors: 
  p1x = x[-c((l-n+2):l)]
  p1y = y[-c((l-n+2):l)]
  p2x = x[-c(1:(n-1))]
  p2y = y[-c(1:(n-1))]
  
  p1 = cbind(p1x, p1y)
  p2 = cbind(p2x, p2y)
  
  
  # calculate distance (in meters) between two adjacent vertices 
  # to find big jumps and remove them from this calculation:
  d = distGeo(p1, p2)
  j = which(d > res*n*2) + floor(n/2)
  
  # calculate bearing (rather than slope account for the distortion of lat lon):
  b = bearing(p1, p2)
  b = c(rep(-999, floor(n/2)), b, rep(-999, floor(n/2))) # make vector original length
  b[j] = -999
  
  
  # calculate spatial jumps again with updated table:
  sj = spatialJumps(tab, lag=1)
  
  if (length(sj) > 0){
    
    #sj = c(1, sj, nrow(tab))
    
    # insert NAs at start, end, and jump in vector:
    for (i in rev(1:length(sj))){
      x = insertRow(x, NA, sj[i])
      y = insertRow(y, NA, sj[i])
      b = insertRow(b, NA, sj[i])
      w = insertRow(w, NA, sj[i])
    }
    
    # length(which(is.na(b)))
    # length(which(is.na(w)))
    
    x = c(NA, x, NA)
    y = c(NA, y, NA)
    b = c(NA, b, NA)
    w = c(NA, w, NA)
    
    
    # get bounds of NAs:
    jNA = which(is.na(b))
    
    
    closeL = jNA[-1] - 1
    closeR = jNA[-length(jNA)] + 1
    farL = jNA[-1] - floor(n/2)
    farR = jNA[-length(jNA)] + floor(n/2)
    
  
    
    # pdfOut = paste0(outDir, "temp.pdf")
    # pdf(pdfOut, width=10, height=10)
    
    # use a linearly shrinking window to calculate bearing at ends of vectors:
    for (i in 1:(length(jNA)-1)){
      
      fL = farL[i]:closeL[i]
      rL = closeR[i]:farR[i]
      
      
      for (ii in 1:length(fL)){
        
        # calculate all points on left sides of jumps:
        L = c((fL[ii]-floor(n/2)), closeL[i])
        b[fL[ii]] = bearing(cbind(x[L[1]], y[L[1]]), cbind(x[L[2]], y[L[2]]))
  
        # handle all points on right sides of vectors:
        R = c(closeR[i], (rL[ii]+floor(n/2)))
        b[rL[ii]] = bearing(cbind(x[R[1]], y[R[1]]), cbind(x[R[2]], y[R[2]]))
        
        # plot(x[c((fL[1]-20):fL[2])], y[c((fL[1]-20):fL[2])])
        # points(x[fL[2]], y[fL[2]], col=4, pch=120)
        # points(x[fL[ii]], y[fL[ii]], col=2)
        # lines(c(x[L[1]], x[L[2]]), c(y[L[1]], y[L[2]]))
        # mtext(paste(i, ii, b[fL[ii]]))
        # 
        # plot(x[c((rL[1]):(rL[2]+20))], y[c((rL[1]):(rL[2]+20))])
        # points(x[rL[1]], y[rL[1]], col=4, pch=120)
        # points(x[rL[ii]], y[rL[ii]], col=2)
        # lines(c(x[R[1]], x[R[2]]), c(y[R[1]], y[R[2]]))
        # mtext(paste(i, ii, b[rL[ii]]))
      }
      
    }
    
    # length(which(is.na(x)))
    # length(which(is.na(b)))
    # dev.off()
    # cmd = paste('open', pdfOut)
    # system(cmd)
    
    # occationally, there are a situations where the GRWL
    # centerline is clipped such that there is only 1 segment, 
    # thus introducing NANs into the bearing calculation above.
    # fill these NANs in with a bearing of 90. 
    b[which(is.nan(b))] = 90
    b[jNA] = NA
    
    
  }
  
  ####
  # interpolate/extrapolate any missing width values:
  # this is pretty ugly code.
  
  wNA = which(w<30)
  if (length(wNA) > 0){
    jw = c(0, which(wNA[-1]-wNA[-length(wNA)]>1), length(wNA))
    
    #jwNA = is.na(w[wNA[jw]+1])
    
    win = 5 # number of measurements to use to fill in missing widths
    
    # for each block of missing width values:
    for (i in 1:(length(jw)-1)){
      
      # find out what is on either side of the block of missing widths:
      wB = wNA[(jw[i]+1):jw[i+1]]
      lwB = length(wB)
      
      L = w[wB[1]-1]
      R = w[wB[lwB]+1]
      
      # if a block of missing widths does not contain a jump on either end,
      
      if (!is.na(L) & !is.na(R)){
        # linear interpolation:
        m = (w[wB[lwB]+1]-w[wB[1]-1])/(lwB+1)
        w[wB] = m * (1:lwB) +  w[wB[1]-1]
        
        # interpolate between missing width values with a cubic spline:
        #spx = c(c(1:win), c(1:win)+lwB)
        #spy = c(w[c((wB[1]-win):(wB[1]-1))], w[c((wB[lwB]+1):(wB[lwB]+win))])
        #spf = splinefun(spx, spy)
        #w[wB] = spf(c((win+1):(win+lwB)))
        
      }else{
        
        # if one side of block contains spatial jump, use data from other 
        # end of block to take a weighted average width:
        if (is.na(L) & !is.na(R)){
          if (wB[lwB] < length(x)-win){
            mW = weighted.mean(w[c((wB[lwB]+1):(wB[lwB]+5))], wt, na.rm=T)
          }else{
            # unusual case where NAs are at the very end of the vector:
            mW = 100
          }
        }
        
        if (!is.na(L) & is.na(R)){
          if (wB[1] > win){
            mW = weighted.mean(w[c((wB[1]-5):(wB[1]-1))], wt, na.rm=T)
          }else{
            mW = 100
          }
        }
        
        w[wB] = rep(mW, lwB)
      }
      
      # if there is a jump on both sides of the NA bloack, fill with 100:
      if (is.na(L) & is.na(R)){
        print("BOTH ENDS OF MISSING WIDTH BLOCK = NA")
        w[wB] = rep(100, lwB)
      }
    }
  }
  
  #w[which(w<30)] # there should be no widhts less than 30]
  
  x = as.vector(na.omit(x))
  y = as.vector(na.omit(y))
  b = as.vector(na.omit(b))
  w = as.vector(na.omit(w))
  
  # j = which(d > res*n*2)
  
  # convert azimuth to quadrant degree coordinate system 
  # (0 degrees to the right, counter clockwise rotation) :
  q = 90-b
  q[q < 0] = q[q < 0] + 360
  
  #g = tan(q*pi/180) # calculate gradient (rise/run)
  #g1x = x - cos(q*pi/180)*1e-4
  #g1y = y - sin(q*pi/180)*1e-4
  #g2x = x + cos(q*pi/180)*1e-4
  #g2y = y + sin(q*pi/180)*1e-4
  
  # test angles to make sure this degree to gradient conversion is correct:
  #s = 1:100
  #polar.plot(s, q[s], main="Test Polar Plot", lwd=3, line.col=4)
  #for (i in s){abline(0, g[s][i])}
  
  p1 = cbind(x + sin(q*pi/180), y - cos(q*pi/180))
  p2 = cbind(x - sin(q*pi/180), y + cos(q*pi/180))
  
  # cross section length:
  # 3x GRWL width (from Yang et al. RivWidthCloud):
  xD = xLength*((w)/distGeo(p1, p2))
  
  # plus an additonal 20% if the channel is multichannel
  # uniqSegID = unique(tab$segmentID)
  # mnChan = rep(1, nrow(tab))
  # for (i in 1:length(uniqSegID)){
  #   bool = uniqSegID[i] == tab$segmentID
  #   mnChan[bool] = mean(tab$nchannels[bool], na.rm=T)
  # }
  # xD = xD*mnChan
  
  o1x = x + sin(q*pi/180)*xD
  o1y = y - cos(q*pi/180)*xD
  o2x = x - sin(q*pi/180)*xD
  o2y = y + cos(q*pi/180)*xD
  
  ####
  # recalculate jumps, this time over a single nextdoor neighbor vertices:
  p1x = x[-l]
  p1y = y[-l]
  p2x = x[-1]
  p2y = y[-1]
  
  p1 = cbind(p1x, p1y)
  p2 = cbind(p2x, p2y)
  
  sj = which(distGeo(p1, p2) > res*2)+1
  
  if (length(sj) < 0){
    for (i in rev(1:length(sj))){
      x = insertRow(x, NA, sj[i])
      y = insertRow(y, NA, sj[i])
    }
  }
  
  
  # keep only a small sample of evenly-sapced cross sections: 
  xI = seq(xSpacing/2, nrow(tab), xSpacing)
  
  
  
  
  # plot:
  # plot(x,y,type='l', lwd=0.4, main=h, asp=1)
  # points(x[sj], y[sj], cex=0.7, col=2, pch=16)
  # points(x[sj-5], y[sj-5], cex=0.7, col=4, pch=17)
  # segments(o1x[xI], o1y[xI], o2x[xI], o2y[xI], col=2, lwd=.1)
  
  ####
  # PLOT as PDF: 
  
  # pdf(pdfPs[h], width=100, height=100)
  # 
  # plot(x, y, type='l', asp=1, lwd=.1, col=1,
  #      xlab="lon", ylab='lat')
  #segments(g1x, g1y, g2x, g2y, col="light blue", lwd=.05)
  #col = rainbow(1000)
  
  #segments(o1x[xI], o1y[xI], o2x[xI], o2y[xI], col=2, lwd=.1)
  #identify(x, y)
  
  #z = which(duplicated(cbind(x, y)))
  #points(x[z], y[z], cex=10, col=4)
  #points(x[z], y[z], cex=.005, col=4)
  
  # dev.off() # close writing pdf file
  

  # WRITE OUT:
  
  # convert bearing to azimuth:
  b[b < 0] = b[b < 0] + 360
  
  # calculate orthogonal to azimuth:
  xDir = b
  xDir[xDir > 360] = xDir[xDir > 360] - 360
  
  tab$azimuth = xDir
  tab$width_m = w
  
  # update original dbf to include azimuth 
  # and extrapolated width data: 
  #write.dbf(tab, dbfPs[h])
  
  
  # Create cross section shapefiles:
  
  # create polygon shapefile:
  X = c(o1x[xI], o2x[xI])
  Y = c(o1y[xI], o2y[xI])
  ID = rep(1:length(o1x[xI]), 2)
  Name = unique(ID)
  
  dd = data.frame(ID=ID, X=X, Y=Y)
  ddTable = data.frame(ID=Name, lat_dd=tab$lat[xI], lon_dd=tab$lon[xI],
                       width_m=tab$width[xI], xDir=xDir[xI])
  ddShapefile = convert.to.shapefile(dd, ddTable, "ID", 3)
  
  
  # write shapefile:
  
  write.shapefile(ddShapefile, outPs[h], arcgis=T)
  
  # copy prj file:
  file.copy(prjPs[h], outPrjPs[h], overwrite=T)
  
  print(paste(h, dbfNs[h], "done run!"))
  
}


outDBFs = list.files(outDir, ".dbf", full.names=T)
nXsec = vector()
for (i in 1:length(outDBFs)){
  nXsec[i] = nrow(foreign::read.dbf(outDBFs[i]))
  print(paste(i, nXsec[i]))
}
sum(nXsec)
max(nXsec)

#1.9M cross sections total
#15.6k most in a single tile

# clipped to MERIT GT ORDER3 FLOWLINES:
outDBFs = list.files("E:/research/2019_04_23_global_CO2_efflux/GIS/out/GRWL/GRWL_Xsections_MERIT_gt_order3_spJoin", ".dbf", full.names=T)
nXsec = vector()
for (i in 1:length(outDBFs)){
  nXsec[i] = nrow(foreign::read.dbf(outDBFs[i]))
  #print(paste(i, nXsec[i]))
}
sum(nXsec)
max(nXsec)

#1.1M cross sections total
#most in a single tile = 9k
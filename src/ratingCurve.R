# Ryan Rigg's original code

library(sf)
library(geosphere)
library(tidyverse)
library(sp)
library(ggplot2)
library(rgeos)
library(dplyr)

if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"rgdal" %in% rownames(installed.packages())){
  install.packages("rgdal")}; require(rgdal)
if (!"shapefiles" %in% rownames(installed.packages())){
  install.packages("shapefiles")}; require(shapefiles)
if (!"RColorBrewer" %in% rownames(installed.packages())){
  install.packages("RColorBrewer")}; require(RColorBrewer)
if (!"zyp" %in% rownames(installed.packages())){
  install.packages("zyp")}; require(zyp)

setwd("~/research/2019_10_15_QfromSpace/git/QfromSpace")

###grwl cross sections w/ grades discharge data & order by ID.  
data_raw = read.dbf("./in/NH15_3x__3spc__spatialJoin.dbf")$dbf
data = data_raw[order(data_raw$ID),] 

###Landsat estimated widths from GEE 2015-present. 
data_val = read.csv("./in/output_NH15_3x_3spc_full_imgcoll.csv")


# usgs_w = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\width_val\\input\\gaugeData\\USGS\\dailyW\\07374000.csv")
# usgs_q = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\width_val\\input\\gaugeData\\USGS\\dailyQ\\07374000.csv")
gageinfo = read.csv("./in/gaugeTable.csv")

##grwl x section length
grwl_l = 3



####Filter out gage info to only NH15 area
gageinfo_lat = gageinfo[gageinfo$LAT >28.5451 & gageinfo$LAT <32.0911,]
gageinfo_lon = gageinfo_lat[gageinfo_lat$LONG > -96.0782 & gageinfo_lat$LONG < -89.7776,]
gageinfo = gageinfo_lon



#### Wow! This is elaborate. Can use geosphere package to calculate 
#### great cirlce distances using lat/lon. 


##Don't read these e&n conversions together-messes up the grwl 'data'##
e = n = rep(NA, nrow(gageinfo))
zone = ceiling(data$lon_dd/6) + 30
data_1 = cbind(gageinfo, e, n)
projString = paste("+proj=utm +zone=", zone, " +datum=WGS84 +ellps=WGS84", sep='')
eCol = which('e'==names(gageinfo))
nCol = which('n'==names(gageinfo))
for (i in 1:nrow(gageinfo)){
  gageinfo[i, c(eCol,nCol)] = project(cbind(gageinfo$LONG[i], gageinfo$LAT[i]), projString[i])
}

# convert grwl xsections from Lat Long to UTM coordinates: 
e = n = rep(NA, nrow(data))
zone = ceiling(data$lon_dd/6) + 30
data = cbind(data, e, n)
projString = paste("+proj=utm +zone=", zone, " +datum=WGS84 +ellps=WGS84", sep='')
eCol = which('e'==names(data))
nCol = which('n'==names(data))
for (i in 1:nrow(data)){
  data[i, c(eCol,nCol)] = project(cbind(data$lon_dd[i], data$lat_dd[i]), projString[i])
}

#How many xsections to consider
nGRWL = 5

#Calculate distance function
distance = function(ge, gn, xe, xn){
  calculate=sqrt(((ge-xe)^2)+((gn-xn)^2))
  return(calculate)
}

#plot and determine nearest xsections
plot(gageinfo$e, gageinfo$n, type="p")
closestDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
distanceDF=as.data.frame(array(NA, c(nrow(as.vector(gageinfo)), nGRWL)))
for(i in 1:nrow(gageinfo)){
  dist=distance(gageinfo$e[i], gageinfo$n[i], data$e, data$n)
  close=order(dist, decreasing=F)[1:nGRWL]
  closeID=data$ID[close]
  ##nearest[gageinfo$SITE_NUM[i]]=order(dist)[1:nGRWL]
  
  points(data$e[close], data$n[close], pch=1, col=i)
  
  closestDF[i,]=closeID
  distanceDF[i,]=close
}

##determine closest xsections to each gage
Site_number_xsections=cbind(gageinfo$SITE_NUM, closestDF)
Site_number_distances=cbind(gageinfo$SITE_NUM, distanceDF)

# sort out rating curve headers:
tab = data
tab = tab[, order(names(tab))]
QrecCols = grep("Q[[:digit:]]", names(tab))
QrecCols = c("Q0_1_34583", "Q10_1_3458", "Q20_1_3458", "Q30_1_3458", "Q40_1_3458", "Q50_1_3458", "Q60_1_3458", "Q70_1_3458", "Q80_1_3458", "Q90_1_3458", "Q100_1_345")
wCols_rev = grep("w[[:digit:]]", names(tab))
wCols = rev(wCols_rev)
##try to sort properly
wCols = c("w100", "w90", "w80", "w70", "w60", "w50", "w40", "w30", "w20", "w10", "w0")
wOccCols = wCols[grep("flag", names(tab)[wCols], invert=T)]
wFlagCols = wCols[grep("flag", names(tab)[wCols])]


# filter:
f = tab$strmOrder >= 0 #& tab$width_m > 100 ##############typically 3
# f = 1:nrow(tab)

pTab = tab[f,]
qTab = pTab[, QrecCols]
wTab = pTab[, wOccCols]
fTab = pTab[, wFlagCols]



# set to NA measurements with values of zero or less and measuremets taken 
# where river water is located at the end of their cross section segments:
rmBoo = wTab<=0 | qTab<=0
qTab[rmBoo] = NA
wTab[rmBoo] = NA




##############################
# DHG:

# plot DHG colored by recurrence interval:
qCols = rainbow(length(QrecCols))
qCols_trans = rainbow(length(QrecCols), alpha=0.1)

aVec = bVec = r2Vec = rep(NA, length(QrecCols))

plot(range(qTab, na.rm=T), range(wTab, na.rm=T),
     main = "DHG",
     xlab = "Q (cms)",
     ylab = "w (m)",
     type="n", log="xy"
)

for (i in 1:length(QrecCols)){
  # set up vars:
  Q = qTab[,i]
  w = wTab[,i]
  
  # remove non positive values:
  boo = Q>0 & w>0
  if (!(T %in% boo)){next}
  
  #rmBoo = is.na(Q) | is.na(w) | Q<=0 | w<=0
  
  
  
  # take least squares linear regression:
  reg = lm(log(w) ~ log(Q))
  aVec[i] = reg$coefficients[[1]]
  bVec[i] = reg$coefficients[[2]]
  r2Vec[i] = summary(reg)$r.squared
  xSeq = seq(min(Q, na.rm=T), max(Q, na.rm=T), length.out=100)
  ySeq = exp(aVec[i])*xSeq^bVec[i]
  
  # plot:
  points(Q, w, pch=16, col=qCols_trans[i], cex=0.8)
  lines(xSeq, ySeq, lwd=1.8, col=1)
  lines(xSeq, ySeq, lwd=1, lty=3, col=qCols[i])
}



legend("bottomright", 
       paste0(names(tab[QrecCols]),
              "  b=", round(bVec, 2),
              "  r2=", round(r2Vec, 2)), 
       pch=16, col=qCols
)




##############################
#### AHG:

qTabLog = log(qTab)
wTabLog = log(wTab)

# qTabLog = qTabLog[, -match(c("Q0", "Q10", "Q90", "Q100"), names(qTabLog))]
# wTabLog = wTabLog[, -match(c("w100", "w090", "w010", "w000"), names(wTabLog))]

# generate exponent and coefficient tables:
regTabLog = as.data.frame(array(NA, c(nrow(qTabLog), 3)))
names(regTabLog) = c("a", "b", "r2")

for (i in seq(nrow(qTabLog))){
  # skip cross sections without at least 3 valid measurements:
  if (length(which((!is.na(as.numeric(wTabLog[i,])))))<=3){next}
  reg = lm(as.numeric(wTabLog[i,]) ~ as.numeric(qTabLog[i,]));
  regTabLog[i,] = c(reg$coefficients, suppressWarnings(summary(reg)$r.squared))
}

# might be a faster way to do it:
# regTabLog = data.frame(t(sapply(seq(nrow(qTabLog)),
#                                 function(x){
#                                   reg = lm(as.numeric(wTabLog[x,])~as.numeric(qTabLog[x,]));
#                                   coef = c(reg$coefficients, summary(reg)$r.squared)
#                                 })))

regTab = regTabLog
regTab$a = exp(regTabLog$a)

summary(regTab)


# plot AHG exponent and coefficient by R2:
pal = colorRampPalette(c("red", "blue"))
rCol = pal(10)[as.numeric(cut(regTabLog$r2, breaks=10))]
plot(regTab$b, regTab$a, log="y", col=rCol)
legend('topright', levels(cut(regTabLog$r2, breaks = 10)), col= pal(10), pch=16)


tab$AHG_a = tab$AHG_b = tab$AHG_r2 = NA

tab$AHG_a[f] = regTab$a
tab$AHG_b[f] = regTab$b
tab$AHG_r2[f] = regTab$r2

##set wd for USGS width measurements and filter to selected gages.     ###################################make sure this is working if "cannot read file directory"
wd =setwd("E:\\research\\GRWL\\GRWL_2015_present\\width_val\\input\\gaugeData\\USGS\\dailyW\\")

##filter site numbers and add variables to be the same as file name. 
Gages_char = as.character(Site_number_xsections$`gageinfo$SITE_NUM`)
Gages_char_p = paste("/0", Gages_char, sep="")
p = paste(Gages_char_p, ".csv", sep="")
trial = paste(wd, p, sep="")

#create function for processing each USGS dataset. ##########################################################
usgs_processing = function(usgs_w){
  usgs_w$measurement_dt <- substr(usgs_w$measurement_dt, 0, 10)
  usgs_w_filter = subset(usgs_w, usgs_w$measured_rating_diff!= "Poor")# &usgs_w$chan_loc_dist <= 500)
  ##Just using the dishcharge values found in the width dataset - slightly different from discharge values.
  w = as.vector(usgs_w_filter$chan_width)
  q = as.vector(usgs_w_filter$discharge_va)
  ##convert to proper units for plotting. 
  w_m = as.numeric(w)*0.3048
  q_cms = as.numeric(q)*0.02832
  q_w = cbind(q_cms, w_m)
  return(q_w)
}

##is error function. 
is.error <- function(
  expr,
  tell=FALSE,
  force=FALSE
)
{
  expr_name <- deparse(substitute(expr))
  test <- try(expr, silent=TRUE)
  iserror <- inherits(test, "try-error")
  if(tell) if(iserror) message("Note in is.error: ", test)
  if(force) if(!iserror) stop(expr_name, " is not returning an error.", call.=FALSE)
  # output:
  iserror
}

##set wd for USGS discharge measurements ####### sometimes you have to run this wd_q twice to get the correct data. 
wd_q =setwd("E:\\research\\GRWL\\GRWL_2015_present\\width_val\\input\\gaugeData\\USGS\\dailyQ\\")
usgs_q_list = paste(wd_q, p, sep="")

#create a function to process q data for simple output. working. 
usgs_q_processing = function(usgs_q){
  q_v = as.vector(usgs_q[,4])
  q_c = as.character(usgs_q[4])
  q_n = as.numeric(q_v)
  q= q_n *0.02832
  usgs_q = cbind(usgs_q, q)
  as.character(usgs_q$datetime)
  return(usgs_q)
}

##read in validation dates to assign to validation cross section measurements. 
validation_dates = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\Landsat_dates_nh15.csv")
dts = as.character(validation_dates$list)
sep_val = unlist(strsplit(dts, "\\,"))
sep_val = as.vector(sep_val)
sep_val = noquote(sep_val)
sep_val = as.character(sep_val)
sep_val = trimws(sep_val)

##read in validation ids to help assign dates for each cross section measurement. 
validation_id = read.csv("E:\\research\\GRWL\\GRWL_2015_present\\Landsat_id_nh15.csv")
id = as.character(validation_id$list)
sep_id = unlist(strsplit(id, "\\,"))
sep_id = as.vector(sep_id)
sep_id = noquote(sep_id)
sep_id = as.character(sep_id)
sep_id = trimws(sep_id)
val_comb = cbind(sep_id, sep_val)

## work on splitting data_val system so landsat ids will match each other. 
syst = data_val$system.index
syst_c = as.character(syst)
syst_1 = substr(syst_c, 0, 21)
data_val$system.index=syst_1

##attach dates to data_val based on corresponding values in val_comb table. ##############################################this data_val for loop won't stop running but if you stop it, it will produce the correct data. 
val_dates= as.vector(nrow(data_val))
for (i in 1:nrow(data_val)){
  ind = match(data_val$system.index, sep_id)
  val_dates = as.character(sep_val[ind])
}
data_val$Date = val_dates
data_val$calc_mean = data_val$mean * data_val$width_m * grwl_l



## loop to create plots with rating curves, usgs data, and landsat estimated widths. ### still producing multiple plots of gages w/ diff xsections.


##create NA columns in data so that we can bypass them in the following loop. 
for (i in 1:nrow(data)){
  notNA = !is.na(qTab[i,])
  Q = as.numeric(qTab[i,notNA])
  w = as.numeric(wTab[i,notNA])
}

##setting outputs of loop. 
xSecq=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecw=as.data.frame(matrix(numeric(), nrow =5, ncol = 11))
xSecIDcol=grep("V", names(Site_number_xsections))
mInd = array(5, dimnames = NULL)
rangedf_1 = as.data.frame(matrix(numeric(), nrow = 1, ncol = 4))
for (i in 1:nrow(Site_number_xsections)){
  for (j in 1:(ncol(Site_number_xsections))){
    xSecIDcol[j] = Site_number_xsections[i,j]
    ##if you uncomment the next two lines it will only produce the gages with 5 lines. Only 4 plots produced. 
    #if(xSecIDcol[1]==Site_number_xsections[i,1] & xSecIDcol[2]==Site_number_xsections[i,2] & xSecIDcol[3]==Site_number_xsections[i,3] ##extra step to plot only proper curves. 
    #&xSecIDcol[4]==Site_number_xsections[i,4]&xSecIDcol[5]==Site_number_xsections[i,5]&xSecIDcol[6]==Site_number_xsections[i,6]){
    xSecID=xSecIDcol[2:6]#[xSecIDcol==Site_number_xsections[i,2:6]]
    mInd = match(xSecID, data$ID)
    xSecw=wTab[mInd,notNA]
    xSecq=qTab[mInd,notNA]
    w_max=max(range(xSecw))
    w_min=min(range(xSecw))
    q_max=max(range(xSecq))
    q_min=min(range(xSecq))
    data_val_subset = subset(data_val, xSecID[1]==data_val$ID | xSecID[2]==data_val$ID| 
                               xSecID[3]==data_val$ID| xSecID[4]==data_val$ID| xSecID[5]==data_val$ID)
    data_val_sub_agg = try(aggregate(data_val_subset$calc_mean, by=list(data_val_subset$Date), FUN=mean))
    #if(!is.error(data_val_sub_agg)){next}else{next}
    rangedf_1 = cbind(q_max, q_min, w_max, w_min)
    
    ##if range is not na, then plot. 
    if(!is.na(rangedf_1)){
      plot(c(rangedf_1[,1], rangedf_1[,2]), c(rangedf_1[,3], rangedf_1[,4]),
           main=paste("USGS Gage:", Site_number_xsections$`gageinfo$SITE_NUM`[i]),
           log="", type="p", xlab="Q (cms)",
           ylab="W (m)",)
      usgs = try(usgs_processing(read.csv(trial[i])))
      
      ##if above line works then add in the usgs width/discharge information. if the usgs data is missing then skip this step. 
      if(!is.error(usgs)){
        points(usgs[,1], usgs[,2])
        
        ###get percentile USGS discharge and width values and plot them. 
        u=as.data.frame(usgs)
        u_w_q = try(quantile(u$w_m, probs = seq(.2, 1, .1), na.rm = TRUE))
        u_q_q = try(quantile(u$q_cms, probs = seq(.2, 1, .1), na.rm = TRUE))
        try(lines(u_q_q, u_w_q, type = "l", lty = 2, col = "red"))
        usgs_q = try(usgs_q_processing(read.csv(usgs_q_list[i], stringsAsFactors = FALSE)))
        
        ##if above line works then add in the usgs discharge information and filter to match with GEE estimated landsat widths for same days. 
        if(!is.error(data_val_sub_agg) & !is.error(usgs_q)){
          usgs_q_ind=which(usgs_q$datetime %in% data_val_sub_agg$Group.1)
          usgs_q_subset= usgs_q[usgs_q_ind,]
          usgs_q_subset[order(usgs_q_subset$datetime),]
          data_val_sub_agg_ind=which(data_val_sub_agg$Group.1 %in% usgs_q_subset$datetime)
          data_val_sub_agg_1 = data_val_sub_agg[data_val_sub_agg_ind,]
          data_val_sub_agg_1[order(data_val_sub_agg_1$Group.1),]
          points(usgs_q_subset$q, data_val_sub_agg_1$x, pch=17, col="blue")
          l=as.data.frame(cbind(usgs_q_subset$q, data_val_sub_agg_1$x))
        } else{next}
      } else{next}
      
      ##add in lines of rating curves to each plot. 
      for (k in 1:nrow(xSecq)){     
        print(xSecq)
        try(lines(as.numeric(xSecq[k,]), as.numeric(xSecw[k,]), type="l"))
        #text(1, 2, xSecID[1:5])
        text(as.numeric(xSecq[k,10]), as.numeric(xSecw[k,10]), xSecID[k])
        try(print(cor.test(u_q_q, as.numeric(xSecq[k,10]))))
        legend("bottomright", 
               legend = c("Landsat widths", "USGS widths"), 
               col = c("blue", 
                       "black"), 
               pch = c(17,01), 
               bty = "n", 
               text.col = "black", 
               horiz = F)# , 
        #inset = c(0.9, 0.25))
      }} else{next}
    # }else{} ##extra step to only keep proper curves. 
  }} 


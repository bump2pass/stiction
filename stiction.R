#libraries used
library(plyr)
library(quantmod)
library(ggplot2)
library(reshape2)
library(easyGgplot2)

#This script takes a wyko file and determines how many peaks here are.
#the file needs to be the ASCII format, simple header, raw data (usually 3MB size)
#the processessing is as follows
#read the data put in array
#do a gross level. do this because in the ss holders, the die are often wildly out of level
#do a find level, do this by finding the glass layer, fitting a model to that, and level to the glass level
#run a peak find on the leveled data by looking athe kernel density of the histogram of heights
#count the number of peaks. if it is not 2 or 3, then stiction declared
#we also check that the total thickness if 28.5-31.5 um

#tuning parameters
#adjust will change how smooth the kernel density estimate is
#is_peak is the percent of the max peak height that we call a peak in the KDE
#npnt is the number of points we randomly select to level the plot initially
adjust <- 1; is_peak <- 0.01; wid <- 1;npnt <- 1000


#if no pmass this value is 2, if pmass it is 3
#it doesn't change the peak finding routine
#just the logic on if there is stiction or not
npeak_expct = 3

#device layer thickness
dev_lay <-  30
dev_lay_tol <- 1.5

#if there are too many na's in the data we need to give a dunno the answer to  'is stuck?'
na_pct <- 0.2 #if more than this say 'don't know'

#if the second peak in the height histogram is within 3.5um of the base then it is pads exclude
#this paramter sets that window
pad_h <- 3.5 #max height for a peak to be considered a pad
pad_l <- 1 #min height for peak to be considered a pad

############### Analysis Functions ########################################
wyko2df <- function(fn){
  #fn is file name
#function to read in wyko data, extract the wavelength and mult parameters
#make a matrix
#note that the data file should be about 3Mb
#set to integer and short header format
vec <- scan(file = fn, n=50 , sep = ",", what = character(), quiet = TRUE)#scan in the header as a vector
mult <- as.numeric(vec[grep('Mult', vec)+3])#search for the Mult reading, a fudge factor in the data
wavelength <- as.numeric(vec[grep('Wavelength', vec)+3]) #get the wavelength reading
df <- read.table(file = fn, header = F, sep = ",", skip = 13, nrow = 736, blank.lines.skip = F, 
                 colClasses = "integer", na.strings = c("Bad", "BAD", "bad"))#read the data
ar <- as.matrix(df)#convet to matrix
dimnames(ar) <- NULL
ar <- ar*wavelength/mult/1000#wavelength is nm/wave so we conver to um/wave by div 1000
ar
}

gross_level <- function(ar, npnt = 1000){
  #function to take tilt out of wyko data
  #sample the matrix at n points
  #compute the slope at that point
  #tilt is teh median
  #npnt is the number of points we sample
  #notice that we don't sample the last row/col
  xpnt <- sample(1:(dim(ar)[1]-1), npnt, replace = TRUE)
  ypnt <- sample(1:(dim(ar)[2]-1), npnt, replace = TRUE)
  xslp <- function(x,y){
    #x,y are indicies in ar
    #the slope is reference to index units
    ar[x+1,y]- ar[x,y]
  }
  yslp <- function(x,y){
    ar[x,y+1]-ar[x,y]
  }
  xslps <- mapply(xslp, xpnt, ypnt)#check doc on mapply, applies xslp over xpnt/ypnt
  yslps <- mapply(yslp, xpnt, ypnt)
  xtilt <- median(xslps, na.rm = TRUE)
  ytilt <- median(yslps, na.rm = TRUE)
  ar_flat <- ar - row(ar)*xtilt - col(ar)*ytilt
  ar_flat #return the flattened data
}

fine_level <- function(ar, adjust, is_peak, wid){
  #do a peak find on the gross leveled data
  #pull out the lowest data, which should be the glass/circuit layer
  #do a linear fit on that
  #use those coeficients to level the data
  hgt  <- pk_fnd(ar, adjust = adjust, is_peak = is_peak)
  df<-melt(ar, varnames = c('xpos', 'ypos'), value.name = 'height')
  df <- df[!is.na(df$height),]
  #select the data around the lowest peak of the histogram of data values
  df0 <- df[(df$height>hgt[1]-wid)&(df$height<hgt[1]+wid),]
  #compute the linear model
  lm1 <- lm(height ~ xpos + ypos, data = df0)
  xtilt <- coef(lm1)[2]
  ytilt <- coef(lm1)[3]
  ar_flat <- ar - row(ar)*xtilt - col(ar)*ytilt
}

level <- function(ar, npnt, adjust, is_peak, wid){
  #compose leveling
  ar <- gross_level(ar, npnt = npnt)
  ar <- fine_level(ar, adjust = adjust, is_peak = is_peak, wid = wid)
  ar
}

pk_fnd <- function(ar, adjust, is_peak){
  #find peaks in matrix ar
  #adjust is used to adjust the bw in density kernel calculation. smaller gives more peaks.
  dar <- density(ar, na.rm = TRUE, adjust = adjust)#compute a smooth density of the histogram
  pk <- findPeaks(dar[['y']])#find the peaks in dar. 
  pkv <- dar[['y']][pk]#pull out peak values
  pks <- pk[which (pkv>(max(pkv)*is_peak))]#use only the ones that are higher than some % of the max peak 
  heights <- dar[['x']][pks]#should just be height of wing, pmass, and base, unless stiction
  heights #return the peak heights
}


################### Testing on Data #####################################
#run the same analysis on all the data
dir1 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\canonical stiction data"
dir2 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\157"
dir3 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\ss holder data no focus no tilt"
#in this data set I picked up the pads for  data set 25. Data set 15 was garbage.
dir4 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\test" 
dir5 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\test2"
#this data set I scanned in 60 parts from a lot. leveled parts by hand
#then scanned as normal (-100-200 um scan)
#the last part in the lot (60) was stuck
#the algorithm correctly id'ed this part
dir6 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\L18A192A"
#this very odd lot had high and low g parts (pmass/nopass in it)
dir7 <- "C:\\Users\\rob.macdonald\\Documents\\R.Projects\\Analyze Wing Wyko Data\\wyko data\\L18A183E"
lf1 <- list.files(dir1, recursive = T, full.names = T)
lf2 <- list.files(dir2, recursive = T, full.names = T)
lf3 <- list.files(dir3, recursive = T, full.names = T)
lf4 <- list.files(dir4, recursive = T, full.names = T)
lf5 <- list.files(dir5, recursive = T, full.names = T)
lf6 <- list.files(dir6, recursive = T, full.names = T)
lf7 <- list.files(dir7, recursive = T, full.names = T)
lf <- c(lf7)#all the files we want to analayze


stiction <- function(fn, adjust, is_peak, wid, npnt, npeak_expct, dev_lay, 
                     dev_lay_tol, na_pct, pad_h, pad_l){
  #fn is a file name we want to check for stiction
  #output is a list with various data
  #stuck if dev layer out of tolerance or wrong number of peaks
  #also not stuck if too many na's in the data 
  #we also need to exclude the 'pads'
  #this will show up as points about 2-3um from the base
  #read data, gross level, fine leve, find peaks, compare to expected # peaks
  ar <- wyko2df(fn)#fn is file name
  #does ar have more na than allowed?
  ar_na <- na_pct < (sum(is.na(ar))/length(ar))
  ar <- level(ar, npnt, adjust, is_peak, wid) #fine level the data
  pks <- pk_fnd(ar, adjust, is_peak) #the peak heights
  n <- length(pks)#the number of peaks
  #check to see if the second peak is a pad
  pads <- pad_h > (pks[2]-pks[1]) & (pks[2]-pks[1])> pad_l
  #if so then add 1 to npeak_expect
  npeak_expct <- npeak_expct + pads
  rnp <- n == npeak_expct #is n = number of peaks we expect? rnp = right number peaks
  TTV <- max(pks) - min(pks) #Total Thickness Variation
  dev <- (TTV < dev_lay + dev_lay_tol) & (TTV > dev_lay - dev_lay_tol) #is TTV about what we expect?
  free <- rnp & dev & !ar_na #return T/F on question, if T, with is Free if False, stuck
  list(ar, pks, rnp, dev, pads, free, ar_na)#result is a list with many of the things we calculated
}

#compute for each file
stuck_wing <- llply(lf, stiction, adjust = adjust, is_peak= is_peak, wid = wid, npnt = npnt, 
                    npeak_expct = npeak_expct, dev_lay = dev_lay, dev_lay_tol = dev_lay_tol,
                    pad_h = pad_h, pad_l = pad_l, na_pct = na_pct)
#and the answer is...
which(!unlist(lapply(stuck_wing, function(ls) ls[[6]])))#print out the files for which we don't have a free wing


################make plots of the results######################################
par(mfrow=c(4,2))
#unpack main result 'free'
free_wing <- unlist(lapply(stuck_wing, function(ls) ls[[6]]))
plot(as.numeric(free_wing))
#check number peaks found
npks <- unlist(lapply(stuck_wing, function(ls) length(ls[[2]]) ))
plot(npks)
#check TTVs found
TTV <- unlist(lapply(stuck_wing, function(ls) max(ls[[2]])-min(ls[[2]]) ))
plot(TTV)
#plot rpn, right number peaks
right_num_peaks <- unlist(lapply(stuck_wing, function(ls) ls[[3]]))
plot(right_num_peaks)
#plot dev, does part have correct device layer thickness
dev <- unlist(lapply(stuck_wing, function(ls) ls[[4]]))
plot(dev)
#plot if we found pads or not
pad_found <-  unlist(lapply(stuck_wing, function(ls) ls[[5]]))
plot(pad_found)
#plot which file had to many na's
to_much_na <- unlist(lapply(stuck_wing, function(ls) ls[[7]]))
plot(to_much_na)
par(mfrow=c(1,1))
which(!free_wing)#print out the files for which we don't have a free wing
which(to_much_na)
which(pad_found)
which(!dev)
############# Useful Function fo checking an wyko data set by hand###############################
hist_raw_vs_level <- function(fn, fun, ...){
  #fn is a file we are looking at which contains a wyko data set
  #plot the effects of leveling ar by fun, 
  #... are all the prams we need to pass to fun
  ar <- wyko2df(fn)
  flt <- fun(fn,...)[[1]]#take out the tilt
  par(mfrow=c(2,2))#make plots
  hist(ar, breaks = 500)
  plot(density(ar, na.rm = T))
  #contour(ar)#plot it
  hist(flt, breaks = 500)
  plot(density(flt, na.rm = T))
  #contour(flt0)
  par(mfrow=c(1,1))
}  

contour_data <- function(fn , fun, ...){
  #for a particular wyko file
  #apply the stiction function fun, with parameters passed to fun ...
  #around each peak found in the leveled data plot a contour plot
  x <- fun(fn, ...) #compute the analysis using fun
  flt <- x[[1]]#pick out the first value of fun, flattend data
  df <- melt(flt, varnames = c('xpos', 'ypos'), value.name = 'height')
  df <- df[!is.na(df$height),]
  #find peaks of flt
  hgt <- x[[2]]
  # plot all the filled contours around each peak +/- 1um
  pmaker <- function(hgt, wid){
    #make a plot for each element of hgt
    df1 <- df[(df$height>hgt-wid)&(df$height<hgt+wid),]
    p <- ggplot(df1, aes(xpos, ypos)) + geom_raster(aes(fill = height))
  }
  #for each peak in hgt, apply pmaker, return a
  p <- lapply(hgt, pmaker, wid )
  
}


#For the files which don't have a free wing, we can study them with these plotting tools
which(!free_wing)
fn <- lf[[16]]#a particular file name , we will check this files behavior
#
hist_raw_vs_level(fn, stiction , adjust = adjust, is_peak= is_peak, wid = wid, npnt = npnt, 
                  npeak_expct = npeak_expct, dev_lay = dev_lay, dev_lay_tol = dev_lay_tol,
                  na_pct = na_pct, pad_h = pad_h, pad_l = pad_l)
p <- contour_data(fn, stiction , adjust = adjust, is_peak= is_peak, wid = wid, npnt = npnt, 
                  npeak_expct = npeak_expct, dev_lay = dev_lay, dev_lay_tol = dev_lay_tol,
                  na_pct = na_pct, pad_h = pad_h, pad_l = pad_l)
ggplot2.multiplot(plotlist = p, cols = 2)



############# Junk Code#############################
#this code is for fun
#calculate pmass depth, wing tip height
#do this for a particular part


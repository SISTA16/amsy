##---------------------------------------------------------------------------------------------
## AMSY
## Original code written by Rainer Froese in January - April 2019
## Additons by Henning Winker:
## 1. Implemented process error on process equation
## 2. Process error is implemented as sigme.R = c(0.05,0.07,0.1,0.15) for Very low, Low, Medium, High.
## 3. Schaefer function now greps sigma.R as defined 
## 4. Observation error implemented as CV (log.sd of CPUE)
## 5. Implemented Kobe prototype with terminal F/Fmsy taken as mean of previous 3 yrs.
## 6. Added save.plot option from CMSY
## 7. Added automatic package installer
## Additions by Gianpaolo Coro:
## 1. improved estimate of prior kq  
## 2. retrospective analysis
## Addition by RF: MVN
##---------------------------------------------------------------------------------------------
# Install required packages if not available
list.of.packages <- c("gplots", "coda","mvtnorm","crayon") #><> mvt not mtv
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(coda) 
library("gplots")
library(mvtnorm) # used for Kobe plot, ignore version error
library(crayon) # to display bold and italics in console

#-----------------------------------------
# Some general settings ----
#-----------------------------------------
# set.seed(999) # use for comparing results between runs
rm(list=ls(all=FALSE)) # clear previous variables etc
options(digits=3) # displays all numbers with three significant digits as default
graphics.off() # close graphics windows from previous sessions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set working directory to source file location

#-----------------------------------------
# Required settings, File names 
#-----------------------------------------
 id_file     <- "EU_Stocks_ID_8.csv" #   name of file containing stock-specific info and settings for the analysis
 outfile     <- paste("Out_",format(Sys.Date(),format="%B%d%Y_"),id_file,sep="") # default name for output file

 #----------------------------------------
# Select stock to be analyzed ----
#----------------------------------------
stocks      <-NA
# If the input files contain more than one stock, specify below the stock to be analyzed
# If the line below is commented out (#), all stocks in the input file will be analyzed
# stocks <- "DIPL.ANN"   #"HL_VL"   #c("anb-78ab")  # c("HH_VL","HL_VL","HLH_VL","LH_VL","LHL_VL","LL_VL") # "anb-78ab"  #c("HL_H","HL_M","HL_L","HL_VL")  #"HH_L"  #"LHL_L"  #"Micr_pou_AD" #"Myxine glutinosa"   #"Eut_gurn_Balt"  # "rjc.27.3a47d"  #"PNSK"   # "Ille_coi_AD" # "WSTM"

 # Read data
 cinfo        <- read.csv(id_file, header=T, dec=".", stringsAsFactors = FALSE)
 cat("File", id_file, "read successfully","\n")
 
#-----------------------------------------
# General settings for the analysis ----
#-----------------------------------------
smooth.cpue  <- T # set to TRUE to apply ksmooth with minimum bandwidth of 3, increased at low r
filter       <- TRUE # set to TRUE for Monte Carlo filtering; if FALSE, incease max.viable to 20000
cor.log.rk   <- -0.607 #-0.871 # empirical value of log r-k correlation in 140 stocks analyzed with BSM
sigma.r      <- c(0.05,0.07,0.1,0.15) # very low, low, medium, high # overall process error for productivity or r
sigma.cpue   <- 0.3 # observation error for cpue 
n.p          <- 50000 # number of r-kq pairs to be analysed; will be doubled if too few viable pairs are found
n.trial      <- 30 # times each year is calculated with new random error terms for r and cpue
min.viable   <- 20 # minimum number of viable r-kq pairs to be accepted for analysis
max.viable   <- 5000 # maximum number of viable r-kq pairs to reduce processing time; set to 20000 if filter==FALSE
creep.graph  <- F # plot graph for effort creep correction, if used
do.plots     <- T # retrospective analysis does not work if FALSE
write.output <- T # set to TRUE if table with results in output file is wanted 
kobe.plot    <- F # HW set to TRUE so produce additional kobe status plot 
save.plots   <- T # set to TRUE to save graphs to JPEG files
close.plots  <- T # set to TRUE to close on-screen plots, to avoid "too many open devices" error in batch-processing;
retros       <- F # retrospective analysis, requires do.plots <- TRUE

#----------------------------------------------
#  FUNCTIONS ----
#----------------------------------------------
# Monte Carlo filtering with Schaefer Function ----
#----------------------------------------------

SchaeferCPUE<-function(yr, cpue, ri, kqi, sigR, filter){
# create matrix for results
    mdat  <- matrix(ncol = (2*nyr+1))
    colnames(mdat) <- c("rv","kqv",paste("c",yr[1:(nyr-1)],sep=""),paste("b",yr[1:nyr],sep=""))
    
    for(i in 1:length(ri)) { # for all submitted r-kq pairs
      
      for(trial in 1:n.trial) { # rerun every r-kq pair several times because error terms per year are random
                                # max one succesful run across all years is returned per trial
        cqt         <- vector()
        cpuet       <- vector()
        FFmsy       <- vector()
        break.flag  <- FALSE
        
        for(t in 1:(nyr-1))  {  # for all years except the last one, for which catch cannot be calculated
          # assign random error terms to surplus production and to cpue  
          err      <- exp(rnorm(1,0,sigR)) # set annual error for productivity  
          cpuet[t] <- cpue[t]*exp(rnorm(1,0,sigma.cpue)) # assign error to cpue
          if(cpuet[t]<=0) {cpuet[t] <- 0.01*kqi[i]} # make sure cpuet is not zero or negative
          
          # calculate catch
          if(cpuet[t]/kqi[i] >= 0.25) {
            cqt[t] <- (cpuet[t] + cpuet[t] * ri[i] * (1-cpuet[t]/kqi[i]))*err - cpue[t+1]  } else {
            cqt[t] <- (cpuet[t] + cpuet[t] * ri[i] * (1-cpuet[t]/kqi[i])*(4*cpuet[t]/kqi[i]))*err - cpue[t+1] } # reduce r linearly below 0.25 kq
 
          # use moving average to smooth hectic catch predictions
          if(t == 2) {cqt[t] <- mean(c(cqt[t-1],cqt[t])) }
          if(t > 2) {cqt[t] <- mean(c(cqt[t-2],cqt[t-1],cqt[t])) }
        
          # calculate MSYq and F/Fmsy
          MSYq     <- ri[i]*kqi[i]/4 
          FFmsy[t] <- 2*cqt[t]/(ri[i]*cpue[t])
      
          if(filter==TRUE) {
          ## Test compatibility of r-kq pairs with general prior popdyn knowledge
          ## If one test fails, break the loop and go to the next trial  
          # (1) Exclude r-kq pair if catch is negative (cqt[t] < 0) 
            mult.kqi    <- ifelse(res=="Very low",-0.06,ifelse(res=="Low",-0.02,0)) #relax rule for Very low and Low resilience
            if(cqt[t] < mult.kqi*kqi[i]) {break.flag<-TRUE;break}
        
          # (2) Exclude r-kq pair if catch exceeds biomass (cqt[t] > cpue[t])
            # if lowest cpue is close to zero, skip this test
            if(min(cpue.raw) > 0.1*max.cpue) {
            # in highly productive species, catch may exceed average annual biomass
              mult.cpue  <- ifelse(res=="High",1.4,ifelse(res=="Medium",1,ifelse(res=="Low",0.5,0.25))) 
              if(cqt[t] > (mult.cpue*cpuet[t])) {break.flag<-TRUE;break} }
             
          # (3) Exclude r-kq pair if catch exceeds MSY   
            # some overshooting of MSY is possible
             mult.msy <- ifelse(res=="Very low",10,ifelse(res=="Low",5,ifelse(res=="Medium",3,2))) 
             if(cqt[t] > mult.msy*MSYq)        {break.flag<-TRUE;break} 

          # (4) Exclude r-k pairs if F/Fmsy is highly unrealistic (negative or much too high)
             FFlow  <- ifelse(res=="Very low",-25,-3) 
             FFhi   <- ifelse(res=="Very low",12,5)   
             if(t > 1 && (FFmsy[t-1] < FFlow || FFmsy[t-1] > FFhi)) {break.flag<-TRUE;break}
          
          # (5) if relative cpue in the year of the B/k prior is outside of the prior range, discard trial
          #  relax rule if lower B/k prior range is <= 0.01    
            if(prior.Bk[1] <= 0.01) { prior.Bk[1] <- 0.0001}
             if(t==Bk.yr.i && (cpuet[Bk.yr.i]/kqi[i] < prior.Bk[1] || cpuet[Bk.yr.i]/kqi[i] > prior.Bk[2])) {break.flag<-TRUE;break }
            
          } # end of condition for filtering
          
         } # end of t-loop through years
        
        # if t-loop was broken and flag==TRUE do not test further, do not plot points, do not store results 
        if(break.flag==TRUE) { next }
        
        # assign error to last cpue and repeat filter (7) if applicable to last year
        cpuet[nyr] <- cpue[nyr]*exp(rnorm(1,0,sigma.cpue)) 
        if(cpuet[nyr]<=0) {cpuet[nyr] <- 0.01*kqi[i]} # make sure cpuet is not zero or negative
        if(filter==TRUE && Bk.yr.i==nyr && (cpuet[nyr]/kqi[i] < prior.Bk[1] || cpuet[nyr]/kqi[i] > prior.Bk[2])) { next }   
        
        # If all tests are passed, add viable r-kq pair and predicted catch to matrix
          mdat     <- rbind(mdat,c(ri[i],kqi[i],cqt[1:(nyr-1)],cpuet[1:(nyr)]))
          
          # plot viable r-kq pairs
          if(do.plots==T) {
          points(x=ri[i],y=kqi[i], col="grey20", pch=".", cex=2.5) }
        
      } # end of trial-loop for trials per r-kq pair
      if(length(mdat[,1])>max.viable) { break} # end searching for viable pairs if n > max.viable
 } # end of i-loop through r-kq pairs
      
 mdat <- na.omit(mdat)
 return(mdat)

} # end of SchaeferCPUE function

#-------------------------------------------------------------
# Function to create multivariate-normal distribution for r-k 
#-------------------------------------------------------------
mvn   <- function(n,mean.log.r,sd.log.r,mean.log.kq,sd.log.kq) {
  cov.log.rk <- cor.log.rk*sd.log.r*sd.log.kq # covariance with empirical correlation and prior variances  covar.log.rk = matrix(NA, ncol=2,nrow=2)   # contract covariance matrix
  covar.log.rk      <- matrix(NA, ncol=2,nrow=2) # covariance matrix
  covar.log.rk[1,1] <- sd.log.r^2                # position [1,1] is variance of log.r
  covar.log.rk[2,2] <- sd.log.kq^2               # position [2,2] is variance of log.k
  covar.log.rk[1,2] = covar.log.rk[2,1] = cov.log.rk     # positions [1,2] and [2,1] are correlations
  mu.log.rk  <- (c(mean.log.r,mean.log.kq))      # vector of log.means
  mvn.log.rk <- rmvnorm(n,mean=mu.log.rk,sigma=covar.log.rk,method="svd") 
  return(mvn.log.rk)
}

#---------------------------------------------
# END OF FUNCTIONS
#---------------------------------------------

#--------------------------------------------
# Create table for output to csv file
#--------------------------------------------
if(write.output==T && substr(id_file,1,3)=="Sim"){ # output for simulated data
  outheaders = data.frame("Stock","r.true", "r.est","r.lcl","r.ucl",
                          "kq.true","kq.est","kq.lcl","kq.ucl",
                          "MSYq.true","MSYq.est","MSYq.lcl","MSYq.ucl",
                          "FFmsy.true","FFmsy.est","FFmsy.lcl","FFmsy.ucl",
                          "BBmsy.true","BBmsy.est","BBmsy.lcl","BBmsy.ucl")
  write.table(outheaders,file=outfile, append = T, sep=",",row.names=F,col.names=F) }


if(write.output==TRUE && is.null(cinfo$MSY.BSM)==F) { # assuming all BSM fields are available
  outheaders = data.frame("Stock","r.BSM","lcl","ucl","r.est","lcl","ucl","k.BSM","lcl","ucl","k.est","lcl","ucl",
                          "BBmsy.BSM","lcl","ucl","BBmsy.est","lcl","ucl",
                          "FFmsy.BSM","lcl","ucl",
                          "FFmsy.est","lcl","ucl") 
  write.table(outheaders,file=outfile, append = T, sep=",",row.names=F,col.names=F)
  }

if(write.output==T && substr(id_file,1,3)!="Sim" && is.null(cinfo$MSY.BSM)==T){
  outheaders = data.frame("Stock","Fmsy.est","Fmsy.lcl","Fmsy.ucl",
                          "FFmsy.est","FFmsy.lcl","FFmsy.ucl",
                          "BBmsy.est","BBmsy.lcl","BBmsy.ucl")
  write.table(outheaders,file=outfile, append = T, sep=",",row.names=F,col.names=F) }

#-----------------------------------------
# Start output to screen
#-----------------------------------------
cat("------------------------------------------------------------\n")
cat("AMSY Analysis,", date(),"\n")

#---------------------------------
# Analyze stock(s)
#---------------------------------
if(is.na(stocks[1])==TRUE){
     stocks         <- as.character(cinfo$Stock) # Analyze stocks in sequence of ID file
    # stocks          <- cinfo$Stock[cinfo$Stock >= "fle-2425"] # Analyze stocks in sequence of ID file
    # stocks         <- sort(as.character(cinfo$Stock)) # Analyze stocks in alphabetic order
}
# analyze one stock after the other
for(stock in stocks) {
  
  #retrospective analysis
  if(retros==T && (cinfo$EndYear[cinfo$Stock==stock]-cinfo$Bk.yr[cinfo$Stock==stock])<3) {
    retros.nyears<-0 #retrospective analysis
    cat("Warning: Retrospective analysis not meaningful and omitted if B/k prior is in the final year(s)\n") } else {
    retros.nyears<-ifelse(retros==T,3,0) #retrospective analysis
  }
  
	FFmsy.retrospective<-list() #retrospective analysis
	BBmsy.retrospective<-list() #retrospective analysis
	years.retrospective<-list() #retrospective analysis
 
 for (retrosp.step in 0:retros.nyears){ #retrospective analysis
 
  cat("------------------------------------------------------------\n")
  cat("Stock ",bold(stock),", ", bold(italic(as.character(cinfo$ScientificName[cinfo$Stock==stock]))),", ",
      as.character(cinfo$EnglishName[cinfo$Stock==stock]),sep="","\n")
  # read data for stock
  cpue_file    <- cinfo$CPUE_File[cinfo$Stock==stock]
  cdat         <- read.csv(cpue_file, header=T, dec=".", stringsAsFactors = FALSE)
  # assign data from cinfo to vectors
  n            <- n.p
  res          <- as.character(cinfo$Resilience[cinfo$Stock==stock])
  res.i        <- which(c("Very low","Low","Medium","High")%in%res) # determines process error strength
  if(length(res.i)==0) {stop("Spelling error in resilience in ID file\n")}
  start.yr     <- as.numeric(cinfo$StartYear[cinfo$Stock==stock])
  end.yr       <- as.numeric(cinfo$EndYear[cinfo$Stock==stock])
  end.yr 	     <- end.yr-retrosp.step #retrospective analysis
  r.low        <- as.numeric(cinfo$r.low[cinfo$Stock==stock])
  r.hi         <- as.numeric(cinfo$r.hi[cinfo$Stock==stock])
  user.log.r   <- ifelse(is.na(r.low)==F & is.na(r.hi)==F,TRUE,FALSE)     
  Bk.yr        <- as.numeric(cinfo$Bk.yr[cinfo$Stock==stock])
  Bk.pr        <- as.character(cinfo$Bk.pr[cinfo$Stock==stock])
  Bk.pr.low    <- as.numeric(cinfo$Bk.pr.low[cinfo$Stock==stock])
  Bk.pr.hi     <- as.numeric(cinfo$Bk.pr.hi[cinfo$Stock==stock])
  e.creep      <- as.numeric(cinfo$e.creep[cinfo$Stock==stock])
  comment      <- as.character(cinfo$Comment[cinfo$Stock==stock])
  Fmsy.ass     <- as.numeric(cinfo$Fmsy.ass[cinfo$Stock==stock])
  Bmsy.ass     <- as.numeric(cinfo$Bmsy.ass[cinfo$Stock==stock])
  source       <- as.character(cinfo$Source[cinfo$Stock==stock])
  

  # check for common errors
  if(length(r.low)==0){
    cat("ERROR: Could not find the stock in the ID input file - check that the stock names match in ID and CPUE files and that commas are used (not semi-colon)")
    return (NA) }
  if(length(cdat$Year[cdat$Stock==stock])==0){
    cat("ERROR: Could not find the stock in the CPUE file - check that the stock names match in ID and CPUE files and that commas are used (not semi-colon)")
    return (NA) }
  if(start.yr < cdat$Year[cdat$Stock==stock][1]){
    cat("ERROR: start year in ID file before first year in CPUE file\n")
    return (NA)}
  if(!(Bk.pr %in% c("Near unexploited","More than half","About half","Small","Very small"))){
    cat("ERROR: Prior for stock size not in: Near unexploited, More than half, About half, Small, Very small\n")
    return (NA)}
  if(Bk.yr < start.yr || Bk.yr > end.yr){
    cat("ERROR: Year for B/k prior outside range of years\n")
    return (NA)}

  #----------------------------------------------------
  # Determine initial ranges for r
  #----------------------------------------------------
  # initial range of r from input file
  if(is.na(r.low)==F & is.na(r.hi)==F) {
    prior.r <- c(r.low,r.hi)
  } else {
    # initial range of r based on resilience
    if(res == "High") {
      prior.r <- c(0.6,1.5)} else if(res == "Medium") {
        prior.r <- c(0.2,0.8)}    else if(res == "Low") {
          prior.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
            prior.r <- c(0.015,0.1)} 
  }
 
  #--------------------------------------
  # extract data on stock
  #--------------------------------------
  yr           <- as.numeric(cdat$Year[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr])
  nyr          <- length(yr) # number of years in the time series
  Bk.yr.i      <- which(yr==Bk.yr)
  cpue.raw     <- as.numeric(cdat$CPUE[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr])

  # get catch from full assessments or simulations if available, for comparison
  if(is.null(cdat$Catch[cdat$Stock==stock][1])==F && is.na(cdat$Catch[cdat$Stock==stock][1])==F){
      C.ass <- cdat$Catch[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr] } else {C.ass <- NA}
  
  # apply correction for effort-creep to commercial(!) CPUE if indicated by user
  if(is.na(e.creep)==FALSE) {
    cpue.cor         <- cpue.raw
    for(i in 1:(length(cpue.raw)-1)) {
      cpue.cor[i+1]  <- cpue.raw[i+1]*(1-e.creep/100)^i # equation for decay in %; first cpue without correction
    }
    if(creep.graph==TRUE) {
      windows(8,6)
      plot(x=yr,y=cpue.raw,ylim=c(0,max(cpue.raw)),type="l",bty="l",xlab="Year",ylab="CPUE")
      lines(x=yr,y=cpue.cor,col="red")
      text(x=yr[length(yr)/2],y=max(cpue.raw),paste(stock," CPUE corrected for effort creep of ",e.creep," %",sep=""),col="red")
    }
    cpue.raw <- cpue.cor
  }
  
  d.cpue.raw   <- max(diff(cpue.raw)/cpue.raw[1:(nyr-1)])
    if(smooth.cpue==T||d.cpue.raw > 1.5) {
      smooth.flag <- TRUE
       bw          <- log(2)/exp(mean(log(prior.r))) # use population doubling time as bandwidth
       bw          <- ifelse(bw < 3,3,bw) # enforce minimum bandwidth of 3
       cpue        <- ksmooth(x=yr,y=cpue.raw,kernel="normal",n.points=length(yr),bandwidth=bw)$y  } else {
        smooth.flag <- FALSE
        cpue <- cpue.raw }
  # use median of 3 largest cpue as max cpue
  max.cpue     <- sort(cpue)[length(cpue)-1] 
  min.cpue     <- min(cpue)
  if(length(Fmsy.ass)>0 && is.null(Fmsy.ass)==F && is.na(Fmsy.ass)==F && is.null(cdat$F[1])==F) {
    FFmsy.ass     <- as.numeric(cdat$F[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Fmsy.ass)
    if(is.null(cdat$FLow[1])==F) {FFmsy.ass.lcl <- as.numeric(cdat$FLow[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Fmsy.ass) }
    if(is.null(cdat$FHi[1])==F) {FFmsy.ass.ucl <- as.numeric(cdat$FHi[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Fmsy.ass) }
  } else {Fmsy.ass <- NA; FFmsy.ass <- NA;FFmsy.ass.lcl <- NA;FFmsy.ass.ucl <- NA}

  if(length(Bmsy.ass)>0 && is.null(Bmsy.ass)==F && is.na(Bmsy.ass)==F && is.null(cdat$B[1])==F) {
    BBmsy.ass     <- as.numeric(cdat$B[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Bmsy.ass)
    if(is.null(cdat$BLow[1])==F) {BBmsy.ass.ucl <- as.numeric(cdat$BLow[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Bmsy.ass) }
    if(is.null(cdat$BHi[1])==F) {BBmsy.ass.lcl <- as.numeric(cdat$BHi[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr]/Bmsy.ass) }
  } else {Bmsy.ass <- NA; BBmsy.ass <- NA;BBmsy.ass.lcl <- NA;BBmsy.ass.ucl <- NA}

  if(length(yr)==0){
    cat("ERROR: Could not find the stock in the CPUE input files - Please check that the stock ID is written correctly")
    return (NA) }
  if(length(yr) != (end.yr-start.yr+1)) {
    cat("ERROR: indicated year range is of different length than years in CPUE file\n")
    return (NA)}

  #----------------------------------------------------
  # Determine ranges for relative biomass
  #----------------------------------------------------
  # initial range of B/k from input file
  if(is.na(Bk.pr.low)==F & is.na(Bk.pr.hi)==F) {
    prior.Bk <- c(Bk.pr.low,Bk.pr.hi)
  } else {
    if(Bk.pr == "Near unexploited") {
      prior.Bk <- c(0.75,1.0)} else if(Bk.pr == "More than half") {
        prior.Bk <- c(0.5,0.85)} else if(Bk.pr == "About half") {
          prior.Bk <- c(0.35,0.65)}    else if(Bk.pr == "Small") {
            prior.Bk <- c(0.15,0.4)}  else { # i.e. Bk.pr== "Very small"
              prior.Bk <- c(0.01,0.2)} 
  }  
  # us relative range of B/k as relative range for kq
  mean.prior.Bk   <- mean(prior.Bk)
  rr.prior.Bk     <- (mean.prior.Bk-prior.Bk[1])/mean.prior.Bk 
 
  prior.kq.low.1  <- (1-rr.prior.Bk)*cpue[Bk.yr.i]/mean.prior.Bk
  prior.kq.hi.1   <- (1+rr.prior.Bk)*cpue[Bk.yr.i]/mean.prior.Bk

   # kq must be > max cpue unless near unexploited
    prior.kq.low.2  <- ifelse(prior.kq.low.1 < max.cpue,ifelse(mean.prior.Bk >= 0.85,0.9*max.cpue,max.cpue),prior.kq.low.1) 
  
   # increase lower kq prior if cpue is small and flat
    if((max(cpue)/min(cpue))<2) {
      prior.kq.low <- ifelse(mean.prior.Bk < 0.3,2*prior.kq.low.2,
                        ifelse(mean.prior.Bk < 0.6,1.5*prior.kq.low.2,prior.kq.low.2))
      } else {prior.kq.low <- prior.kq.low.2 }
  
    # kq.hi at least 30-50% larger than kq.low, depending on Bk prior  
    if(mean.prior.Bk >= 0.6) {
      prior.kq.hi.2  <- ifelse(prior.kq.hi.1 < (1.3*prior.kq.low),1.3*prior.kq.low,prior.kq.hi.1) } else {
        prior.kq.hi.2  <- ifelse(prior.kq.hi.1 < (1.5*prior.kq.low),1.5*prior.kq.low,prior.kq.hi.1) }
  
   # if upper prior kq is too hi, limit to 3 times lower range 
      prior.kq.hi   <- ifelse(prior.kq.hi.2 > (3*prior.kq.low),3*prior.kq.low,prior.kq.hi.2)
     
   prior.kq           <- c(prior.kq.low,prior.kq.hi)
  
  #------------------------------------------------------------------
  # Sampling of r-k space 
  #------------------------------------------------------------------
  # turn numerical ranges into log-normal distributions 
  
    mean.log.r=mean(log(prior.r))
    sd.log.r=(log(prior.r[2])-log(prior.r[1]))/4  # assume range covers 4 SD
   
    mean.log.kq <- mean(log(prior.kq))
    sd.log.kq   <- (log(prior.kq[2])-log(prior.kq[1]))/4 # assume range covers 4 SD
 
    mvn.log.rk <- mvn(n=n,mean.log.r=mean.log.r,sd.log.r=sd.log.r,mean.log.kq=mean.log.kq,sd.log.kq=sd.log.kq)
    ri1    <- exp(mvn.log.rk[,1])
    kqi1   <- exp(mvn.log.rk[,2])
   
  #------------------------------------------------------------------
  # print prior info on screen
  #------------------------------------------------------------------
  cat("CPUE data for years ",yr[1]," - ",yr[nyr],", CPUE range ",min.cpue," - ",max(cpue),", smooth = ",smooth.flag,sep="","\n")
  cat("Prior for r                   = ",res,", ", r.low," - ",r.hi,sep="","\n")
  if(is.na(r.low)==T) {
     cat("Used prior range for r        = ", prior.r[1]," - ",prior.r[2],sep="","\n") } else {
       cat("Used prior range for r        = ", quantile(ri1,0.01)," - ",quantile(ri1,0.99),sep="","\n") }
  cat("Prior for ",Bk.yr," stock status   = ", Bk.pr,", ",Bk.pr.low," - ",Bk.pr.hi,sep="","\n") 
  cat("Used ",Bk.yr," prior B/B0 range    = ",prior.Bk[1]," - ",prior.Bk[2],", prior B/Bmsy = ",2*prior.Bk[1]," - ",2*prior.Bk[2],sep="","\n")
  cat("Used prior range for kq       = ",prior.kq[1]," - ",prior.kq[2]," [original range = ",prior.kq.low.1," - ",prior.kq.hi.1,"]\n",sep="") 
  if(is.na(Fmsy.ass)==F) {cat("Assessment Fmsy               =",Fmsy.ass,"\n")}
  if(is.na(FFmsy.ass[1])==F) {cat("Assessment F/Fmsy             =",FFmsy.ass[nyr-1],ifelse(is.na(FFmsy.ass.lcl[1])==F,
                            paste(",",format(FFmsy.ass.lcl[nyr-1],digits = 2),"-",format(FFmsy.ass.ucl[nyr-1],digits=2),"")),
                            "(",yr[nyr-1],")\n")}
  if(is.na(Bmsy.ass)==F) {cat("Assessment proxy Bmsy         =",Bmsy.ass,"\n")}
  if(is.na(BBmsy.ass[1])==F) {cat("Assessment proxy B/Bmsy       =",BBmsy.ass[nyr],ifelse(is.na(BBmsy.ass.lcl[1])==F,
                            paste(",",format(BBmsy.ass.lcl[nyr],digits = 2),"-",format(BBmsy.ass.ucl[nyr],digits=2),"")),
                            "(",yr[nyr],")\n")
                            cat("Source:",source,"\n")}
  
  if(is.null(cinfo$MSY.BSM[cinfo$Stock==stock])==F && is.na(cinfo$MSY.BSM[cinfo$Stock==stock])==F) { # assume all info from BSM analysis is available
    cat("BSM r                         =", cinfo$r.BSM[cinfo$Stock==stock],",",cinfo$r.BSM.lcl[cinfo$Stock==stock],
        "-",cinfo$r.BSM.ucl[cinfo$Stock==stock],"\n")
    cat("BSM k                         =", cinfo$k.BSM[cinfo$Stock==stock]*1000,",",cinfo$k.BSM.lcl[cinfo$Stock==stock]*1000,
        "-",cinfo$k.BSM.ucl[cinfo$Stock==stock]*1000,"\n")
    cat("BSM MSY                       =", cinfo$MSY.BSM[cinfo$Stock==stock]*1000,",",cinfo$MSY.BSM.lcl[cinfo$Stock==stock]*1000,
        "-",cinfo$MSY.BSM.ucl[cinfo$Stock==stock]*1000, "\n")
    cat("BSM last B/Bmsy               =", cinfo$B_Bmsy[cinfo$Stock==stock],",",cinfo$B_Bmsy.lcl[cinfo$Stock==stock],
        "-",cinfo$B_Bmsy.ucl[cinfo$Stock==stock], "\n")
    cat("BSM last F/Fmsy               =", cinfo$F_Fmsy[cinfo$Stock==stock],",",cinfo$F_Fmsy.lcl[cinfo$Stock==stock],
        "-",cinfo$F_Fmsy.ucl[cinfo$Stock==stock], "\n")
  }
  cat("Comment:",comment,"\n")
  if(is.na(source)==F) {cat("Source:",source,"\n")}
  cat("\n")
  #-----------------------------------------------------------------
  # Plot CPUE data and prior CPUE_msy
  #-----------------------------------------------------------------
  if(close.plots==T) {graphics.off()} # close previous plots, e.g. in batch processing
  if(do.plots==T) {
  
  # check for operating system, open separate window for graphs if Windows
  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(14,9)}
  
  par(mfrow=c(2,3))
  # (a): plot CPUE with prior for CPUEmsy
  plot(x=yr, y=cpue, 
       ylim=c(0,max(ifelse(substr(id_file,1,3)=="Sim",1.1*cinfo$true.MSYq,0),1.2*max(cpue),
                    1.2*max(cpue.raw),prior.kq[2])),
       type ="l", bty="l", main=paste("(a) CPUE", stock), xlab="Year", 
       ylab="CPUE", lwd=0.5, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  lines(x=yr,y=cpue.raw,lwd=2)
  # arrows(x0=Bk.yr,y0=(1+cv.prior.Bk)*cpue[Bk.yr.i],x1=Bk.yr,y1=(1-cv.prior.Bk)*cpue[Bk.yr.i],length=0.05,angle = 90, code=3)# plot prior range for B/k
  arrows(x0=Bk.yr,y0=prior.kq[1],x1=Bk.yr,y1=prior.kq[2],length=0.05,angle = 90, code=3,col="blue")# plot prior range for B/k

  lines(x=c(yr[1],yr[nyr]),y=c(prior.kq.low/2,prior.kq.low/2),lty="dotted")
  lines(x=c(yr[1],yr[nyr]),y=c(prior.kq.hi/2,prior.kq.hi/2),lty="dotted")
  text(x=yr[nyr-as.integer(0.15*nyr)],y=(prior.kq.low+prior.kq.hi)/4,"Bmsy_q")

  # (b): plot r-k graph 
   plot(x=ri1, y=kqi1, xlim = c(0.95*quantile(ri1,0.001),1.2*quantile(ri1,0.999)), 
                       ylim = c(0.95*quantile(kqi1,0.001),1.2*quantile(kqi1,0.999)),
        log="xy", xlab="r", ylab="kq", main="(b) Finding viable r-kq", pch=".", cex=2, bty="l", 
        col="grey90", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
   lines(x=c(prior.r[1],prior.r[2],prior.r[2],prior.r[1],prior.r[1]), # plot original prior range
         y=c(prior.kq[1],prior.kq[1],prior.kq[2],prior.kq[2],prior.kq[1]),
         lty="dotted") 
  } # end of do.plots loop
  
  #---------------------------------------------------------------------
  # Call AMSY-Schaefer function to filter r-kq space for viable r-kq pairs
  #---------------------------------------------------------------------
  cat("Monte Carlo filtering of r-kq space with",n,"points and",n.trial,"error patterns. \n")
  MCA1 <-  SchaeferCPUE(yr=yr, cpue=cpue, ri=ri1, kqi=kqi1, sigR=sigma.r[res.i], filter=filter) #><> correct PE input
 
  n.viable <- length(MCA1[,"rv"])
  cat("Viable r-kq pairs =",n.viable,"\n") 
  
  if((n.viable)<min.viable) {
    cat("Too few r-kq pairs after filtering, repeating analysis with 2 times more pairs, extended prior ranges, and increased smoothing:\n")
    mvn.log.rk2 <- mvn(n=2*n,mean.log.r=mean.log.r,sd.log.r=1.2*sd.log.r,mean.log.kq=mean.log.kq,sd.log.kq=1.2*sd.log.kq)
    ri2  <- exp(mvn.log.rk2[,1])
    kqi2 <- exp(mvn.log.rk2[,2])
    points(x=ri2, y=kqi2,pch=".",cex=2,col="grey90") # plot new potential r-k pairs
    lines(x=c(prior.r[1],prior.r[2],prior.r[2],prior.r[1],prior.r[1]), # re-plot original prior range
          y=c(prior.kq[1],prior.kq[1],prior.kq[2],prior.kq[2],prior.kq[1]),
          lty="dotted")
    points(MCA1,col="grey20", pch=".", cex=2.5) # re-plot points found so far
    text(x=prior.r[2],y=1.15*prior.kq[2],"extended")
    cpue <- ksmooth(x=yr,y=cpue.raw,kernel="normal",n.points=length(yr),bandwidth=5)$y
    
    MCA2 <-  SchaeferCPUE(yr=yr, cpue=cpue, ri=ri2, kqi=kqi2, sigR=sigma.r[res.i], filter=filter) 
    MCA  <- rbind(MCA1,MCA2)
    n.viable <- length(MCA[,"rv"])
    cat("Viable r-kq pairs =",n.viable,"\n") } else {MCA <- MCA1}
  
  if((n.viable)<10) {
    cat("Too few r-kq pairs after filtering, doing analysis without filters:\n")
    MCA <-  SchaeferCPUE(yr=yr, cpue=cpue, ri=ri2, kqi=kqi2, sigR=sigma.r[res.i], filter=FALSE) 
    text(x=prior.r[1],y=0.9*prior.kq[1],"no filters")
    }
  
  rv       <- MCA[,"rv"]
  kqv      <- MCA[,"kqv"]
  
  MSYqv     <- rv * kqv / 4
  MSYq.est  <- median(MSYqv)    
  MSYq.lcl  <- as.numeric(quantile(MSYqv,0.025))      
  MSYq.ucl  <- as.numeric(quantile(MSYqv,0.975))      
  
  n.v        <- length(MSYqv)
  
  kqv.est    <- median(kqv)
  kqv.lcl  <- as.numeric(quantile(kqv,0.025))
  kqv.ucl  <- as.numeric(quantile(kqv,0.975))

  rv.est     <- 4*MSYq.est/kqv.est    # rv corresponding to median(kqv)
  rv.lcl   <- as.numeric(quantile(rv,0.025))
  rv.ucl   <- as.numeric(quantile(rv,0.975))

  cqt.sel           <- matrix(nrow=length(rv),ncol=nyr-1)
  colnames(cqt.sel) <- c(yr[1:nyr-1])
  for(j in 1:(nyr-1)) {
    cqt.sel[,j]     <- MCA[,j+2]}
  cqt.median        <- apply(cqt.sel,2,median)
  cqt.lcl           <- apply(cqt.sel,2,quantile,probs=0.025)
  cqt.ucl           <- apply(cqt.sel,2,quantile,probs=0.975)
  
  cpuet.sel           <- matrix(nrow=length(rv),ncol=nyr)
  colnames(cpuet.sel) <- c(yr[1:nyr])
  for(j in 1:(nyr)) {
    cpuet.sel[,j]     <- MCA[,j+2+nyr-1]}
  cpuet.median        <- apply(cpuet.sel,2,median)
  cpuet.lcl           <- apply(cpuet.sel,2,quantile,probs=0.025)
  cpuet.ucl           <- apply(cpuet.sel,2,quantile,probs=0.975)
   
  BBmsy.end       <- cpuet.median[nyr]/(kqv.est/2)
  BBmsy.end.lcl   <- cpuet.lcl[nyr]/(kqv.est/2)
  BBmsy.end.ucl   <- cpuet.ucl[nyr]/(kqv.est/2)
  
  
  Ft            <- cqt.median[1:(nyr-1)]/cpuet.median[1:(nyr-1)]
  FFmsy         <- Ft/(0.5*rv.est)
  FFmsy.end     <- FFmsy[nyr-1]
  Ft.lcl        <- cqt.lcl[1:(nyr-1)]/cpuet.median[1:(nyr-1)]
  FFmsy.lcl     <- Ft.lcl/(0.5*rv.est)
  FFmsy.end.lcl <- FFmsy.lcl[nyr-1]
  Ft.ucl        <- cqt.ucl[1:(nyr-1)]/cpuet.median[1:(nyr-1)]
  FFmsy.ucl     <- Ft.ucl/(0.5*rv.est)
  FFmsy.end.ucl <- FFmsy.ucl[nyr-1]
  
  if(substr(id_file,1,3)=="Sim"){ # if dealing with simulated data, get the "true" values
   MSYq.true      <- cinfo$true.MSYq[cinfo$Stock==stock]
   r.true         <- as.numeric(cinfo$true.r[cinfo$Stock==stock])
   kq.true        <- as.numeric(cinfo$true.kq[cinfo$Stock==stock])
   cqt.true       <- as.numeric(cdat$Catch[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr])
   Ft.true        <- cqt.true/cpuet.median
   FFmsy.true     <- Ft.true/(0.5*cinfo$true.r[cinfo$Stock==stock])
   FFmsy.end.true <- FFmsy.true[nyr-1]
   BBmsy.end.true <- cinfo$true.Bk.end[cinfo$Stock==stock]*2
  } else {FFmsy.true <- NA; FFmsy.end.true <- NA; cqt.true <- NA}
  
  cat("\n Results:",
      "\n viable r-kq pairs   = ",n.v,
      "\n median kq           = ",kqv.est,", ",kqv.lcl," - ",kqv.ucl,
      "\n median MSYq         = ",MSYq.est,", ",MSYq.lcl," - ",MSYq.ucl,
      "\n r (4 MSYq/kq)       = ",rv.est,", ",rv.lcl," - ",rv.ucl,
      "\n Fmsy (r/2)          = ",rv.est/2,", ",rv.lcl/2," - ",rv.ucl/2,
      "\n F/Fmsy              = ",FFmsy.end,", ",FFmsy.end.lcl," - ",FFmsy.end.ucl," (",yr[nyr-1],")",
      ifelse(substr(id_file,1,3)=="Sim",paste(", true:",format(FFmsy.end.true,digits = 3)),""),
      "\n B/Bmsy              = ",BBmsy.end,", ",BBmsy.end.lcl," - ",BBmsy.end.ucl," (",yr[nyr],")",
          ifelse(substr(id_file,1,3)=="Sim",paste(", true:",format(BBmsy.end.true,digits = 3)),""),
      "\n",sep="")
 
  # -----------------------------------------
  # Plot results 
  # -----------------------------------------
  if(do.plots==T) {
  
  # (b): continued.. 
  # Add estimated r-k with confidence limits, for comparison with prior range
   points(x=rv.est,y=kqv.est,col="red3",lwd=2) # most probably r-kq pair estimate
   lines(x=c(rv.lcl,rv.ucl),y=c(kqv.est,kqv.est),col="red3",lwd=2) # confidence limits
   lines(x=c(rv.est,rv.est),y=c(kqv.lcl,kqv.ucl),col="red3",lwd=2) # confidence limits
  #  Add true r-kq to simulated data in plot B:
   if(substr(id_file,1,3)=="Sim") points(x=cinfo$true.r[cinfo$Stock==stock],y=cinfo$true.kq[cinfo$Stock==stock],col="blue",cex=2,lwd=2)
   if(substr(id_file,1,4)=="CMSY") points(x=cinfo$r.BSM[cinfo$Stock==stock],y=1000*cinfo$k.BSM[cinfo$Stock==stock],col="blue",cex=2,lwd=2)
    
  # (c): Analysis of viable r-k plot 
  # ----------------------------
    plot(x=MCA[,"rv"], y=MCA[,"kqv"], 
       xlim=c(min(c(rv,cinfo$true.r[cinfo$Stock==stock],rv.lcl),na.rm=T),max(c(rv,cinfo$true.r[cinfo$Stock==stock],rv.ucl),na.rm=T)),
       ylim=c(min(c(kqv,cinfo$true.kq[cinfo$Stock==stock],kqv.lcl),na.rm=T),max(c(kqv,cinfo$true.kq[cinfo$Stock==stock],kqv.ucl),na.rm=T)),
       pch=".", cex=2.5, col="grey20", log="xy", bty="l",
       xlab="r", ylab="kq", main="(c) Analysis of viable r-kq",  cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
 
#  points(x=rv,y=kqv,pch=".", cex=2.5, col="indianred1") # selected points for deriving estimates
  points(x=rv.est,y=kqv.est,col="red3",lwd=2) # most probably r-kq pair estimate
  lines(x=c(rv.lcl,rv.ucl),y=c(kqv.est,kqv.est),col="red3",lwd=2) # confidence limits
  lines(x=c(rv.est,rv.est),y=c(kqv.lcl,kqv.ucl),col="red3",lwd=2) # confidence limits
  if(substr(id_file,1,3)=="Sim") points(x=cinfo$true.r[cinfo$Stock==stock],y=cinfo$true.kq[cinfo$Stock==stock],col="blue",cex=2,lwd=2) # true r-kq pair

  # plot r-k from BSM if available
  if(is.null(cinfo$r.BSM[cinfo$Stock==stock])==F) {
    abline(v=cinfo$r.BSM[cinfo$Stock==stock],lty="dashed",col="blue")
  }
  
  # plot r-k confidence limits if CPUE test file is available
  if(substr(id_file,1,4)=="CMSY") {
    points(x=cinfo$r.BSM[cinfo$Stock==stock],y=cinfo$k.BSM[cinfo$Stock==stock]*1000,col="blue",lwd=2) # most probably r-kq pair estimate
    lines(x=c(cinfo$r.BSM.lcl[cinfo$Stock==stock],cinfo$r.BSM.ucl[cinfo$Stock==stock]),
          y=c(cinfo$k.BSM[cinfo$Stock==stock]*1000,cinfo$k.BSM[cinfo$Stock==stock]*1000),col="blue",lwd=2) # confidence limits
    lines(x=c(cinfo$r.BSM[cinfo$Stock==stock],cinfo$r.BSM[cinfo$Stock==stock]),
          y=c(cinfo$k.BSM.lcl[cinfo$Stock==stock]*1000,cinfo$k.BSM.ucl[cinfo$Stock==stock]*1000),col="blue",lwd=2) # confidence limits
  }

  # (d) Pred. rel. catch plot 
  #--------------------
  # get data from full assessments if available
  C_MSY.ass <- NA
  if(is.na(C.ass[1])==F && is.na(Fmsy.ass)==F && is.na(Bmsy.ass)==F) {
    MSY.ass   <- Fmsy.ass*Bmsy.ass
    C_MSY.ass <- C.ass/MSY.ass
  }
  # get data from BSM if available
  if(is.null(cinfo$MSY.BSM[cinfo$Stock==stock])==F) { 
    C_MSY.ass <- C.ass/(cinfo$MSY.BSM[cinfo$Stock==stock]*1000)
  }
  
  # determine height of y-axis in plot
  max.y  <- max(c(cqt.ucl/MSYq.est,1.1,cqt.true/MSYq.est,C_MSY.ass), na.rm=T)
  # Main plot of relative CMSY catch up to last year because the schaefer equation is not reliable in nyr
  plot(x=yr[1:(nyr-1)],y=cqt.median[1:(nyr-1)]/MSYq.est, lwd=2, xlab="Year", ylab="C/MSY", type="l",
       ylim=c(0,max.y), bty="l", main="(d) Catch/MSY", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  lines(x=yr[2:(nyr-1)], y=cqt.lcl[2:(nyr-1)]/MSYq.est,type="l",lty="dotted")
  lines(x=yr[2:(nyr-1)], y=cqt.ucl[2:(nyr-1)]/MSYq.est,type="l",lty="dotted")
  lines(x=c(yr[1],yr[nyr-1]),y=c(1,1), lty="dashed", lwd=1.5)   # line indicating MSY
  text(x=yr[nyr-as.integer(0.1*nyr)],y=1.1, "MSY")
 
   # plot true catch * q from simulations
  if(substr(id_file,1,3)=="Sim") {
    cqt.true <- as.numeric(cdat$Catch[cdat$Stock==stock & cdat$Year >= start.yr & cdat$Year <= end.yr])
    lines(x=yr[1:(nyr-1)], y=cqt.true[1:(nyr-1)]/MSYq.true,col="blue") } 

  # plot catch/MSY from BSM of full assessments
  if(is.na(C.ass[1])==F && is.null(cinfo$r.BSM[cinfo$Stock==stock])==F) {
    lines(x=yr[1:(nyr-1)], y=C_MSY.ass[1:(nyr-1)],col="blue")
  } 

# plot (e): F/Fmsy
#---------------
  max.y <- max(c(1.2,FFmsy.ucl,FFmsy.true,FFmsy.ass.ucl),na.rm=T)
  plot(x=yr[1:(nyr-1)],y=FFmsy, 
       ylim=c(0,max.y), 
       lwd=2, xlab="Year", ylab="F/Fmsy", type="l",
       bty="l", main="(e) F/Fmsy", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  lines(x=yr[2:(nyr-1)],y=FFmsy.lcl[2:(nyr-1)],lty="dotted") # lcl of F/Fmsy
  lines(x=yr[2:(nyr-1)],y=FFmsy.ucl[2:(nyr-1)],lty="dotted") # ucl of F/Fmsy
  lines(x=c(yr[1],yr[nyr-1]),y=c(1,1),lty="dashed") # indicates F/Fmsy = 1
  if(substr(id_file,1,3)=="Sim") {lines(x=yr[1:(nyr-1)],y=FFmsy.true[1:(nyr-1)], col="blue")} # true F/Fmsy from simulations
  if(is.na(FFmsy.ass[1])==F) {lines(x=yr[1:(nyr-1)],y=FFmsy.ass[1:(nyr-1)], col="blue")} # F/Fmsy from assessments
  if(is.na(FFmsy.ass.ucl[1])==F) {lines(x=yr[1:(nyr-1)],y=FFmsy.ass.ucl[1:(nyr-1)], lty="dotted",col="blue")} # F/Fmsy from assessments
  if(is.na(FFmsy.ass.lcl[1])==F) {lines(x=yr[1:(nyr-1)],y=FFmsy.ass.lcl[1:(nyr-1)], lty="dotted",col="blue")} # F/Fmsy from assessments
  
  # plot F/Fmsy from BSM
  if(is.null(cinfo$F_Fmsy[cinfo$Stock==stock])==F) {
    points(x=yr[nyr-1],y=cinfo$F_Fmsy[cinfo$Stock==stock],col="blue")
    lines(x=c(yr[nyr-1],yr[nyr-1]),y=c(cinfo$F_Fmsy.lcl[cinfo$Stock==stock],
                                   cinfo$F_Fmsy.ucl[cinfo$Stock==stock]),col="blue")
  }
  
  # plot F/Fmsy from BSM if CMSY file
  if(is.na(C.ass[1])==F && is.null(cinfo$r.BSM[cinfo$Stock==stock])==F && substr(id_file,1,4)=="CMSY") {
    F.BSM      <- C.ass[1:(nyr-1)] / cpue[1:(nyr-1)]
    FFmsy.BSM  <- F.BSM / (cinfo$r.BSM[cinfo$Stock==stock]/2)
    lines(x=yr[1:(nyr-1)], y=FFmsy.BSM[1:(nyr-1)],col="blue")
  }
    
# plot (f): B/Bmsy
#---------------
  Bkt        <- cpuet.median/(0.5*kqv.est)
  Bmsy.true  <- cinfo$true.kq[cinfo$Stock==stock]/kqv.est
  if(is.na(C.ass[1])==F && is.null(cinfo$r.BSM[cinfo$Stock==stock])==F && substr(id_file,1,4)=="CMSY") {
    BBmsy.BSM  <- cpue.raw / (cinfo$k.BSM[cinfo$Stock==stock]*1000 / 2) } else { BBmsy.BSM <- NA }
  max.y      <- max(c(Bkt, kqv.ucl/kqv.est,Bmsy.true,cpue.raw/(kqv.est/2),cpuet.ucl/(0.5*kqv.est),BBmsy.ass.ucl,BBmsy.BSM),na.rm=T)
  
  plot(x=yr,y=Bkt,type="l",ylim=c(0,max.y),
       bty="l",lwd=2,xlab="Year",ylab="B/Bmsy",
       main="(f) B/Bmsy", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
  lines(x=yr[2:nyr],y=cpuet.lcl[2:nyr]/(0.5*kqv.est),lty="dotted")
  lines(x=yr[2:nyr],y=cpuet.ucl[2:nyr]/(0.5*kqv.est),lty="dotted")
  lines(x=c(yr[1],yr[nyr]),y=c(1,1),lty="dashed")
  lines(x=c(yr[1],yr[nyr]),y=c(kqv.lcl/kqv.est,kqv.lcl/kqv.est),lty="dotted")
  lines(x=c(yr[1],yr[nyr]),y=c(kqv.ucl/kqv.est,kqv.ucl/kqv.est),lty="dotted")
  lines(x=c(yr[1],yr[nyr]),y=c(0.5,0.5),lty="longdash",col="red")
  if(substr(id_file,1,3)=="Sim") {
    lines(x=c(yr[1],yr[nyr]),y=c(Bmsy.true,Bmsy.true),lty="dashed",col="blue")
    lines(x=yr,y=cpue.raw/(cinfo$true.kq[1]/2),col="blue")}
  if(is.na(BBmsy.ass[1])==F) {lines(x=yr,y=BBmsy.ass,col="blue")}
  if(is.na(BBmsy.ass.ucl[1])==F) {lines(x=yr,y=BBmsy.ass.ucl,lty="dotted",col="blue")}
  if(is.na(BBmsy.ass.lcl[1])==F) {lines(x=yr,y=BBmsy.ass.lcl,lty="dotted",col="blue")}
  
  # plot B/Bmsy from BSM
  if(is.null(cinfo$B_Bmsy[cinfo$Stock==stock])==F) {
    points(x=yr[nyr],y=cinfo$B_Bmsy[cinfo$Stock==stock],col="blue")
    lines(x=c(yr[nyr],yr[nyr]),y=c(cinfo$B_Bmsy.lcl[cinfo$Stock==stock],
                                   cinfo$B_Bmsy.ucl[cinfo$Stock==stock]),col="blue")
  }
  # plot B/Bmsy from BSM if CMSY file
  if(is.na(C.ass[1])==F && is.null(cinfo$r.BSM[cinfo$Stock==stock])==F && substr(id_file,1,4)=="CMSY") {
    lines(x=yr, y=BBmsy.BSM,col="blue")
  }
  
  #if(is.na(C.ass[1])==F && is.null(cinfo$r.BSM[cinfo$Stock==stock])==F) {
   # lines(x=yr, y=BBmsy.BSM,col="blue")
  #}

  if (save.plots==TRUE & do.plots==TRUE) {
    jpgfile<-paste(stock,"_AMSY.jpg",sep="")
    dev.copy(jpeg,jpgfile,
             width = 1024, 
             height = 768, 
             units = "px", 
             pointsize = 18,
             quality = 95,
             res=80,
             antialias="cleartype")
    dev.off()
  }
  
  if(kobe.plot==T){
    # open window for plot of four panels
    if(grepl("win",tolower(Sys.info()['sysname']))) {windows(7,7)}
    par(mfrow=c(1,1))  
    # make margins narrower
    par(mar=c(4.1,4.1,2.1,2.1))
    
    
      bbmsy = (cpuet.sel[,nyr]/(0.5*kqv.est))
      ffmsy = ((apply(cqt.sel[,(nyr-4):(nyr-1)],1,median)/cpuet.sel[,nyr])/(0.5*rv.est))
      log.bbmsy = log(bbmsy[ffmsy>0]) # Prevents NA warning
      log.ffmsy = log(ffmsy[ffmsy>0]) # Prevents NA warning
      
      # get mean after all the CMSY subsetting (can't match with biomass sbmsetting)
      mu.kobe = c(median(log.ffmsy),median(log.bbmsy))
      # Get covariance of the 2 vectors
      cov.kobe = cov(cbind(log.ffmsy,log.bbmsy)) 
      # Generate 10000 new random deviates from a MVN
      log.kobe.mvn = rmvnorm(10000 ,mean = mu.kobe,sigma = cov.kobe)
      kobe.mvn = exp(log.kobe.mvn)
      # Generate 10000 new random deviates from a MVN
      x.F_Fmsy =exp(log.kobe.mvn[,1])
      y.b_bmsy =exp(log.kobe.mvn[,2])
    
    kernelF <- ci2d(y.b_bmsy,x.F_Fmsy,nbins=151,factor=2.2,ci.levels=c(0.50,0.80,0.75,0.90,0.95),show="none",col=1,xlab= ifelse(harvest.label=="Fmsy",expression(paste(F/F[MSY])),expression(paste(H/H[MSY]))),ylab=expression(paste(B/B[MSY])))
    
    max.y   <- max(c(2, quantile(x.F_Fmsy,0.96),na.rm =T))
    max.x    <- max(max(2,quantile(y.b_bmsy,0.999)))
    
    
    #Create plot
    plot(1000,1000,type="b", xlim=c(0,max.x), ylim=c(0,max.y),lty=3,xlab="",ylab="", bty="l",  cex.main = 2, cex.lab = 1.35, cex.axis = 1.35,xaxs = "i",yaxs="i")
    mtext(expression(B/B[MSY]),side=1, line=2.7, cex=1.4)
    mtext(expression(F/F[MSY]),side=2, line=2.7, cex=1.4)
    c1 <- c(-1,100)
    c2 <- c(1,1)
    
    # extract interval information from ci2d object
    # and fill areas using the polygon function
    zb2 = c(0,1)
    zf2  = c(1,100)
    zb1 = c(1,100)
    zf1  = c(0,1)
    polygon(c(zb1,rev(zb1)),c(0,0,1,1),col="green",border=0)
    polygon(c(zb2,rev(zb2)),c(0,0,1,1),col="yellow",border=0)
    polygon(c(1,100,100,1),c(1,1,100,100),col="orange",border=0)
    polygon(c(0,1,1,0),c(1,1,100,100),col="red",border=0)
    
    polygon(kernelF$contours$"0.95",lty=2,border=NA,col="cornsilk4")
    polygon(kernelF$contours$"0.8",border=NA,lty=2,col="grey")
    polygon(kernelF$contours$"0.5",border=NA,lty=2,col="cornsilk2")
    points(Bkt,c(FFmsy,median(x.F_Fmsy)),pch=16,cex=1)
    
    lines(c1,c2,lty=3,lwd=0.7)
    lines(c2,c1,lty=3,lwd=0.7)
    lines(Bkt,c(FFmsy,median(x.F_Fmsy)), lty=1,lwd=1.)
    points(Bkt[1],c(FFmsy,median(x.F_Fmsy))[1],col=1,pch=22,bg="white",cex=1.5)
    points(Bkt[which(yr==median(yr))],c(FFmsy,median(x.F_Fmsy))[which(yr==median(yr))],col=1,pch=21,bg="white",cex=1.5)
    points(Bkt[nyr],c(FFmsy,median(x.F_Fmsy))[nyr],col=1,pch=24,bg="white",cex=1.5)
    # Get Propability
    Pr.green = sum(ifelse(y.b_bmsy>1 & x.F_Fmsy<1,1,0))/length(y.b_bmsy)*100
    Pr.red = sum(ifelse(y.b_bmsy<1 & x.F_Fmsy>1,1,0))/length(y.b_bmsy)*100
    Pr.yellow = sum(ifelse(y.b_bmsy<1 & x.F_Fmsy<1,1,0))/length(y.b_bmsy)*100
    Pr.orange = sum(ifelse(y.b_bmsy>1 & x.F_Fmsy>1,1,0))/length(y.b_bmsy)*100
    
    sel.years = c(median(yr))
    
    legend('topright', 
           c(paste(start.yr),paste(median(yr)),paste(end.yr),"50% C.I.","80% C.I.","95% C.I.",paste0(round(c(Pr.red,Pr.yellow,Pr.orange,Pr.green),1),"%")), 
           lty=c(1,1,1,rep(-1,8)),pch=c(22,21,24,rep(22,8)),pt.bg=c(rep("white",3),"cornsilk2","grey","cornsilk4","red","yellow","orange","green"), 
           col=1,lwd=1.1,cex=1.1,pt.cex=c(rep(1.3,3),rep(1.7,3),rep(2.2,4)),bty="n",y.intersp = 1.)  
  } # End Kobe
    
    if (save.plots==TRUE & kobe.plot==TRUE) 
    {
      jpgfile<-paste(stock,"_KOBE.jpg",sep="")
      dev.copy(jpeg,jpgfile,
               width = 1024*0.7, 
               height = 1024*0.7, 
               units = "px", 
               pointsize = 18,
               quality = 95,
               res=80,
               antialias="cleartype")
      dev.off()
    }  

  FFmsy.retrospective[[retrosp.step+1]]<-FFmsy #retrospective analysis
  BBmsy.retrospective[[retrosp.step+1]]<-Bkt #retrospective analysis
  years.retrospective[[retrosp.step+1]]<-yr #retrospective analysis
 
  } # end of do.plots loop
} #retrospective analysis - end loop

#retrospective analysis plots
if (retros.nyears>0 && do.plots==T){

  if(grepl("win",tolower(Sys.info()['sysname']))) {windows(12,7)}
    par(mfrow=c(1,2))  
  
	allyears<-years.retrospective[[1]]
	nyrtotal<-length(allyears)
	plot(x=allyears[1:(nyrtotal-1)],y=FFmsy.retrospective[[1]], main=as.character(stock), ylim=c(0,max(FFmsy.retrospective[[1]],na.rm=T)), lwd=2, xlab="Year", ylab="F/Fmsy", type="l", bty="l",  cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
	lines(x=allyears[1:(nyrtotal-2)],y=FFmsy.retrospective[[2]], type = "o", pch=15, col="red")
	lines(x=allyears[1:(nyrtotal-3)],y=FFmsy.retrospective[[3]], type = "o", pch=16, col="green") 
	lines(x=allyears[1:(nyrtotal-4)],y=FFmsy.retrospective[[4]], type = "o", pch=17, col="blue") 
	legend("bottomleft", legend = c("Reference",allyears[nyrtotal-2],allyears[nyrtotal-3],allyears[nyrtotal-4]), 
       col=c("black","red", "green", "blue"), lty=1, pch=c(-1,15,16,17))

	plot(x=allyears[1:(nyrtotal)],y=BBmsy.retrospective[[1]],main=as.character(stock), ylim=c(0,max(BBmsy.retrospective[[1]],na.rm=T)), lwd=2, xlab="Year", ylab="B/Bmsy", type="l", bty="l",cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
	lines(x=allyears[1:(nyrtotal-1)],y=BBmsy.retrospective[[2]], type = "o", pch=15, col="red")
	lines(x=allyears[1:(nyrtotal-2)],y=BBmsy.retrospective[[3]], type = "o", pch=16, col="green") 
	lines(x=allyears[1:(nyrtotal-3)],y=BBmsy.retrospective[[4]], type = "o", pch=17, col="blue") 
	legend("bottomleft", legend = c("Reference",allyears[nyrtotal-1],allyears[nyrtotal-2],allyears[nyrtotal-3]), 
		col=c("black","red", "green", "blue"), lty=1, pch=c(-1,15,16,17))
 } #retrospective analysis plots - end


#---------------------------------------	
# write results to file
#---------------------------------------	
	if(write.output==TRUE && substr(id_file,1,3)=="Sim") {
	  output = data.frame(as.character(stock), 
	                      r.true,rv.est,rv.lcl,rv.ucl,
	                      kq.true,kqv.est,kqv.lcl,kqv.ucl,
	                      cinfo$true.MSYq[cinfo$Stock==stock],
	                      MSYq.est,MSYq.lcl,MSYq.ucl,
	                      FFmsy.end.true,
	                      FFmsy.end,FFmsy.end.lcl,FFmsy.end.ucl,
	                      BBmsy.end.true,
	                      BBmsy.end,BBmsy.end.lcl,BBmsy.end.ucl )
	  
	  write.table(output, file=outfile, append = T, sep = ",", 
	              dec = ".", row.names = FALSE, col.names = FALSE)
	} 
	if(write.output==TRUE && is.null(cinfo$MSY.BSM[cinfo$Stock==stock])==F) { # assuming all BSM fields are available
	  output = data.frame(as.character(stock), 
	                      cinfo$r.BSM[cinfo$Stock==stock],cinfo$r.BSM.lcl[cinfo$Stock==stock],cinfo$r.BSM.ucl[cinfo$Stock==stock],
	                      rv.est,rv.lcl,rv.ucl,
	                      cinfo$k.BSM[cinfo$Stock==stock]*1000,cinfo$k.BSM.lcl[cinfo$Stock==stock]*1000,cinfo$k.BSM.ucl[cinfo$Stock==stock]*1000,
	                      kqv.est,kqv.lcl,kqv.ucl,
	                      cinfo$B_Bmsy[cinfo$Stock==stock],cinfo$B_Bmsy.lcl[cinfo$Stock==stock],cinfo$B_Bmsy.ucl[cinfo$Stock==stock],
	                      BBmsy.end,BBmsy.end.lcl,BBmsy.end.ucl,
	                      cinfo$F_Fmsy[cinfo$Stock==stock],cinfo$F_Fmsy.lcl[cinfo$Stock==stock],cinfo$F_Fmsy.ucl[cinfo$Stock==stock],
	                      FFmsy[nyr-1],FFmsy.lcl[nyr-1],FFmsy.ucl[nyr-1])
	  
	  write.table(output, file=outfile, append = T, sep = ",",dec = ".", row.names = FALSE, col.names = FALSE)
	} # end of BSM option to write results to file
	
	
	if(write.output==T && substr(id_file,1,3)!="Sim" && is.null(cinfo$MSY.BSM)==T){ # regular assessment
	  output = data.frame(as.character(stock),rv.est/2,rv.lcl/2,rv.ucl/2,
	                      kqv.est,kqv.lcl,kqv.ucl,
	                      FFmsy[nyr-1],FFmsy.lcl[nyr-1],FFmsy.ucl[nyr-1],
	                      cpue[nyr]/kqv.est,cpue[nyr]/kqv.ucl,cpue[nyr]/kqv.lcl)
	  write.table(output, file=outfile, append = T, sep = ",", 
	              dec = ".", row.names = FALSE, col.names = FALSE)
	} # end of regular output

} # end loop on stocks



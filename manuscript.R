# Title: Data Analysis for Raman & SERS characterization of Extracellular Vesicles 
# Date: 20.03.2019
# Author: Eric Boateng

# Load Packages

library(hyperSpec)
library(plotrix)
library(baseline)
library(viridisLite)


# set working directory and load temporary files
getwd () # get working directory
path <- "~/Downloads/IPHT-Data/Man"
setwd(path) # Set working directory and read/search for files with .spc extension
# Get the temporary working directory 
tmp <- Sys.glob ("~/Downloads/IPHT-Data/Man/Dried/*.[sS][pP][cC]")


################################################################################
# Function to read all list of hyperSpec temporary files                       #
for (i in 1:length(tmp)){                                                      #
  tmp[[i]] <-  list((paste(tmp[[i]])))                                         #
}                                                                              #
#                                                                              #
#                                                                              #
# Read all hyperSpec files in temporary folder                                 #
files <- sapply (tmp, read.spc.Kaiser)                                         #
                                                                               #
#Rename temporary Data                                                         #
                                                                               #
a1   <- files[1]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/a1_avg.SPC`     #
a2   <- files[2]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/a2_avg.SPC`     #
a4   <- files[3]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/a4_avg.SPC`     #
a5   <- files[4]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/a5_avg.SPC`     #
a7   <- files[5]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/a7_avg.SPC`     #
c1   <- files[7]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/c1_avg.SPC`     #
c2   <- files[8]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/c2_avg.SPC`     #
c4   <- files[9]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/c4_avg.SPC`     #
c5   <- files[10]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/c5_avg.SPC`    #
b7   <- files[6]$`/Users/boateng/Downloads/IPHT-Data/Man/Dried/b7_avg.SPC`     #

EV5         <- rbind (a1, c1)
EV12        <- rbind (a2, c2)
EV120       <- rbind (a4, c4)
EV120_PEG   <- rbind (a5, c5)
FC          <- rbind (a7, b7)






################################################################################
#                                                                              #
#                 Data Analysis for the Dried Samples                          #  
#                                                                              #
################################################################################
# Read all hyperSpec files in temporary folder 
files <-Sys.glob("~/Downloads/IPHT-Data/Man/Diff_dried/EV5/*.[sS][pP][cC]")
data <- sapply (files, read.spc.Kaiser)

# combine all spectra into one data file
data <- collapse (data)
data <- data[,, c(400~1720)]  # select only fingerprint region
labels(data, "spc") <- expression ("Raman Intensity / a.u")
labels(data, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data, col = matlab.dark.palette (nrow(data)),
   title.args= list(main = "FC"))     # plot Raw data



# Data Preprocessing
# Baseline correction # Alternative 1 (ALS)
baselineGUI(data[[]])
bsl <- baseline (data[[]], method = "als", lambda = 6, p=0.1)
spc.bl <- data - bsl@baseline
plot(spc.bl, zeroline=NA)    # plot subtracted baseline


# Alternative 2 (Modified polynomial fitting)
#baselineGUI(spc.out.cr[[]])
# bl <- baseline(data[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
# spc.bl <- data - bl@baseline
# plot(spc.bl, zeroline=NA)  # plot subtracted baseline



# Spectra smoothing (Noise Removal)
spc.sm <-spc.loess(spc.bl[,,], newx=400:1720)

# plot smoothed spectra
spc.identify(formatter = spc.label.wlonly,plot(spc.sm,col = c(1,2),zeroline=NA,
  title.args = list(main = "Raman spectra of EV5 [smoothed]")))

box()
#legend
legend("topleft", text.font=2, legend=c("Control","Cancer"),
  col = c(1,2), box.lty=1, lty=1,cex=1.1)


# Spectra Normalization (Normalization around amide I band)
spc.cor <- spc.sm / rowMeans (spc.sm [,, 1660])

# plot normalized spectra
plot(spc.cor, col = c(1,2),zeroline=NA,
  title.args=list(main="Raman spectra of BSA [Normalized]",
  ylab = expression(paste("Raman Intensity / arb. units")),
  xlab = expression(paste("Wavenumber /", cm^-1, ""))))


# plot spectra (overlay)
spc.identify(formatter = spc.label.wlonly,plot(spc.cor, col = c(1,2),zeroline=NA,
   title.args= list(main = "Raman spectra of dried FC", ylab="")))

box()
#legend
legend("topleft", text.font=2, legend=c("Control","Cancer"),
   col = c(1,2), box.lty=1, lty=1,cex=1.1)

# Calculating Difference spectra (control-cancer)
Spectrum_1 <- spc.cor [1]
Spectrum_2 <- spc.cor [2]
Diff <- Spectrum_1 - Spectrum_2
#plot (Spectrum_1)
#plot (Spectrum_2,add=T,col=2)

#plot Difference spectrum
spc.identify(formatter = spc.label.wlonly,plot(Diff,wl.range = c (400 ~ 1720),col=1,
   plot.args = list(ylim = c (-0.3, 0.3)), zeroline="3"))

box()

legend("bottomleft", text.font=2, legend=c("Difference"), 
   col = 1, box.lty=1, lty=1,cex=1.1)
################################################################################







################################################################################
#                                                                              #
#                Data Analysis for the Droplet Samples (SERS)                  #  
#                                                                              #
################################################################################

setwd("~/Downloads/IPHT-Data/Man/Droplet")
tmp <- Sys.glob ("~/Downloads/IPHT-Data/Man/Droplet/*.[sS][pP][cC]")
# Read all hyperSpec files in temporary folder                                 
files <- sapply (tmp, read.spc.Kaiser)                                         

#Rename temporary Data                                                         
a1   <- files[1]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/a1_avg.SPC`     
a2   <- files[2]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/a2_avg.SPC`     
a4   <- files[3]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/a4_avg.SPC`     
a5   <- files[4]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/a5_avg.SPC`     
a6   <- files[5]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/a6_avg.SPC`     
b1   <- files[6]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/b1_avg.SPC`     
b2   <- files[7]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/b2_avg.SPC`     
b4   <- files[8]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/b4_avg.SPC`     
b5   <- files[9]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/b5_avg.SPC`    
b6   <- files[10]$`/Users/boateng/Downloads/IPHT-Data/Man/Droplet/b6_avg.SPC`     


EV5         <- rbind (a1, b1)
EV12        <- rbind (a2, b2)
EV120       <- rbind (a4, b4)
EV120_PEG   <- rbind (a5, b5)
FC          <- rbind (a6, b6)


data <- EV5[,, c(400~1720)]
labels(data, "spc") <- expression ("Raman Intensity / a.u")
labels(data, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data, col = matlab.dark.palette (nrow (data)),
  title.args= list(main = "SERS spectra of EV5", 
  xlab = expression(paste("Wavenumber /", cm^-1, ""))))


# Baseline Correction (Modified polynomial fitting)
#baselineGUI(spc.out.cr[[]])
bl <- baseline(data[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl <- data - bl@baseline
plot(data)
plot(spc.bl, zeroline=NA)    # plot subtracted baseline


# Spectra smoothing (Noise Removal)
spc.sm <-spc.loess(spc.bl[,,], newx=400:1720)

# plot smoothed spectra
spc.identify(formatter = spc.label.wlonly,plot(spc.sm,stacked=T,col = c(1,2),zeroline=NA,
  title.args = list(main = "Raman spectrum of EV5 sample [smoothed]")))

#plot spectra (overlay)
spc.identify(formatter = spc.label.wlonly,plot(spc.sm, col = c(1,2),zeroline=NA,
   title.args= list(main = "SERS spectra of EV5")))

box()
#legend
legend("topleft", text.font=2, legend=c("Control","Cancer"),
       col = c(1,2), box.lty=1, lty=1,cex=1.1)

# Calculating Difference spectra (cancer-control)
Spectrum_1 <- Spectrum [1]
Spectrum_2 <- Spectrum [2]
Diff <- Spectrum_2 - Spectrum_1

# Spectrum <- rbind (Spectrum_1*3, Spectrum_2)
spc.identify(formatter = spc.label.wlonly,plot(Spectrum,zeroline=NA,col = c(1,2), stacked=T,wl.range = c (400 ~ 1720),
   title.args= list(main = "SERS spectra of EV5 droplet", xlab = expression(paste("Wavenumber /", cm^-1, "")))))
box()
legend("topleft", text.font=2, legend=c("Control","Cancer"), 
   col = c(1,2), box.lty=1, lty=1,cex=1.0)

#plot Difference spectrum
spc.identify(formatter = spc.label.wlonly,plot(Diff,wl.range = c (400 ~ 1720),col=1,
   zeroline="3" ,title.args= list(xlab = expression(paste("Wavenumber /", cm^-1, "")))))
legend("topleft", text.font=2, legend=c("Difference"), 
   col = 1, box.lty=1, lty=1,cex=1.0)

abline(v=c(713, 854, 1004, 1132, 1238, 1393),col="gray", lty=2) # add vertical gray dashed lines on indicated wavenumbers




################################################################################
#                                                                              #
#                 Mean & Std of Raman and SERS samples                         #  
#                                                                              #
################################################################################
# Read all hyperSpec files in temporary folder 
files1 <-Sys.glob("~/Downloads/IPHT-Data/Man/Droplet/a1/*.[sS][pP][cC]")
data1 <- sapply (files1, read.spc.Kaiser)

# combine all spectra into one data file
data1 <- collapse (data1)
data1 <- data1[,, c(350~1800)]  # select only fingerprint region
labels(data1, "spc") <- expression ("Raman Intensity / a.u")
labels(data1, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data1, col = matlab.dark.palette (nrow(data1)),
   title.args= list(main = "FC"))


# Baseline (Modified polynomial fitting)
bl1 <- baseline(data1[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl1 <- data1 - bl1@baseline
plot(spc.bl1, zeroline=NA)

# Spectra smoothing (Noise Removal)
A1corr <-spc.loess(spc.bl1[,,], newx=350:1800)
################################################################################
# Read all hyperSpec files in temporary folder 
files2 <-Sys.glob("~/Downloads/IPHT-Data/Man/Droplet/a2/*.[sS][pP][cC]")
data2 <- sapply (files2, read.spc.Kaiser)

# combine all spectra into one data file
data2 <- collapse (data2)
data2 <- data2[,, c(350~1800)]  # select only fingerprint region
labels(data2, "spc") <- expression ("Raman Intensity / a.u")
labels(data2, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data2, col = matlab.dark.palette (nrow(data2)),
   title.args= list(main = "FC"))   


# Baseline (Modified polynomial fitting)
bl2 <- baseline(data2[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl2 <- data2 - bl2@baseline
plot(spc.bl2, zeroline=NA)

# Spectra smoothing (Noise Removal)
A2corr <-spc.loess(spc.bl2[,,], newx=350:1800)
################################################################################
files3 <-Sys.glob("~/Downloads/IPHT-Data/Man/Droplet/a4/*.[sS][pP][cC]")
data3 <- sapply (files3, read.spc.Kaiser)

# combine all spectra into one data file
data3 <- collapse (data3)
data3 <- data3[,, c(350~1800)]  # select only fingerprint region
labels(data3, "spc") <- expression ("Raman Intensity / a.u")
labels(data3, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data3, col = matlab.dark.palette (nrow(data3)),
   title.args= list(main = "FC"))


# Baseline (Modified polynomial fitting)
bl3 <- baseline(data3[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl3 <- data3 - bl3@baseline
plot(spc.bl3, zeroline=NA)

# Spectra smoothing (Noise Removal)
A4corr <-spc.loess(spc.bl3[,,], newx=350:1800)
################################################################################
files4 <-Sys.glob("~/Downloads/IPHT-Data/Man/Droplet/a6/*.[sS][pP][cC]")
data4 <- sapply (files4, read.spc.Kaiser)

# combine all spectra into one data file
data4 <- collapse (data4)
data4 <- data4[,, c(350~1800)]  # select only fingerprint region
labels(data4, "spc") <- expression ("Raman Intensity / a.u")
labels(data4, ".wavelength") <- expression(paste("Wavenumber /", cm^-1, ""))
plotspc(data4, col = matlab.dark.palette (nrow(data4)),
   title.args= list(main = "FC"))


# Baseline (Modified polynomial fitting)
bl4 <- baseline(data4[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl4 <- data4 - bl4@baseline
plot(spc.bl4, zeroline=NA)

# Spectra smoothing (Noise Removal)
A7corr <-spc.loess(spc.bl4[,,], newx=350:1800)
################################################################################

# Normalization
A1corr <- A1corr / rowMeans (A1corr [,, 350 ~ 1800])
A2corr <- A2corr / rowMeans (A2corr [,, 350 ~ 1800])
A4corr <- A4corr / rowMeans (A4corr [,, 350 ~ 1800])
A7corr <- A7corr / rowMeans (A7corr [,, 350 ~ 1800])

################################################################################

A1corr <- A1corr[,,350~1800]
A2corr <- A2corr[,,350~1800]
A4corr <- A4corr[,,350~1800]
A7corr <- A7corr[,,350~1800]


# Collapse data
A1corr$sample <- factor ("EV5")
A2corr$sample <- factor ("EV12")
A4corr$sample <- factor ("EV120")
A7corr$sample <- factor ("FC")


all <- collapse(A1corr,A2corr,A4corr,A7corr)
spc.agg <- aggregate(all, by = list (all$sample), mean_pm_sd)
cols <- c(1,2,3,4)

# Normalization
#spc.cor <- spc.agg / rowMeans (spc.agg [,, 350 ~ 1800])

# plot sapectra mean +/- standard deviation 
spc.identify(formatter = spc.label.wlonly,plot(spc.agg, col= cols,stacked = ".aggregate",
   fill = ".aggregate",zeroline=NA, #wl.range = c (400 ~ 1720),
   title.args = list(main = "SERS spectra [Mean \u00B1 std] of control sample",
   ylab = expression(paste("")))))

# Vertical lines for SERS cancer
abline(v=c(382,394,524,600,713,854,1004,1132,1238,1277,1314,1393,1450,1560,1589,1615,1655),col="gray", lty=2)
# Vertical lines for SERS control
abline(v=c(382,394,478,564,622,656,713,803,854,899,1003,1030,1132,1236,1343,1393,1450,1589,1615,1653),col="gray", lty=2)
# Vertical lines for Raman control & cancer
abline(v=c(457,505,621,643,758,830,851,937,1003,1032,1126,1157,1207,1243,1315,1340,1450,1528,1555,1606,1657),col="gray", lty=2)

box()
legend("topleft",horiz = F, text.font=2, legend=c("EV5","EV12","EV120","FC"), 
   col = cols, box.lty=1, lty=1, cex=0.9)




################################################################################
#                                                                              #
#                 Quantifying SERS Difference spectra                          #  
#                                                                              #
################################################################################

setwd("~/Downloads/IPHT-Data/Man")  # set working directory

# Read all hyperSpec files in temporary folder 
files <-Sys.glob("~/Downloads/IPHT-Data/Man/Diff_Droplet/EV5/*.[sS][pP][cC]")
data  <- sapply (files, read.spc.Kaiser)
data  <- collapse (data)
data  <- data[,, c(350~1800)]
labels(data, "spc") <- expression ("Raman Intensity")
plotspc(data, col = matlab.dark.palette (nrow (data)),
   title.args= list(main = "EV5", 
   xlab = expression(paste("Wavenumber /", cm^-1, ""))))

# Baseline (Modified polynomial fitting)
bl <- baseline(data[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl <- data - bl@baseline
plot(spc.bl)
Spectrum <-spc.loess(spc.bl[,,350~1800], newx=350:1800)

# Difference Spectra(normal-cancer)
Spectrum_1 <- Spectrum[1]
Spectrum_2 <- Spectrum[2]
dif1 <- Spectrum_1 - Spectrum_2          #difference spectrum for EV5
################################################################################
files <-Sys.glob("~/Downloads/IPHT-Data/Man/Diff_Droplet/EV12/*.[sS][pP][cC]")
data  <- sapply (files, read.spc.Kaiser)
data  <- collapse (data)
data  <- data[,, c(350~1800)]
labels(data, "spc") <- expression ("Raman Intensity")
plotspc(data, col = matlab.dark.palette (nrow (data)),
  title.args= list(main = "EV5", 
  xlab = expression(paste("Wavenumber /", cm^-1, ""))))

# Baseline (Modified polynomial fitting)
bl <- baseline(data[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl <- data - bl@baseline
plot(spc.bl)
Spectrum <-spc.loess(spc.bl[,,350~1800], newx=350:1800)

# Difference Spectra(normal-cancer)
Spectrum_1 <- Spectrum[1]
Spectrum_2 <- Spectrum[2]
dif2 <- Spectrum_1 - Spectrum_2         #difference spectrum for EV12
################################################################################
files <-Sys.glob("~/Downloads/IPHT-Data/Man/Diff_Droplet/EV120/*.[sS][pP][cC]")
data  <- sapply (files, read.spc.Kaiser)
data  <- collapse (data)
data  <- data[,, c(350~1800)]
labels(data, "spc") <- expression ("Raman Intensity")
plotspc(data, col = matlab.dark.palette (nrow (data)),
   title.args= list(main = "EV5", 
   xlab = expression(paste("Wavenumber /", cm^-1, ""))))

# Baseline (Modified polynomial fitting)
bl <- baseline(data[[]], method = "modpolyfit", deg=6, tol=0.01, rep=100)
spc.bl <- data - bl@baseline
plot(spc.bl)
Spectrum <-spc.loess(spc.bl[,,350~1800], newx=350:1800)

# Difference Spectra(cancer-control)
Spectrum_1 <- Spectrum[1]
Spectrum_2 <- Spectrum[2]
dif3 <- Spectrum_1 - Spectrum_2         #difference spectrum for EV120
################################################################################

# Sum of difference bands
abline(v=c(382, 394, 713, 854, 1004, 1132, 1238, 1393, 1559,1589),col="gray", lty=2)
spc.identify(formatter = spc.label.wlonly,plot(dif1, col=1,zeroline="3"))

# combined difference spectra (EV5, EV12, EV120)
Diff = rbind(dif1,dif2,dif3)

#plot Difference spectrum
spc.identify(formatter = spc.label.wlonly,plot(Diff,wl.range = c (350 ~ 1800),col =c(2,3,4),
   #plot.args = list(ylim = c (-8000, 68000)),
   zeroline="3" ,title.args= list(xlab = expression(paste("Wavenumber /", cm^-1, "")))))

box()
legend("topleft", text.font=2, legend=c("Difference"), 
   col = 1, box.lty=1, lty=1,cex=1.0)


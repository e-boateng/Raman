# Data analysis for Raman & SERS characterization of extracellular vesicles

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

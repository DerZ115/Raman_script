start_time = format(Sys.time(), '%X') # get starttime of analysis

### USER SETTINGS ###

VERSION_NR <- 001                    # change for multiple analysis (only for documentation)

work_dir <- ("C:/Users/Daniel/Desktop/Raman_script/")     # Working directory (name of project folder)

data_folder <- "Data/"              # Name of your Data folder in your working directory
groups <- c("E","Z")                 # Filenames (unique part)
ending <- ".TXT"                     # Filetype
legend <- c("control","induced")     # Naming in publication

plot_original <- T                   # plot original spectra: yes = T | no = F
Legend_pos <- "topright"             # position of legend in plots "topleft" or "topright"
num_x <- 1                           # number of plots in one single window: 1 or 2

max_range <- "1600"                  # reduce maximum wavenumber to...
min_range <- "300"                   # reduce minimum wavenumber to...

remove_area <- F                     # remove Area: yes = T | no = F
area <- c("1250", "1450")            # exclude values from-to wavenumber

select_outlier <- T                  # select outlier during PCA: yes = T |no = F
PC <- 8                              # Number of components in variance plots
PCs <- 4                             # Number of PCs in overview plots
PC_PCA <- 6                          # numper of PCs calculated in PCA

plot_loadings <- T                   # show loading plots
loadings <- 3                        # number of loadings plotted simultaneously
x_PC <- 1                            # PC on x-axis of score plots
y_PC <- 2                            # PC on y-axis of score plots

### Load packages
Packages <- c("dplyr","IDPmisc","prospectr","dendextend","baseline",
              "pls","plotrix","knitr","ggplot2","gridExtra","ggpubr",
              "ChemoSpec", "matrixStats", "stringr")

for (p in 1:length(Packages)) {
  require(Packages[p], character.only = TRUE)
}

### Function for import
# work_dir: Path to files
# data_folder: folder with data
# groups: String vector of unique groupnames
# ending: Filetype
# Maximum 8 Groups possible 

Import.data <- function(work_dir, data_folder, groups, ending) {
  
  location <- paste(work_dir, data_folder, sep="")
  
  Files_full <- list.files(path = location, pattern = ending)
  Files_full <- str_replace(Files_full, ending, "")
  x <- c()
  for (group in groups) {
    x <- c(x,grep(pattern = group, x = Files_full))
  }
  Files <- paste(Files_full[x], ending, sep="")
  
  import <- function(data) {
    setwd(location)
    df <- try(read.csv(data, header = FALSE, 
                       sep = ",",              
                       dec = "." )) 
    setwd(work_dir)
    return(df)  
  }
  
  # read files
  Raw.list <- lapply(Files,import)
  #combine list of dataframes
  Raw.data <- do.call("cbind", Raw.list)
  # remove wavenumbers
  Spectra <- as.data.frame(t(Raw.data[,-c(which(colnames(Raw.data) == "V1"))]))
  # extract wavenumbers and set to column names
  Wavenumber <- as.numeric(Raw.data[,1])
  colnames(Spectra) <- Wavenumber
  # check for groupvector
  if(groups[1] != Files[1]){
    # create groupvector
    groups_v <- c(1:length(Files))
    for (i in 1:length(groups)){
      Pos <- c(grep(pattern = as.character(groups[i]), x = str_replace(Files, ending, ""))) # check for parts of filenames
      groups_v[Pos] <- groups[i] # fill in vector
    }
  }else{
    groups_v <- sub(pattern=ending, x=groups, replacement="")
  }
  
  #create a list of groups, wavenumbers and spectra
  OriginalData <- list(Wavenumber,Spectra,as.factor(groups_v))
  names(OriginalData) <- c("Wavenumber","Spectra","Groups")
  rownames(OriginalData$Spectra) <- 1:length(OriginalData$Groups)
  
  return(OriginalData)
}

### Import data
Data <- Import.data(work_dir = work_dir, data_folder = data_folder, groups = groups, ending = ending)

### Function for Spectra plots
# Liste = list with Wavenumber, Spectra, Groups
# Spektren = which Spectra should be plotted
# Wellenzahl = Vektor with wavenumbers
# area = which area should be plotted
# Code = vector for colours 
# Bereich = select wavenumbers

plot.spectra <- function(Liste,
                         Spektren,
                         Wellenzahl = Liste$Wavenumber,
                         area = c(max(Wellenzahl),min(Wellenzahl)), 
                         Code = Liste$Groups, 
                         Bereich = Sample,
                         main = main){
  
  if(Spektren == "Spectra_d"){          # check for derivatives
    Wellenzahl = Liste$Wavenumber_d     
    Bereich = Sample_d
  }
  
  matplot(x = Wellenzahl[Bereich],
          y = t(as.matrix(Liste[[Spektren]][,Bereich])),
          lty = 1, 
          type = "l", 
          col = Code,
          main = main,
          xlab = "Wavenumber [1/cm]",
          ylab = "Intensity [-]",
          font = 2, 
          font.lab = 2,  ###font = 2 ist Fett gedruckt
          lab = c(20,15,10), 
          xlim = area,
          ylim = c(min(Liste[[Spektren]][,Bereich]),
                   max(Liste[[Spektren]][,Bereich])), 
          bty = "l", 
          family = "sans", 
          xaxs = "i")
  grid(lwd = 0.8)
}

### Plot of the original spectra
if (plot_original == T) {
  
  plot.spectra(Liste = Data,
               Spektren = "Spectra",
               Bereich = c(1:length(Data$Wavenumber)),
               area = c(max(Data$Wavenumber),min(Data$Wavenumber)),
               main = paste("Original Spectra: ","max Wavenumber = ",
                            max(Data$Wavenumber),", min Wavenumber = ",
                            min(Data$Wavenumber), 
                            sep=""))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
}

# Plot with reduced range
plot.spectra(Liste = Data,
             Spektren = "Spectra", 
             Bereich = c(1:length(Data$Wavenumber)),
             main = "Reduce Wavenumber Range")

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

abline(v = as.numeric(max_range), lty = "dashed")
abline(v = as.numeric(max(Data$Wavenumber)), lty = "dashed")
abline(v = as.numeric(min(Data$Wavenumber)), lty = "dashed")
abline(v = as.numeric(min_range), lty = "dashed")

arrows(x0 = as.numeric(max_range), 
       y0 = 0, 
       x1 = as.numeric(max(Data$Wavenumber)), 
       y1 = 0, 
       col = "blue", 
       lwd = 2, 
       length = 0.1, 
       code = 3)
arrows(x0 = as.numeric(min(Data$Wavenumber)), 
       y0 = 0, 
       x1 = as.numeric(min_range), 
       y1 = 0, 
       col = "blue", 
       lwd = 2, 
       length = 0.1, 
       code = 3)

Data$Spectra_min_max <- subset(Data$Spectra, 
                                      select = -c(which(colnames(Data$Spectra)==max_range):
                                                    which(colnames(Data$Spectra)==max(Data$Wavenumber)),
                                                  which(colnames(Data$Spectra)==min(Data$Wavenumber)):
                                                    which(colnames(Data$Spectra)==min_range)))

Data$Wavenumber_min_max <- as.numeric(colnames(Data$Spectra_min_max))

plot.spectra(Liste = Data,
             Spektren = "Spectra_min_max",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             Wellenzahl = Data$Wavenumber_min_max, 
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)), 
             main = paste("Reduced Wavenumber Range from", max_range, "to", min_range))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### Plot with excluded area
if (remove_area == T) {
  plot.spectra(Liste = Data, 
               Spektren = "Spectra_min_max",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               Wellenzahl = Data$Wavenumber_min_max,
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = "Exclude Wavenumber Area")
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  abline(v = as.numeric(area[1]), 
         lty = "dashed")
  abline(v = as.numeric(area[2]), 
         lty = "dashed")
  
  arrows(x0 = as.numeric(area[1]), 
         y0 = 0, 
         x1 = as.numeric(area[2]), 
         y1 = 0, 
         col = "blue", 
         lwd = 2, 
         length = 0.1, 
         code = 3)
  
  Data$Spectra_area <- subset(Data$Spectra_min_max, 
                                     select = -c(which(colnames(Data$Spectra_min_max)==area[1]):
                                                   which(colnames(Data$Spectra_min_max)==area[2])))
  
  Data$Wavenumber_area <- as.numeric(colnames(Data$Spectra_area))
  
  plot.spectra(Liste = Data,
               Spektren = "Spectra_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               Wellenzahl = Data$Wavenumber_area, 
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)), 
               main = paste("Spectra with excluded Area from", area[2], "to", area[1]))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
}

### Data pre-treatment
# Data types
Type_rw <- c("Spectra_SNV_rw", "Spectra_det_rw", "Spectra_als_rw", "Spectra_fil_rw","Spectra_irl_rw", 
          "Spectra_med_rw", "Spectra_mod_rw", "Spectra_pea_rw", "Spectra_rbe_rw", "Spectra_rol_rw")

Type_a <- c("Spectra_SNV_area", "Spectra_det_area", "Spectra_als_area", "Spectra_fil_area","Spectra_irl_area", 
          "Spectra_med_area", "Spectra_mod_area", "Spectra_pea_area", "Spectra_rbe_area", "Spectra_rol_area")
# Data names for plots
Names <- c("SNV", "Detrend", "Asymmetric Least Squares", "FillPeaks", "Iterative Restricted Least Squares", 
           "Median window", "Modified polynomial fitting", "Simultaneous Peak Detection and Baseline Correction", 
           "Robust Baseline Estimation", "Rolling ball")

### Standard Normal Variate - REDUCED WAVENUMBER                                                                   #
par(mfrow = c(num_x,1))
Data$Spectra_SNV_rw <- standardNormalVariate(Data$Spectra_min_max)

# calculate mean of spectra
Spectra_SNV_mean_rw <- data.frame()
for (name in groups) {
  Spectra_SNV_mean_rw <- rbind(Spectra_SNV_mean_rw, t(colMeans(as.matrix(Data$Spectra_SNV_rw[Data$Groups==name,]))))
}
rownames(Spectra_SNV_mean_rw) <- groups
Data$Spectra_SNV_mean_rw <- Spectra_SNV_mean_rw

plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_SNV_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Standard Normal Variate, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_SNV_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

# calculate median of spectra
Spectra_SNV_median_rw <- data.frame()
for (name in groups) {
  Spectra_SNV_median_rw <- rbind(Spectra_SNV_median_rw, t(colMedians(as.matrix(Data$Spectra_SNV_rw[Data$Groups==name,]))))
}
rownames(Spectra_SNV_median_rw) <- groups
Data$Spectra_SNV_median_rw <- Spectra_SNV_median_rw

plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_SNV_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Standard Normal Variate, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_SNV_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### Standard Normal Variate - REDUCED WAVENUMBER/EXCLUDED AREA
if (remove_area == T) {
  
  Data$Spectra_SNV_area <- standardNormalVariate(Data$Spectra_area)
  
  # calculate mean of spectra
  Spectra_SNV_mean_area <- data.frame()
  for (name in groups) {
    Spectra_SNV_mean_area <- rbind(Spectra_SNV_mean_area, t(colMeans(as.matrix(Data$Spectra_SNV_area[Data$Groups==name,]))))
  }
  rownames(Spectra_SNV_mean_area) <- groups
  Data$Spectra_SNV_mean_area <- Spectra_SNV_mean_area
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_SNV_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Standard Normal Variate, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_SNV_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  # calculate median of spectra
  Spectra_SNV_median_area <- data.frame()
  for (name in groups) {
    Spectra_SNV_median_area <- rbind(Spectra_SNV_median_area, t(colMedians(as.matrix(Data$Spectra_SNV_area[Data$Groups==name,]))))
  }
  rownames(Spectra_SNV_median_area) <- groups
  Data$Spectra_SNV_median_area <- Spectra_SNV_median_area
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_SNV_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Standard Normal Variate, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_SNV_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  par(mfrow = c(1,1))
}

### Detrend - REDUCED WAVENUMBER
par(mfrow = c(num_x,1))

Data$Spectra_det_rw <- detrend(Data$Spectra_min_max, wav = Data$Wavenumber_min_max)

# calculate mean of spectra
Spectra_det_mean_rw <- data.frame()
for (name in groups) {
  Spectra_det_mean_rw <- rbind(Spectra_det_mean_rw, t(colMeans(as.matrix(Data$Spectra_det_rw[Data$Groups==name,]))))
}
rownames(Spectra_det_mean_rw) <- groups
Data$Spectra_det_mean_rw <- Spectra_det_mean_rw

plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_det_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Detrend, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_det_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

# calculate median of spectra
Spectra_det_median_rw <- data.frame()
for (name in groups) {
  Spectra_det_median_rw <- rbind(Spectra_det_median_rw, t(colMeans(as.matrix(Data$Spectra_det_rw[Data$Groups==name,]))))
}
rownames(Spectra_det_median_rw) <- groups
Data$Spectra_det_median_rw <- Spectra_det_median_rw

plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_det_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Detrend, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_det_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### Detrend - REDUCED WAVENUMBER/EXCLUDED AREA
if (remove_area == T) {
  

  Data$Spectra_det_area <- detrend(Data$Spectra_area, wav = Data$Wavenumber_area)
  
  # calculate mean of spectra
  Spectra_det_mean_area <- data.frame()
  for (name in groups) {
    Spectra_det_mean_area <- rbind(Spectra_det_mean_area, t(colMeans(as.matrix(Data$Spectra_det_area[Data$Groups==name,]))))
  }
  rownames(Spectra_det_mean_area) <- groups
  Data$Spectra_det_mean_area <- Spectra_det_mean_area
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_det_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Detrend, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_det_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  # calculate median of spectra
  Spectra_det_median_area <- data.frame()
  for (name in groups) {
    Spectra_det_median_area <- rbind(Spectra_det_median_area, t(colMeans(as.matrix(Data$Spectra_det_area[Data$Groups==name,]))))
  }
  rownames(Spectra_det_median_area) <- groups
  Data$Spectra_det_median_area <- Spectra_det_median_area
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_det_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Detrend, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_det_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  par(mfrow = c(1,1))
}

### Asymmetric Least Squares - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_als_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                   method = "als",  
                   lambda = 5, 
                   p = 0.05, 
                   maxit = 20)

# extract data from baseline object
Data$Spectra_als_orig_rw <- as.data.frame(getSpectra(bc_als_rw))
colnames(bc_als_rw@baseline) <- colnames(bc_als_rw@spectra)
Data$Spectra_als_bl_rw <- as.data.frame(getBaseline(bc_als_rw))
Data$Spectra_als_rw <- as.data.frame(getCorrected(bc_als_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_als_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_als_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_als_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_als_mean_bl_rw <- rbind(Spectra_als_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_als_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_als_mean_bl_rw) <- groups
Data$Spectra_als_mean_bl_rw <- Spectra_als_mean_bl_rw


### calculate median of baseline

Spectra_als_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_als_median_bl_rw <- rbind(Spectra_als_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_als_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_als_median_bl_rw) <- groups
Data$Spectra_als_median_bl_rw <- Spectra_als_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_als_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_als_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_rw)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_als_mean_rw <- data.frame()
for (name in groups) {
  Spectra_als_mean_rw <- rbind(Spectra_als_mean_rw, t(colMeans(as.matrix(Data$Spectra_als_rw[Data$Groups==name,]))))
}
rownames(Spectra_als_mean_rw) <- groups
Data$Spectra_als_mean_rw <- Spectra_als_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_als_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_als_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_als_median_rw <- data.frame()
for (name in groups) {
  Spectra_als_median_rw <- rbind(Spectra_als_median_rw, t(colMedians(as.matrix(Data$Spectra_als_rw[Data$Groups==name,]))))
}
rownames(Spectra_als_median_rw) <- groups
Data$Spectra_als_median_rw <- Spectra_als_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_als_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_als_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### Asymmetric Least Squares - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {

  # calculate baseline with "baseline" package
  bc_als_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                        method = "als",  
                        lambda = 7, 
                        p = 0.05, 
                        maxit = 20)
  
  # extract data from baseline object
  Data$Spectra_als_orig_area <- as.data.frame(getSpectra(bc_als_area))
  colnames(bc_als_area@baseline) <- colnames(bc_als_area@spectra)
  Data$Spectra_als_bl_area <- as.data.frame(getBaseline(bc_als_area))
  Data$Spectra_als_area <- as.data.frame(getCorrected(bc_als_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_als_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_als_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_als_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_als_mean_bl_area <- rbind(Spectra_als_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_als_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_als_mean_bl_area) <- groups
  Data$Spectra_als_mean_bl_area <- Spectra_als_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_als_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_als_median_bl_area <- rbind(Spectra_als_median_bl_area, t(colMedians(as.matrix(Data$Spectra_als_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_als_median_bl_area) <- groups
  Data$Spectra_als_median_bl_area <- Spectra_als_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_als_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Asymmetric Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_area)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_als_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Asymmetric Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_area)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_als_mean_area <- data.frame()
  for (name in groups) {
    Spectra_als_mean_area <- rbind(Spectra_als_mean_area, t(colMeans(as.matrix(Data$Spectra_als_area[Data$Groups==name,]))))
  }
  rownames(Spectra_als_mean_area) <- groups
  Data$Spectra_als_mean_area <- Spectra_als_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_als_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Asymmetric Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_als_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_als_median_area <- data.frame()
  for (name in groups) {
    Spectra_als_median_area <- rbind(Spectra_als_median_area, t(colMedians(as.matrix(Data$Spectra_als_area[Data$Groups==name,]))))
  }
  rownames(Spectra_als_median_area) <- groups
  Data$Spectra_als_median_area <- Spectra_als_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_als_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Asymmetric Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_als_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}

### FillPeaks - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_fil_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "fillPeaks",  
                      lambda=5, 
                      hwi=6, 
                      it=5, 
                      int=200)

# extract data from baseline object
Data$Spectra_fil_orig_rw <- as.data.frame(getSpectra(bc_fil_rw))
colnames(bc_fil_rw@baseline) <- colnames(bc_fil_rw@spectra)
Data$Spectra_fil_bl_rw <- as.data.frame(getBaseline(bc_fil_rw))
Data$Spectra_fil_rw <- as.data.frame(getCorrected(bc_fil_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_fil_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_fil_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_fil_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_fil_mean_bl_rw <- rbind(Spectra_fil_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_fil_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_fil_mean_bl_rw) <- groups
Data$Spectra_fil_mean_bl_rw <- Spectra_fil_mean_bl_rw


### calculate median of baseline

Spectra_fil_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_fil_median_bl_rw <- rbind(Spectra_fil_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_fil_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_fil_median_bl_rw) <- groups
Data$Spectra_fil_median_bl_rw <- Spectra_fil_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_fil_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_fil_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_fil_mean_rw <- data.frame()
for (name in groups) {
  Spectra_fil_mean_rw <- rbind(Spectra_fil_mean_rw, t(colMeans(as.matrix(Data$Spectra_fil_rw[Data$Groups==name,]))))
}
rownames(Spectra_fil_mean_rw) <- groups
Data$Spectra_fil_mean_rw <- Spectra_fil_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_fil_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_fil_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

# calculate median of corrected reduced spectra

Spectra_fil_median_rw <- data.frame()
for (name in groups) {
  Spectra_fil_median_rw <- rbind(Spectra_fil_median_rw, t(colMedians(as.matrix(Data$Spectra_fil_rw[Data$Groups==name,]))))
}
rownames(Spectra_fil_median_rw) <- groups
Data$Spectra_fil_median_rw <- Spectra_fil_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_fil_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_fil_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### FillPeaks - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {

  # calculate baseline with "baseline" package
  bc_fil_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "als",  
                          lambda = 7, 
                          p = 0.05, 
                          maxit = 20)
  
  # extract data from baseline object
  Data$Spectra_fil_orig_area <- as.data.frame(getSpectra(bc_fil_area))
  colnames(bc_fil_area@baseline) <- colnames(bc_fil_area@spectra)
  Data$Spectra_fil_bl_area <- as.data.frame(getBaseline(bc_fil_area))
  Data$Spectra_fil_area <- as.data.frame(getCorrected(bc_fil_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_fil_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_fil_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_fil_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_fil_mean_bl_area <- rbind(Spectra_fil_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_fil_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_mean_bl_area) <- groups
  Data$Spectra_fil_mean_bl_area <- Spectra_fil_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_fil_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_fil_median_bl_area <- rbind(Spectra_fil_median_bl_area, t(colMedians(as.matrix(Data$Spectra_fil_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_median_bl_area) <- groups
  Data$Spectra_fil_median_bl_area <- Spectra_fil_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_fil_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_fil_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_fil_mean_area <- data.frame()
  for (name in groups) {
    Spectra_fil_mean_area <- rbind(Spectra_fil_mean_area, t(colMeans(as.matrix(Data$Spectra_fil_area[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_mean_area) <- groups
  Data$Spectra_fil_mean_area <- Spectra_fil_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_fil_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_fil_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_fil_median_area <- data.frame()
  for (name in groups) {
    Spectra_fil_median_area <- rbind(Spectra_fil_median_area, t(colMedians(as.matrix(Data$Spectra_fil_area[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_median_area) <- groups
  Data$Spectra_fil_median_area <- Spectra_fil_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_fil_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_fil_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}


### Iterative Restricted Least Squares - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN                                #
# calculate baseline with "baseline" package
bc_irl_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "irls", 
                      lambda1 = 5, 
                      lambda2 = 9, 
                      maxit = 200, 
                      wi = 0.05)

# extract data from baseline object
Data$Spectra_irl_orig_rw <- as.data.frame(getSpectra(bc_irl_rw))
colnames(bc_irl_rw@baseline) <- colnames(bc_irl_rw@spectra)
Data$Spectra_irl_bl_rw <- as.data.frame(getBaseline(bc_irl_rw))
Data$Spectra_irl_rw <- as.data.frame(getCorrected(bc_irl_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_irl_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_irl_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_irl_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_irl_mean_bl_rw <- rbind(Spectra_irl_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_irl_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_irl_mean_bl_rw) <- groups
Data$Spectra_irl_mean_bl_rw <- Spectra_irl_mean_bl_rw


### calculate median of baseline

Spectra_irl_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_irl_median_bl_rw <- rbind(Spectra_irl_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_irl_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_irl_median_bl_rw) <- groups
Data$Spectra_irl_median_bl_rw <- Spectra_irl_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_irl_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_irl_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_irl_mean_rw <- data.frame()
for (name in groups) {
  Spectra_irl_mean_rw <- rbind(Spectra_irl_mean_rw, t(colMeans(as.matrix(Data$Spectra_irl_rw[Data$Groups==name,]))))
}
rownames(Spectra_irl_mean_rw) <- groups
Data$Spectra_irl_mean_rw <- Spectra_irl_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_irl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_irl_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_irl_median_rw <- data.frame()
for (name in groups) {
  Spectra_irl_median_rw <- rbind(Spectra_irl_median_rw, t(colMedians(as.matrix(Data$Spectra_irl_rw[Data$Groups==name,]))))
}
rownames(Spectra_irl_median_rw) <- groups
Data$Spectra_irl_median_rw <- Spectra_irl_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_irl_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_irl_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### Iterative Restricted Least Squares - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {

  # calculate baseline with "baseline" package
  bc_irl_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "irls", 
                          lambda1 = 5, 
                          lambda2 = 9, 
                          maxit = 200, 
                          wi = 0.05)
  
  # extract data from baseline object
  Data$Spectra_irl_orig_area <- as.data.frame(getSpectra(bc_irl_area))
  colnames(bc_irl_area@baseline) <- colnames(bc_irl_area@spectra)
  Data$Spectra_irl_bl_area <- as.data.frame(getBaseline(bc_irl_area))
  Data$Spectra_irl_area <- as.data.frame(getCorrected(bc_irl_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_irl_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_irl_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_irl_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_irl_mean_bl_area <- rbind(Spectra_irl_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_irl_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_mean_bl_area) <- groups
  Data$Spectra_irl_mean_bl_area <- Spectra_irl_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_irl_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_irl_median_bl_area <- rbind(Spectra_irl_median_bl_area, t(colMedians(as.matrix(Data$Spectra_irl_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_median_bl_area) <- groups
  Data$Spectra_irl_median_bl_area <- Spectra_irl_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_irl_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_irl_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_irl_mean_area <- data.frame()
  for (name in groups) {
    Spectra_irl_mean_area <- rbind(Spectra_irl_mean_area, t(colMeans(as.matrix(Data$Spectra_irl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_mean_area) <- groups
  Data$Spectra_irl_mean_area <- Spectra_irl_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_irl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Iterative Restricted Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_irl_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_irl_median_area <- data.frame()
  for (name in groups) {
    Spectra_irl_median_area <- rbind(Spectra_irl_median_area, t(colMedians(as.matrix(Data$Spectra_irl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_median_area) <- groups
  Data$Spectra_irl_median_area <- Spectra_irl_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_irl_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Iterative Restricted Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_irl_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}



### Median window - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_med_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "medianWindow", 
                      hwm = 300)

# extract data from baseline object
Data$Spectra_med_orig_rw <- as.data.frame(getSpectra(bc_med_rw))
colnames(bc_med_rw@baseline) <- colnames(bc_med_rw@spectra)
Data$Spectra_med_bl_rw <- as.data.frame(getBaseline(bc_med_rw))
Data$Spectra_med_rw <- as.data.frame(getCorrected(bc_med_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_med_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_med_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_med_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_med_mean_bl_rw <- rbind(Spectra_med_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_med_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_med_mean_bl_rw) <- groups
Data$Spectra_med_mean_bl_rw <- Spectra_med_mean_bl_rw


### calculate median of baseline

Spectra_med_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_med_median_bl_rw <- rbind(Spectra_med_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_med_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_med_median_bl_rw) <- groups
Data$Spectra_med_median_bl_rw <- Spectra_med_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_med_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Median Window, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_med_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Median Window, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_med_mean_rw <- data.frame()
for (name in groups) {
  Spectra_med_mean_rw <- rbind(Spectra_med_mean_rw, t(colMeans(as.matrix(Data$Spectra_med_rw[Data$Groups==name,]))))
}
rownames(Spectra_med_mean_rw) <- groups
Data$Spectra_med_mean_rw <- Spectra_med_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_med_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Median Window, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_med_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_med_median_rw <- data.frame()
for (name in groups) {
  Spectra_med_median_rw <- rbind(Spectra_med_median_rw, t(colMedians(as.matrix(Data$Spectra_med_rw[Data$Groups==name,]))))
}
rownames(Spectra_med_median_rw) <- groups
Data$Spectra_med_median_rw <- Spectra_med_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_med_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Median Window, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_med_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### Median Window - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {
  
  # calculate baseline with "baseline" package
  bc_med_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "medianWindow", 
                          hwm = 300)
    
  # extract data from baseline object
  Data$Spectra_med_orig_area <- as.data.frame(getSpectra(bc_med_area))
  colnames(bc_med_area@baseline) <- colnames(bc_med_area@spectra)
  Data$Spectra_med_bl_area <- as.data.frame(getBaseline(bc_med_area))
  Data$Spectra_med_area <- as.data.frame(getCorrected(bc_med_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_med_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_med_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_med_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_med_mean_bl_area <- rbind(Spectra_med_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_med_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_med_mean_bl_area) <- groups
  Data$Spectra_med_mean_bl_area <- Spectra_med_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_med_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_med_median_bl_area <- rbind(Spectra_med_median_bl_area, t(colMedians(as.matrix(Data$Spectra_med_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_med_median_bl_area) <- groups
  Data$Spectra_med_median_bl_area <- Spectra_med_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_med_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Median Window, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_med_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Median Window, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_med_mean_area <- data.frame()
  for (name in groups) {
    Spectra_med_mean_area <- rbind(Spectra_med_mean_area, t(colMeans(as.matrix(Data$Spectra_med_area[Data$Groups==name,]))))
  }
  rownames(Spectra_med_mean_area) <- groups
  Data$Spectra_med_mean_area <- Spectra_med_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_med_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Median Window, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_med_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_med_median_area <- data.frame()
  for (name in groups) {
    Spectra_med_median_area <- rbind(Spectra_med_median_area, t(colMedians(as.matrix(Data$Spectra_med_area[Data$Groups==name,]))))
  }
  rownames(Spectra_med_median_area) <- groups
  Data$Spectra_med_median_area <- Spectra_med_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_med_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Median Window, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_med_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}

### Modified Polynomial Fitting - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_mod_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "modpolyfit", 
                      deg = 8)

# extract data from baseline object
Data$Spectra_mod_orig_rw <- as.data.frame(getSpectra(bc_mod_rw))
colnames(bc_mod_rw@baseline) <- colnames(bc_mod_rw@spectra)
Data$Spectra_mod_bl_rw <- as.data.frame(getBaseline(bc_mod_rw))
Data$Spectra_mod_rw <- as.data.frame(getCorrected(bc_mod_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_mod_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_mod_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_mod_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_mod_mean_bl_rw <- rbind(Spectra_mod_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_mod_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_mod_mean_bl_rw) <- groups
Data$Spectra_mod_mean_bl_rw <- Spectra_mod_mean_bl_rw


### calculate median of baseline

Spectra_mod_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_mod_median_bl_rw <- rbind(Spectra_mod_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_mod_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_mod_median_bl_rw) <- groups
Data$Spectra_mod_median_bl_rw <- Spectra_mod_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_mod_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_mod_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_mod_mean_rw <- data.frame()
for (name in groups) {
  Spectra_mod_mean_rw <- rbind(Spectra_mod_mean_rw, t(colMeans(as.matrix(Data$Spectra_mod_rw[Data$Groups==name,]))))
}
rownames(Spectra_mod_mean_rw) <- groups
Data$Spectra_mod_mean_rw <- Spectra_mod_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_mod_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_mod_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_mod_median_rw <- data.frame()
for (name in groups) {
  Spectra_mod_median_rw <- rbind(Spectra_mod_median_rw, t(colMedians(as.matrix(Data$Spectra_mod_rw[Data$Groups==name,]))))
}
rownames(Spectra_mod_median_rw) <- groups
Data$Spectra_mod_median_rw <- Spectra_mod_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_mod_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_mod_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### Modified Polynomial Fitting - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {
  
  # calculate baseline with "baseline" package
  bc_mod_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "modpolyfit", 
                          deg = 8)
  
  # extract data from baseline object
  Data$Spectra_mod_orig_area <- as.data.frame(getSpectra(bc_mod_area))
  colnames(bc_mod_area@baseline) <- colnames(bc_mod_area@spectra)
  Data$Spectra_mod_bl_area <- as.data.frame(getBaseline(bc_mod_area))
  Data$Spectra_mod_area <- as.data.frame(getCorrected(bc_mod_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_mod_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_mod_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_mod_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_mod_mean_bl_area <- rbind(Spectra_mod_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_mod_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_mean_bl_area) <- groups
  Data$Spectra_mod_mean_bl_area <- Spectra_mod_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_mod_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_mod_median_bl_area <- rbind(Spectra_mod_median_bl_area, t(colMedians(as.matrix(Data$Spectra_mod_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_median_bl_area) <- groups
  Data$Spectra_mod_median_bl_area <- Spectra_mod_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_mod_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Modified Polynomial Fitting, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_mod_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Modified Polynomial Fitting, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_mod_mean_area <- data.frame()
  for (name in groups) {
    Spectra_mod_mean_area <- rbind(Spectra_mod_mean_area, t(colMeans(as.matrix(Data$Spectra_mod_area[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_mean_area) <- groups
  Data$Spectra_mod_mean_area <- Spectra_mod_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_mod_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Modified Polynomial Fitting, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_mod_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_mod_median_area <- data.frame()
  for (name in groups) {
    Spectra_mod_median_area <- rbind(Spectra_mod_median_area, t(colMedians(as.matrix(Data$Spectra_mod_area[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_median_area) <- groups
  Data$Spectra_mod_median_area <- Spectra_mod_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_mod_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Modified Polynomial Fitting, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_mod_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}


### Simultaneous Peak Detection and Baseline Correction - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
options(warn=-1) # suppress errors
# calculate baseline with "baseline" package
bc_pea_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "peakDetection", 
                      left=30,
                      right=30,
                      lwin=50, 
                      rwin=50)

# extract data from baseline object
Data$Spectra_pea_orig_rw <- as.data.frame(getSpectra(bc_pea_rw))
colnames(bc_pea_rw@baseline) <- colnames(bc_pea_rw@spectra)
Data$Spectra_pea_bl_rw <- as.data.frame(getBaseline(bc_pea_rw))
Data$Spectra_pea_rw <- as.data.frame(getCorrected(bc_pea_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_pea_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_pea_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_pea_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_pea_mean_bl_rw <- rbind(Spectra_pea_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_pea_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_pea_mean_bl_rw) <- groups
Data$Spectra_pea_mean_bl_rw <- Spectra_pea_mean_bl_rw


### calculate median of baseline

Spectra_pea_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_pea_median_bl_rw <- rbind(Spectra_pea_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_pea_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_pea_median_bl_rw) <- groups
Data$Spectra_pea_median_bl_rw <- Spectra_pea_median_bl_rw



### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_pea_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_pea_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_pea_mean_rw <- data.frame()
for (name in groups) {
  Spectra_pea_mean_rw <- rbind(Spectra_pea_mean_rw, t(colMeans(as.matrix(Data$Spectra_pea_rw[Data$Groups==name,]))))
}
rownames(Spectra_pea_mean_rw) <- groups
Data$Spectra_pea_mean_rw <- Spectra_pea_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_pea_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_pea_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_pea_median_rw <- data.frame()
for (name in groups) {
  Spectra_pea_median_rw <- rbind(Spectra_pea_median_rw, t(colMedians(as.matrix(Data$Spectra_pea_rw[Data$Groups==name,]))))
}
rownames(Spectra_pea_median_rw) <- groups
Data$Spectra_pea_median_rw <- Spectra_pea_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_pea_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_pea_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))
options(warn=0)

### Simultaneous Peak Detection and Baseline Correction - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
options(warn=-1)
if (remove_area == T) {
  
  # calculate baseline with "baseline" package
  bc_pea_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "peakDetection")#, 
                          #left=300, 
                          #right=300, 
                          #lwin=50, 
                          #rwin=50)
  
  # extract data from baseline object
  Data$Spectra_pea_orig_area <- as.data.frame(getSpectra(bc_pea_area))
  colnames(bc_pea_area@baseline) <- colnames(bc_pea_area@spectra)
  Data$Spectra_pea_bl_area <- as.data.frame(getBaseline(bc_pea_area))
  Data$Spectra_pea_area <- as.data.frame(getCorrected(bc_pea_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_pea_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_pea_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_pea_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_pea_mean_bl_area <- rbind(Spectra_pea_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_pea_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_mean_bl_area) <- groups
  Data$Spectra_pea_mean_bl_area <- Spectra_pea_mean_bl_area
  
  
  ### calculate median of baseline
  
  Spectra_pea_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_pea_median_bl_area <- rbind(Spectra_pea_median_bl_area, t(colMedians(as.matrix(Data$Spectra_pea_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_median_bl_area) <- groups
  Data$Spectra_pea_median_bl_area <- Spectra_pea_median_bl_area
  
  
  
  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_pea_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_pea_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_pea_mean_area <- data.frame()
  for (name in groups) {
    Spectra_pea_mean_area <- rbind(Spectra_pea_mean_area, t(colMeans(as.matrix(Data$Spectra_pea_area[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_mean_area) <- groups
  Data$Spectra_pea_mean_area <- Spectra_pea_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_pea_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Simultaneous Peak Detection and Baseline Correction, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_pea_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_pea_median_area <- data.frame()
  for (name in groups) {
    Spectra_pea_median_area <- rbind(Spectra_pea_median_area, t(colMedians(as.matrix(Data$Spectra_pea_area[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_median_area) <- groups
  Data$Spectra_pea_median_area <- Spectra_pea_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_pea_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Simultaneous Peak Detection and Baseline Correction, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_pea_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}
options(warn=0)

### Robust Baseline Estimation - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_rbe_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "rfbaseline", 
                      span=10, 
                      NoXP=100)

# extract data from baseline object
Data$Spectra_rbe_orig_rw <- as.data.frame(getSpectra(bc_rbe_rw))
colnames(bc_rbe_rw@baseline) <- colnames(bc_rbe_rw@spectra)
Data$Spectra_rbe_bl_rw <- as.data.frame(getBaseline(bc_rbe_rw))
Data$Spectra_rbe_rw <- as.data.frame(getCorrected(bc_rbe_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_rbe_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_rbe_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_rbe_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_rbe_mean_bl_rw <- rbind(Spectra_rbe_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_rbe_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_rbe_mean_bl_rw) <- groups
Data$Spectra_rbe_mean_bl_rw <- Spectra_rbe_mean_bl_rw


### calculate median of baseline

Spectra_rbe_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_rbe_median_bl_rw <- rbind(Spectra_rbe_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_rbe_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_rbe_median_bl_rw) <- groups
Data$Spectra_rbe_median_bl_rw <- Spectra_rbe_median_bl_rw

### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_rbe_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_rbe_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### calculate mean of corrected reduced spectra

Spectra_rbe_mean_rw <- data.frame()
for (name in groups) {
  Spectra_rbe_mean_rw <- rbind(Spectra_rbe_mean_rw, t(colMeans(as.matrix(Data$Spectra_rbe_rw[Data$Groups==name,]))))
}
rownames(Spectra_rbe_mean_rw) <- groups
Data$Spectra_rbe_mean_rw <- Spectra_rbe_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_rbe_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_rbe_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_rbe_median_rw <- data.frame()
for (name in groups) {
  Spectra_rbe_median_rw <- rbind(Spectra_rbe_median_rw, t(colMedians(as.matrix(Data$Spectra_rbe_rw[Data$Groups==name,]))))
}
rownames(Spectra_rbe_median_rw) <- groups
Data$Spectra_rbe_median_rw <- Spectra_rbe_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_rbe_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_rbe_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

### Robust Baseline Estimation - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {
  
  # calculate baseline with "baseline" package
  bc_rbe_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "rfbaseline", 
                          span=NULL, 
                          NoXP=1000)
  
  # extract data from baseline object
  Data$Spectra_rbe_orig_area <- as.data.frame(getSpectra(bc_rbe_area))
  colnames(bc_rbe_area@baseline) <- colnames(bc_rbe_area@spectra)
  Data$Spectra_rbe_bl_area <- as.data.frame(getBaseline(bc_rbe_area))
  Data$Spectra_rbe_area <- as.data.frame(getCorrected(bc_rbe_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_rbe_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
  
  
  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_rbe_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area
  
  
  ### calculate mean of baseline
  
  Spectra_rbe_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_rbe_mean_bl_area <- rbind(Spectra_rbe_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_rbe_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_mean_bl_area) <- groups
  Data$Spectra_rbe_mean_bl_area <- Spectra_rbe_mean_bl_area

  ### calculate median of baseline
  
  Spectra_rbe_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_rbe_median_bl_area <- rbind(Spectra_rbe_median_bl_area, t(colMedians(as.matrix(Data$Spectra_rbe_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_median_bl_area) <- groups
  Data$Spectra_rbe_median_bl_area <- Spectra_rbe_median_bl_area

  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_rbe_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Robust Baseline Estimation, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_rbe_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Robust Baseline Estimation, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")

  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_rbe_mean_area <- data.frame()
  for (name in groups) {
    Spectra_rbe_mean_area <- rbind(Spectra_rbe_mean_area, t(colMeans(as.matrix(Data$Spectra_rbe_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_mean_area) <- groups
  Data$Spectra_rbe_mean_area <- Spectra_rbe_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_rbe_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Robust Baseline Estimation, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_rbe_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_rbe_median_area <- data.frame()
  for (name in groups) {
    Spectra_rbe_median_area <- rbind(Spectra_rbe_median_area, t(colMedians(as.matrix(Data$Spectra_rbe_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_median_area) <- groups
  Data$Spectra_rbe_median_area <- Spectra_rbe_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_rbe_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Robust Baseline Estimation, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_rbe_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}

### Rolling ball - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN
# calculate baseline with "baseline" package
bc_rol_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                      method = "rollingBall", 
                      wm=20, 
                      ws=20)

# extract data from baseline object
Data$Spectra_rol_orig_rw <- as.data.frame(getSpectra(bc_rol_rw))
colnames(bc_rol_rw@baseline) <- colnames(bc_rol_rw@spectra)
Data$Spectra_rol_bl_rw <- as.data.frame(getBaseline(bc_rol_rw))
Data$Spectra_rol_rw <- as.data.frame(getCorrected(bc_rol_rw))

par(mfrow = c(num_x,1))

### calculate mean of reduced spectra

Spectra_orig_mean_rw <- data.frame()
for (name in groups) {
  Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_rol_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_mean_rw) <- groups
Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw

### calculate median of reduced spectra

Spectra_orig_median_rw <- data.frame()
for (name in groups) {
  Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_rol_orig_rw[Data$Groups==name,]))))
}
rownames(Spectra_orig_median_rw) <- groups
Data$Spectra_orig_median_rw <- Spectra_orig_median_rw

### calculate mean of baseline

Spectra_rol_mean_bl_rw <- data.frame()
for (name in groups) {
  Spectra_rol_mean_bl_rw <- rbind(Spectra_rol_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_rol_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_rol_mean_bl_rw) <- groups
Data$Spectra_rol_mean_bl_rw <- Spectra_rol_mean_bl_rw


### calculate median of baseline

Spectra_rol_median_bl_rw <- data.frame()
for (name in groups) {
  Spectra_rol_median_bl_rw <- rbind(Spectra_rol_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_rol_bl_rw[Data$Groups==name,]))))
}
rownames(Spectra_rol_median_bl_rw) <- groups
Data$Spectra_rol_median_bl_rw <- Spectra_rol_median_bl_rw

### combine mean of reduced spectra and baseline
Data$Spectra_orig_and_bl_mean_rw <- data.frame()
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_rol_mean_bl_rw)

# plot mean of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Mean of original spectra and calculated baseline with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### combine median of reduced spectra and baseline
Data$Spectra_orig_and_bl_med_rw <- data.frame()
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_rol_median_bl_rw)

# plot median of reduced spectra and baseline
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_orig_and_bl_med_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Median of original spectra and calculated baseline with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(c(groups,groups)))

legend(Legend_pos,
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate mean of corrected reduced spectra

Spectra_rol_mean_rw <- data.frame()
for (name in groups) {
  Spectra_rol_mean_rw <- rbind(Spectra_rol_mean_rw, t(colMeans(as.matrix(Data$Spectra_rol_rw[Data$Groups==name,]))))
}
rownames(Spectra_rol_mean_rw) <- groups
Data$Spectra_rol_mean_rw <- Spectra_rol_mean_rw

# plot mean of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_rol_mean_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (mean) corrected with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_rol_mean_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

### calculate median of corrected reduced spectra

Spectra_rol_median_rw <- data.frame()
for (name in groups) {
  Spectra_rol_median_rw <- rbind(Spectra_rol_median_rw, t(colMedians(as.matrix(Data$Spectra_rol_rw[Data$Groups==name,]))))
}
rownames(Spectra_rol_median_rw) <- groups
Data$Spectra_rol_median_rw <- Spectra_rol_median_rw

# plot median of corrected reduced spectra
plot.spectra(Liste = Data,
             Wellenzahl = Data$Wavenumber_min_max,
             Spektren = "Spectra_rol_median_rw",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
             main = paste("Spectra (median) corrected with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
             Code = as.factor(rownames(Data$Spectra_rol_median_rw)))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

par(mfrow = c(1,1))

# Rolling ball - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
if (remove_area == T) {
  
  # calculate baseline with "baseline" package
  bc_rol_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                          method = "rollingBall", 
                          wm=100, 
                          ws=100)
  
  # extract data from baseline object
  Data$Spectra_rol_orig_area <- as.data.frame(getSpectra(bc_rol_area))
  colnames(bc_rol_area@baseline) <- colnames(bc_rol_area@spectra)
  Data$Spectra_rol_bl_area <- as.data.frame(getBaseline(bc_rol_area))
  Data$Spectra_rol_area <- as.data.frame(getCorrected(bc_rol_area))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced/excluded area spectra
  
  Spectra_orig_mean_area <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_rol_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_area) <- groups
  Data$Spectra_orig_mean_area <- Spectra_orig_mean_area

  ### calculate median of reduced/excluded area spectra
  
  Spectra_orig_median_area <- data.frame()
  for (name in groups) {
    Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_rol_orig_area[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_area) <- groups
  Data$Spectra_orig_median_area <- Spectra_orig_median_area

  ### calculate mean of baseline
  
  Spectra_rol_mean_bl_area <- data.frame()
  for (name in groups) {
    Spectra_rol_mean_bl_area <- rbind(Spectra_rol_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_rol_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_mean_bl_area) <- groups
  Data$Spectra_rol_mean_bl_area <- Spectra_rol_mean_bl_area

  ### calculate median of baseline
  
  Spectra_rol_median_bl_area <- data.frame()
  for (name in groups) {
    Spectra_rol_median_bl_area <- rbind(Spectra_rol_median_bl_area, t(colMedians(as.matrix(Data$Spectra_rol_bl_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_median_bl_area) <- groups
  Data$Spectra_rol_median_bl_area <- Spectra_rol_median_bl_area

  ### combine mean of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_mean_area <- data.frame()
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
  Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_rol_mean_bl_area)
  
  # plot mean of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Mean of original spectra and calculated baseline with Rolling ball, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced/excluded area spectra and baseline
  Data$Spectra_orig_and_bl_med_area <- data.frame()
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
  Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_rol_median_bl_area)
  
  # plot median of reduced/excluded area spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_orig_and_bl_med_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Median of original spectra and calculated baseline with Rolling ball, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")

  ### calculate mean of corrected reduced/excluded area spectra
  
  Spectra_rol_mean_area <- data.frame()
  for (name in groups) {
    Spectra_rol_mean_area <- rbind(Spectra_rol_mean_area, t(colMeans(as.matrix(Data$Spectra_rol_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_mean_area) <- groups
  Data$Spectra_rol_mean_area <- Spectra_rol_mean_area
  
  # plot mean of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_rol_mean_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (mean) corrected with Rolling ball, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_rol_mean_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced/excluded area spectra
  
  Spectra_rol_median_area <- data.frame()
  for (name in groups) {
    Spectra_rol_median_area <- rbind(Spectra_rol_median_area, t(colMedians(as.matrix(Data$Spectra_rol_area[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_median_area) <- groups
  Data$Spectra_rol_median_area <- Spectra_rol_median_area
  
  # plot median of corrected reduced/excluded area spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_area,
               Spektren = "Spectra_rol_median_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
               main = paste("Spectra (median) corrected with Rolling ball, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
               Code = as.factor(rownames(Data$Spectra_rol_median_area)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
}

# PCA - reduced wavenumber 
par(mfrow = c(1,1))

Type = Type_rw
Sample <- c(1:length(Data$Wavenumber_min_max))

for (Stats in 1:length(Type)) {
  print("---------")
  print(Names[Stats])
  PCA <- prcomp(Data[[Type[Stats]]][!is.na(Data$Groups),Sample],
                center = TRUE,
                rank. = PC_PCA, 
                scale. = TRUE )
  
  # remove "NAs" for correct colours in the plot
  PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
  
  print(summary(PCA))
  
  Data[[paste("Var_", Type[Stats])]] <- PCA$sdev^2/sum(PCA$sdev^2)
  
  plot(x = 1:PC, 
       y = cumsum(Data[[paste("Var_", Type[Stats])]][1:PC]*100),
       type = "b", 
       main = paste("Method: ",Names[Stats], sep=""),
       xlab = "Number of Principal Components",
       ylab = "Explained Variance [%]",
       col = "dimgrey", 
       pch = 21, 
       bg = "darkgrey",
       font = 2, 
       font.lab = 2,
       lab = c(10,20,20),
       xaxt = "n")
  axis(1,xaxp=c(1,PC,PC-1))
  
  grid(lwd = 0.8)
  
  pairs(PCA$x[,1:PCs],
        col = Data$Groups,
        main = paste("Method: ",Names[Stats], sep=""),
        pch = 19)
  
  par(xpd = T) 
  
  legend("bottomright", 
         legend = legend, 
         pch = 19, 
         col = unique(Data$Groups), 
         inset = 0.05,
         cex = 0.7,
         horiz = T, 
         bty = "n")
  
  par(xpd = F)
  
  plot(x = PCA$x[,x_PC], 
       y = PCA$x[,y_PC],
       xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats])]][x_PC]*100,1),"%]", sep = " "),
       ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats])]][y_PC]*100,1),"%]", sep = " "),
       main = paste("Method: ",Names[Stats], sep=""),
       col = Data$Groups[!is.na(Data$Groups)], 
       pch = 19, 
       font.lab = 2)
  
  grid(lwd = 0.8)
  
  abline(h = 0, v = 0)
  
  legend("bottomleft", 
         legend = legend, 
         pch = 20, 
         col = unique(Data$Groups[!is.na(Data$Groups)]), 
         inset = 0.01,
         bty = "n")
  
  # remove outlier
  
  Data[paste(Type[Stats],"_outl",sep="")] <- Data[paste(Type[Stats])]
  Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]] <- Data$Groups
  
  if(select_outlier == T){
    print(paste(Stats, "Select outlier: ESC to exit"))
    outlier = identify( x = PCA$x[,x_PC], y = PCA$x[,y_PC])
    print(outlier)
    Data[[paste(Type[Stats],"_outl",sep="")]][outlier,] = NA
    Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][outlier] = NA
    #print(Data$Groups)
    #print(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])
  }
  
  if (plot_loadings == T) {
    
    DataLoadings <- data.frame(WaveNum=Data$Wavenumber_min_max)
    
    for (loadingNR in (1:loadings)) {
      #print(paste("PC",loadingNR,sep=""))
      DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
      colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
    }
    
    DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
    
    matplot(x = DataLoadingsRev[1],
            y = DataLoadingsRev[-1],
            type = "l", 
            main = paste("Loadings: ",Names[Stats], sep=""), 
            ylab = "",
            xlab = "Wavenumber [1/cm]",
            xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)),
            xaxs = "i")
    abline(0,0, col="black", lwd=1.5)
    legend("topleft", 
           legend = colnames(DataLoadings)[2:(loadings+1)], 
           pch = 20, 
           col = c(1:loadings), 
           inset = 0.01,
           bty = "n")
    grid(lwd = 0.8)
    
  }
}

### PCA - reduced wavenumber - without outlier
if (select_outlier == T) {
  for (Stats in 1:length(Type)) {
    print(Names[Stats])
    PCA <- prcomp(na.omit(Data[[paste(Type[Stats],"_outl", sep="")]][!is.na(Data$Groups),Sample]),
                  center = TRUE,
                  rank. = PC_PCA, 
                  scale. = TRUE )
    
    # remove "NAs" for correct colours in the plot
    PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
    #print(Data$Groups)
    #print(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])])
    print(summary(PCA))
    
    Data[[paste("Var_", Type[Stats],"_outl", sep="")]] <- PCA$sdev^2/sum(PCA$sdev^2)
    
    plot(x = 1:PC, 
         y = cumsum(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][1:PC]*100),
         type = "b", 
         main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
         xlab = "Number of Principal Components",
         ylab = "Explained Variance [%]",
         col = "dimgrey", 
         pch = 21, 
         bg = "darkgrey",
         font = 2, 
         font.lab = 2,
         lab = c(10,20,20),
         xaxt = "n")
    axis(1,xaxp=c(1,PC,PC-1))
    
    grid(lwd = 0.8)
    
    pairs(PCA$x[,1:PCs],
          #col = Data$Groups,
          col = Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])],
          main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
          pch = 19)
    
    par(xpd = T) 
    
    legend("bottomright", 
           legend = legend, 
           pch = 19, 
           col = unique(Data$Groups), 
           inset = 0.05,
           cex = 0.7,
           horiz = T, 
           bty = "n")
    
    par(xpd = F)
    
    plot(x = PCA$x[,x_PC], 
         y = PCA$x[,y_PC],
         xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][x_PC]*100,1),"%]", sep = " "),
         ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][y_PC]*100,1),"%]", sep = " "),
         main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
         #col = Data$Groups[!is.na(Data$Groups)],
         col = Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])],
         pch = 19, 
         font.lab = 2)
    
    grid(lwd = 0.8)
    
    abline(h = 0, v = 0)
    
    legend("bottomleft", 
           legend = legend, 
           pch = 20, 
           col = unique(Data$Groups[!is.na(Data$Groups)]), 
           inset = 0.01,
           bty = "n")
    
    if (plot_loadings == T) {
      #plot.load(main=paste("Method: ",Names[Stats], " (Outl)", sep=""))
      
      DataLoadings <- data.frame(WaveNum=Data$Wavenumber_min_max)
      for (loadingNR in (1:loadings)) {
        #print(paste("PC",loadingNR,sep=""))
        DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
        colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
      }
      
      DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
      
      matplot(x = DataLoadingsRev[1],
              y = DataLoadingsRev[-1],
              type = "l", 
              main=paste("Loadings: ",Names[Stats], " (Outl)", sep=""),
              xlab = "Wavenumber [1/cm]",
              ylab = "", 
              #xaxt ="n",
              xaxs = "i",
              xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)))
      
      abline(0,0, col="black", lwd=1.5)
      legend("topleft", 
             legend = colnames(DataLoadings)[2:(loadings+1)], 
             pch = 20, 
             col = c(1:loadings), 
             inset = 0.01,
             bty = "n")
      grid(lwd = 0.8)
      
    }
  }
}


### PCA - reduced wavenumber/area 
if (remove_area == T) {

  Type = Type_a
  Sample <- c(1:length(Data$Wavenumber_area))
  
  for (Stats in 1:length(Type)) {
    print(Names[Stats])
    PCA <- prcomp(Data[[Type[Stats]]][!is.na(Data$Groups),Sample],
                  center = TRUE,
                  rank. = PC_PCA, 
                  scale. = TRUE )
    
    # remove "NAs" for correct colours in the plot
    PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
    
    print(summary(PCA))
    
    Data[[paste("Var_", Type[Stats])]] <- PCA$sdev^2/sum(PCA$sdev^2)
    
    plot(x = 1:PC, 
         y = cumsum(Data[[paste("Var_", Type[Stats])]][1:PC]*100),
         type = "b", 
         main = paste("Method: ",Names[Stats], sep=""),
         xlab = "Number of Principal Components",
         ylab = "Explained Variance [%]",
         col = "dimgrey", 
         pch = 21, 
         bg = "darkgrey",
         font = 2, 
         font.lab = 2,
         lab = c(10,20,20),
         xaxt = "n")
    axis(1,xaxp=c(1,PC,PC-1))
    
    grid(lwd = 0.8)
    
    pairs(PCA$x[,1:PCs],
          col = Data$Groups,
          main = paste("Method: ",Names[Stats], sep=""),
          pch = 19)
    
    par(xpd = T) 
    
    legend("bottomright", 
           legend = legend, 
           pch = 19, 
           col = unique(Data$Groups), 
           inset = 0.05,
           cex = 0.7,
           horiz = T, 
           bty = "n")
    
    par(xpd = F)
    
    plot(x = PCA$x[,x_PC], 
         y = PCA$x[,y_PC],
         xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats])]][x_PC]*100,1),"%]", sep = " "),
         ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats])]][y_PC]*100,1),"%]", sep = " "),
         main = paste("Method: ",Names[Stats], sep=""),
         col = Data$Groups[!is.na(Data$Groups)], 
         pch = 19, 
         font.lab = 2)
    
    grid(lwd = 0.8)
    
    abline(h = 0, v = 0)
    
    legend("bottomleft", 
           legend = legend, 
           pch = 20, 
           col = unique(Data$Groups[!is.na(Data$Groups)]), 
           inset = 0.01,
           bty = "n")
    
    # remove outlier
    Data[paste(Type[Stats],"_outl",sep="")] <- Data[paste(Type[Stats])]
    Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]] <- Data$Groups
    
    if(select_outlier == T){
      print(paste(Stats, "Select outlier: ESC to exit"))
      outlier = identify( x = PCA$x[,x_PC], y = PCA$x[,y_PC])
      print(outlier)
      Data[[paste(Type[Stats],"_outl",sep="")]][outlier,] = NA
      Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][outlier] = NA
      #print(Data$Groups)
      #print(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])
      
    }
    
    if (plot_loadings == T) {
      
      DataLoadings <- data.frame(WaveNum=Data$Wavenumber_area)
      
      for (loadingNR in (1:loadings)) {
        #print(paste("PC",loadingNR,sep=""))
        DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
        colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
      }
      
      DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
      
      matplot(x = DataLoadingsRev[1],
              y = DataLoadingsRev[-1],
              type = "l", 
              main = paste("Loadings: ",Names[Stats], sep=""), 
              ylab = "",
              xlab = "Wavenumber [1/cm]",
              xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)),
              xaxs = "i")
      abline(0,0, col="black", lwd=1.5)
      legend("topleft", 
             legend = colnames(DataLoadings)[2:(loadings+1)], 
             pch = 20, 
             col = c(1:loadings), 
             inset = 0.01,
             bty = "n")
      grid(lwd = 0.8)
      
    }
  }
}

### PCA - reduced wavenumber/area - without outlier
if (remove_area == T) {
  
  if (select_outlier == T) {
    for (Stats in 1:length(Type)) {
      print(Names[Stats])
      PCA <- prcomp(na.omit(Data[[paste(Type[Stats],"_outl", sep="")]][!is.na(Data$Groups),Sample]),
                    center = TRUE,
                    rank. = PC_PCA, 
                    scale. = TRUE )
      
      # remove "NAs" for correct colours in the plot
      PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
      #print(Data$Groups)
      #print(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])])
      print(summary(PCA))
      
      Data[[paste("Var_", Type[Stats],"_outl", sep="")]] <- PCA$sdev^2/sum(PCA$sdev^2)
      
      plot(x = 1:PC, 
           y = cumsum(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][1:PC]*100),
           type = "b", 
           main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
           xlab = "Number of Principal Components",
           ylab = "Explained Variance [%]",
           col = "dimgrey", 
           pch = 21, 
           bg = "darkgrey",
           font = 2, 
           font.lab = 2,
           lab = c(10,20,20),
           xaxt = "n")
      axis(1,xaxp=c(1,PC,PC-1))
      
      grid(lwd = 0.8)
      
      pairs(PCA$x[,1:PCs],
            #col = Data$Groups,
            col = Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])],
            main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
            pch = 19)
      
      par(xpd = T) 
      
      legend("bottomright", 
             legend = legend, 
             pch = 19, 
             col = unique(Data$Groups), 
             inset = 0.05,
             cex = 0.7,
             horiz = T, 
             bty = "n")
      
      par(xpd = F)
      
      plot(x = PCA$x[,x_PC], 
           y = PCA$x[,y_PC],
           xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][x_PC]*100,1),"%]", sep = " "),
           ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][y_PC]*100,1),"%]", sep = " "),
           main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
           #col = Data$Groups[!is.na(Data$Groups)],
           col = Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])],
           pch = 19, 
           font.lab = 2)
      
      grid(lwd = 0.8)
      
      abline(h = 0, v = 0)
      
      legend("bottomleft", 
             legend = legend, 
             pch = 20, 
             col = unique(Data$Groups[!is.na(Data$Groups)]), 
             inset = 0.01,
             bty = "n")
      
      if (plot_loadings == T) {
        
        DataLoadings <- data.frame(WaveNum=Data$Wavenumber_area)
        for (loadingNR in (1:loadings)) {
          #print(paste("PC",loadingNR,sep=""))
          DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
          colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
        }
        
        DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
        
        matplot(x = DataLoadingsRev[1],
                y = DataLoadingsRev[-1],
                type = "l", 
                main=paste("Loadings: ",Names[Stats], " (Outl)", sep=""),
                xlab = "Wavenumber [1/cm]",
                ylab = "", 
                #xaxt ="n",
                xaxs = "i",
                xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)))
        
        abline(0,0, col="black", lwd=1.5)
        legend("topleft", 
               legend = colnames(DataLoadings)[2:(loadings+1)], 
               pch = 20, 
               col = c(1:loadings), 
               inset = 0.01,
               bty = "n")
        grid(lwd = 0.8)
        
      }
    }
  }
}

end_time = format(Sys.time(), '%X') # get endtime of analysis
# print time difference between starttime and endtime
print(paste("Time used for analysis:", round(as.difftime(end_time, units = "mins")-as.difftime(start_time, units = "mins"),digits=2),"minutes"))


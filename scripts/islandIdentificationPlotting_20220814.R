######################################
### INSTALL PACKAGES & LOAD FUNCTIONS
.libPaths()
packages_needed <- c("ggplot2", "tidyverse")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

# if running on a cluster change to YES
CLUSTER <- "NO"

################################################################################

################################################################################
# Specify working directory and input files
## EVERYTHING IN THE SECTION WILL NEED TO BE UPDATED TO FIT YOUR NEEDS

if(CLUSTER == "YES"){
  args = commandArgs(trailingOnly=TRUE)
  WORKDIR <- args[1]
  DATADIR <- args[2]
  CHROMFILE <- args[3]
  METAFILE <- args[4]
  OUTFILE <- args[5]
  N <- as.numeric(args[6])
  FOCALPOP <- args[7]
  NLOCI_ISLAND <- as.numeric(args[8])
  STEP_SIZE <- as.numeric(args[9])
  WIND_SIZE <- as.numeric(args[10])
  FST_THRESHOLD <- as.numeric(args[11])
  
} else {
  
  # set main working directory - at one point in the code we have to change the working directory so we just want
  # the ability to switch back to our main working directory 
  WORKDIR <- "/fs/cbsubscb16/storage/rkc/"
  
  # Set directory for where your angsd fst output files are located. These should be files that contain data for
  # the whole genome not by chromosome
  DATADIR <- "angsd/fst/"
  
  # Tab-delimited table that has chrom name (NC_XXXXXX.X) and simplified name ( Chr 1, etc.)
  CHROMFILE <- "sample_lists/chrom_meta_data.txt" 
  
  # Full path to sampling location metadata with two columns: pop (populations in the data set) and region 
  # (general region for the populations) this second column may just be a single region depending on your sampling
  METAFILE <- tibble(Region = rep("Alaska",5), Loc = c("AI", "EastBering", "GOA", "NorthBering", "SEAK"))
  
  # Directory to where figures should be written to
  OUTFILE <- "/figures/fst/"
  
  # proportion of loci to subset per chromosome 
  # (I do one order of magnitude smaller than what I want for the final percent of loci genome-wide)
  # (so for top 0.1% of all loci I do 0.01% for N)
  N <- 0.01
  
  # For plotting example of loci called in top % 
  FOCALPOP <- "GOA"
  
  ## THESE PARAMETERS MAY NEED CHANGING DEPENDING THE CHARACTERISTICS OF YOUR ISLANDS ##
  # Set three parameters: 1) number of loci in top 0.1% that would constitute an interesting island, 2) the step size along the genome to search
  # and 3) the window size to look in
  NLOCI_ISLAND <- 7
  STEP_SIZE <- 5000
  WIND_SIZE <- 10000
  
  # set the threshold for subsetting the fst values
  FST_THRESHOLD <- 0.25
}

################################################################################

################################################################################
# Read in files

chrom_df <- read.table(paste0(WORKDIR, CHROMFILE), header = TRUE)
meta_df <- METAFILE

# Specify the order of some factors for plotting later
meta_df$Loc <- factor(meta_df$Loc, levels = meta_df$Loc)
#meta_df$region <- factor(meta_df$region, levels = unique(meta_df$region))
POPLIST <- unique(meta_df$Loc)

# THIS MAY NEED UPDATING DEPENDING ON YOUR COLOR PALETTE NEEDS
mypalette <- c("navyblue","cornflowerblue","#006d2c", "#31a354", "#74c476")
meta_df$mypalette <- mypalette
meta_df$mypalette <- factor(meta_df$mypalette, levels = meta_df$mypalette)

################################################################################

################################################################################
# Part 1: Create a concatenated dataframe and save it as a text file
# This part of the code is just parsing the angsd files, concatenating them, and adding some information that will make plotting easier. 

setwd(DATADIR)
ALLDATAFILENAMES <- Sys.glob("*")
allData_file_list <- as.list(ALLDATAFILENAMES)

fst_df <- allData_file_list %>%
  set_names(nm = ALLDATAFILENAMES) %>%
  map_dfr(
    ~ read_delim(.x, skip = 1,col_types = cols(), col_names = c("region", "chr", "midpos", "nsites", "fst"), delim = "\t"),
    .id = "comparison"
  )

head(fst_df)
setwd(WORKDIR)

# Append information about simplified chromosome names 
fst_df <- fst_df %>% 
  filter(chr %in% chrom_df$chr)
fst_df <- left_join(fst_df, chrom_df, by = "chr")
fst_df$chr_num <- factor(fst_df$chr_num, levels = chrom_df$chr_num)
fst_df$marker <- paste0(fst_df$chr, "_", fst_df$midpos)

# end Part 1
################################################################################

################################################################################
# Part 2: Identify top XX% of high FST loci

# 2.a. Subsetting - subset top XX% percent of high FST loci 
head(fst_df)

df.top.perChrom <- NULL
for(i in 1:length(unique(fst_df$chr))){
  top.chrom <- fst_df[fst_df$chr == unique(fst_df$chr)[i],]
  df.top <- top.chrom[top.chrom$fst > quantile(top.chrom$fst,prob=1-N/100),]
  df.top <- df.top[df.top$fst > FST_THRESHOLD,]
  df.top.perChrom <- rbind(df.top.perChrom, df.top) 
}

## Use this to troubleshoot what the value of N should be for the final percent of genome-wide markers
nMarkers <- length(unique(df.top.perChrom$marker))
nUniqueMarkers <- length(unique(fst_df$marker))
percentMarkers <- round((length(unique(df.top.perChrom$marker))/nUniqueMarkers)*100, digits = 3) # percent of unique markers
percentMarkersAcrossComp <- round((length(unique(df.top.perChrom$marker))/nrow(fst_df))*100, digits = 3) # percent of all markers

df.top.perChrom$comparison_marker <- paste0(df.top.perChrom$comparison, "_", df.top.perChrom$marker)

# 2.b. Plotting - subset for a focal pop to plot what loci are in the top XX%
###### THIS PART MAY HAVE TO CHANGE DEPENDING ON YOUR NAMING SYNTAX #######
plot_df <- fst_df[grep(FOCALPOP, fst_df$comparison),] %>%
  mutate(temp1 = gsub(".fst.txt", "", comparison)) %>%
  mutate(temp2 = gsub(FOCALPOP, "", temp1)) %>%
  mutate(POP2 = gsub("-", "", temp2))

# Make the POP2 column in the dataframe into a factor with a specific order
plot_df$POP2 <- factor(plot_df$POP2, levels = meta_df$Loc)

plot_df$comparison_marker <- paste0(plot_df$comparison, "_", plot_df$marker)
plot_df$topMarker <- "no"
plot_df$topMarker[plot_df$comparison_marker %in% unique(df.top.perChrom$comparison_marker)] <- "yes"

plot_df <- plot_df %>%
  mutate(midpos_Mb = midpos/1000000)

# get rid of negative fst values
plot_df$fst_neg0 <- plot_df$fst
plot_df$fst_neg0[plot_df$fst_neg0 < 0] <- 0

theme_set(
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
    axis.text.y = element_text(angle = 0, size = 7),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

for(i in 1:length(unique(plot_df$chr_num))){
  plotChrom <- unique(plot_df$chr_num)[i]
  manhplot <- ggplot() +
    geom_point(data = plot_df[plot_df$chr_num == plotChrom,],
               aes(x = midpos_Mb, y = fst_neg0, color = topMarker, alpha = topMarker)) +
    scale_color_manual(values = c("black", "red")) +
    scale_alpha_manual(values = c(0.5, 1)) +
    facet_grid(POP2~., scales = "free_x") +
    ylab(expression(italic(F[ST]))) +
    xlab("Chromosome position (Mb)") +
    ggtitle(paste(FOCALPOP, plotChrom, "top", percentMarkers,"% of high FST markers")) + ## TITLE WILL NEED CHANGING
  # set tick mark spacing
    scale_y_continuous(breaks = c(0.0, 0.5, 1)) +
    scale_x_continuous(breaks = seq(from=0, to=max(plot_df$midpos_Mb[plot_df$chr_num == plotChrom]), by = 1))
  ggsave(paste0(WORKDIR, OUTFILE, FOCALPOP,"_", plotChrom, "_fst_top_", N,"Markers.jpg"), plot = manhplot, width = 10, height = 10, units = "in")
}

################################################################################

################################################################################
# Part 3: Identify interesting windows

# step through genome to identify interesting islands

# the only problem is the final window .. either the window size can be constant between chromosomes except the final
# window on each chromosome because of variation in chromosome size or the window size can vary depending on the size of the chromosome
# which would give all windows within a chromosome the same size, but would differ between chromosomes


island_windows_df <- NULL
chromosomes <- unique(fst_df$chr)
comparisons <- unique(fst_df$comparison)
# step through each chromosome
chr_count <- 1
for(i in 1:length(chromosomes)){ 
  start_chr <- 1
  end_chr <- max(fst_df$midpos[fst_df$chr_num == paste0("chr_",i)]) # last snp location on the chromosome
  step_size <- STEP_SIZE 
  wind_size <- WIND_SIZE
  starts <- seq(from = start_chr, to = end_chr, by = step_size) # vector of start locations for the steps
  ends <- starts + wind_size  # vector of end locations for the steps 
  ends <- ends[1:(length(which(ends < end_chr))+1)] # only need those windows that end before the end of the chromosome
  starts <- starts[1:length(ends)] # grab matching starts for final end vector
  comp_count <- 0
  start_time <- Sys.time()
  # step through each pairwise comparison
  for(c in 1:length(comparisons)){
    focal_comparison_df <- fst_df[fst_df$comparison == comparisons[c] & fst_df$chr == chromosomes[i],]
    focal_comparison_df$comparison_marker <- paste0(focal_comparison_df$comparison, "_", focal_comparison_df$marker)
    # step through each window along the chromosome for that pairwise comparison
    for(j in 1:length(starts)){ 
      focal_window_df <- focal_comparison_df[focal_comparison_df$midpos >= starts[j] & focal_comparison_df$midpos < ends[j],]
      if(sum(focal_window_df$comparison_marker %in% df.top.perChrom$comparison_marker) >= NLOCI_ISLAND){
        # output a unique identifier for the window that is chromosome_windStart and the comparison its found in
        write.table(t(c(paste0(chromosomes[i], "_", starts[j]), chromosomes[i], starts[j], ends[j], comparisons[c])), paste0(WORKDIR, "islands_fstcutoff_", FST_THRESHOLD, "_nloci", NLOCI_ISLAND, "_stepSize", STEP_SIZE, "_windSize", WIND_SIZE, ".txt"),
                    col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
      } 
    } # end for loop stepping through each window along the chromosome for that pairwise comparison
    comp_count <- comp_count + 1
  } # end for loop stepping through each pairwise comparison
  end_time <- Sys.time()
  (end_time - start_time)
  print(paste(chromosomes[i], chr_count))
  chr_count <- chr_count + 1
} # end for loop stepping through each chromosome

################################################################################


################################################################################
### PLOTTING 

# upload file that has called island information
ISLANDFILE <- paste0("islands_fstcutoff_", FST_THRESHOLD, "_nloci", NLOCI_ISLAND, "_stepSize", STEP_SIZE, "_windSize", WIND_SIZE, ".txt")
island_df <- as.data.frame(read.table(paste0(WORKDIR, ISLANDFILE), header = FALSE))
colnames(island_df) <- c("window", "chr", "start", "end", "comparison")
island_df$island <- "Yes"

windows_df <- NULL
for(i in 1:length(unique(chrom_df$chr))){ 
  start_chr <- 1
  end_chr <- max(fst_df$midpos[fst_df$chr_num == paste0("chr_",i)]) # last snp location on the chromosome
  step_size <- STEP_SIZE 
  wind_size <- WIND_SIZE
  starts <- seq(from = start_chr, to = end_chr, by = step_size) # vector of start locations for the steps
  ends <- starts + wind_size  # vector of end locations for the steps 
  ends <- ends[1:(length(which(ends < end_chr))+1)] # only need those windows that end before the end of the chromosome
  starts <- starts[1:length(ends)] # grab matching starts for final end vector
  chromosome_window_df <- cbind(paste0(unique(chrom_df$chr)[i], "_", starts), unique(chrom_df$chr)[i], starts, ends)
  windows_df <- rbind(windows_df, chromosome_window_df)
}
windows_df <- as.data.frame(windows_df)[, 1:4]
windows_df$starts <- as.numeric(windows_df$starts)
windows_df$ends <- as.numeric(windows_df$ends)
head(windows_df)
colnames(windows_df) <- c("window", "chr", "starts", "ends")
plot_df <- left_join(windows_df, island_df[c("window", "island", "comparison")], by = "window")

plot_df$island[is.na(plot_df$island)] <- "No"
plot_df$island[plot_df$island == "No"] <- 0.1
plot_df$island[plot_df$island == "Yes"] <- 0.5
plot_df$island <- as.factor(plot_df$island)
plot_df$starts <- as.numeric(plot_df$starts)
plot_df <- plot_df %>%
  mutate(starts_Mb = starts/1000000)

for(i in 1:length(unique(plot_df$chr))){
  plotChrom <- unique(plot_df$chr)[i]
  chrom_num <- chrom_df$chr_num[i]
  manplot <- ggplot() +
    geom_point(data = plot_df[plot_df$chr == plotChrom,], 
               aes(x = starts_Mb, y = island, col = island ), size = 3) +
    scale_color_manual(values = c("black", "red")) +
    ylab("Island yes or no") +
    xlab("Chromosome position (Mb)") +
    ggtitle(paste(chrom_num, "-", NLOCI_ISLAND, "loci - step size", STEP_SIZE, "- window size", WIND_SIZE)) + ## TITLE WILL NEED CHANGING
    theme(panel.border = element_blank(), text = element_text(size = 25),
          axis.line = element_line(colour = "black"))+
    # set tick mark spacing
    # scale_y_continuous(breaks = c(0.0, 0.5, 1)) +
    scale_x_continuous(breaks = seq(from=0, to=max(plot_df$starts_Mb[plot_df$chr == plotChrom]), by = 2))
  
  ## FILENAME MIGHT NEED CHANGING
  ggsave(paste0(WORKDIR, OUTFILE,"_overallFst_", chrom_num, "_fst_top_", N,"Markers_", NLOCI_ISLAND, "loci_", STEP_SIZE, "kbstep_", WIND_SIZE,"kbwind.jpg"), plot = manplot, width = 10, height = 10, units = "in")
}

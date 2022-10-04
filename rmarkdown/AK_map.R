## Plot Alaska site map in separate window. 

color = c("#F8766D", "#E58700", "#C99800", "#A3A500", "#6BB100", "#00BA38", "#00BF7D", "#00C0AF", "#00BCD8", "#00B0F6", "#619CFF", "#B983FF", "#E76BF3", "#FD61D1", "#FF67A4")
BSEA_coords <- sample_table_BSEA_GOA %>% 
  dplyr::select(population, StartLonDD, StartLatDD)
BSEA_coords <- unique(BSEA_coords) %>%  
  arrange(population)

maps::map("worldHires","USA", xlim=c(-180,-140),ylim=c(54.5,71.2), col="gray90", fill=TRUE) #plot the region of USA containing all sites  
points(BSEA_coords$StartLonDD, BSEA_coords$StartLatDD, pch=c(rep(c(15,16,17,18),7), 15, 16), cex=1.5, col = alpha(color,0.7))  #plot my sample sites 
text(BSEA_coords$StartLonDD, BSEA_coords$StartLatDD, labels = BSEA_coords$population, pos = 2, offset = 0.5)

library(scales)

#display ggplot2 default hex color codes from 1 to 8
for(i in 1:15){
  print(hue_pal()(i))
}
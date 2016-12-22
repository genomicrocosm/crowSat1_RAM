setwd("/Users/jackdaw/Dropbox/Heterochromatin/Submissions/GenomeResearch/supplementary_data")
#Before reading in the cmap file, get rid of all headers starting with '#', except
# '#h' 

om<-read.table("S_Up_H32.cmap", header=T) # load the optical map
om<-read.table("D_Ko_C29.cmap", header=T)

# determine the window size and create a list with all contig IDs
contig_list<-levels(as.factor(om$CMapId))
window_size<-20000

# this loop calculates the nick density per 20 kb window and saves it in a separate 
# table

output<-data.frame()
for (i in 1:length(contig_list)){
  c<-subset(om, CMapId==contig_list[i]) # set up the contig for the loop
  window_nr<-c[1,2] %/% window_size
  
 bed<-matrix(nrow=window_nr, ncol=2)
 from<-c[2,6] # use the first nick label position as startig point
 for (j in 1:window_nr){
   bed[j,1]<-from
   to<-from+window_size
   bed[j,2]<-to
   from<-to
 }
 
 suboutput<-matrix(nrow=window_nr, ncol=4)
 for (k in 1:length(bed[,1])){
   counts<-length(c$LabelChannel[c$Position > bed[k,1] & c$Position < bed[k,2]])
   density<-counts/window_size
   suboutput[k,1]<-as.numeric(contig_list[i])
   suboutput[k,2]<-bed[k,1]
   suboutput[k,3]<-bed[k,2]
   suboutput[k,4]<-density
 }
  output<-rbind(output,suboutput)
}

density<-output

#set a label density threshold for detecting repetitve regions. In this case it is 0.0004, which corresponds to one label every 2.5 kb. 

# first mark the windows which have a high density and save them in a seperate vector
outlier_windows<-c()
for (i in 1:length(density$V1)){
  if (density$V4[i] > 0.0004){
    outlier_windows[i]<-1
  }  
  else {
    outlier_windows[i]<-0
  }
}

all_outliers<-(cbind(density,outlier_windows))

#then find the spots in the map where 5 consecutive windows show a higher label density
# and get the corresponding contigs. 

outlier_clusters<-c()
for (i in 5:length(all_outliers$outlier_windows)){
if ((all_outliers$outlier_windows[i] == 1) & (all_outliers$outlier_windows[i-4]==all_outliers$outlier_windows[i])) {
  outlier_clusters[i]<-1
  }
else {
  outlier_clusters[i]<-0
}
}


repeat_contigs<-unique(all_outliers$V1[which(outlier_clusters==1)])
repeat_contigs

H32_repeat_contigs<-repeat_contigs
C29_repeat_contigs<-repeat_contigs

# read in repetitve contigs from 'between-nick-site' method

H32_repeat_contis_frag_length<-scan("gm_rep_contigs_frag_length_S_Up_H32")
C29_repeat_contis_frag_length<-scan("gm_rep_contigs_frag_length_D_Ko_C29")
repeat_contigs<-H32_repeat_contis_frag_length
repeat_contigs<-C29_repeat_contis_frag_length

# now read in the xmap alignment file (OM vs LR and SR assembly) and check where the repetitive OM contigs align

xmap<-read.table("S_Up_H32_map_vs_illumina.xmap", header=T)
xmap<-read.table("S_Up_H32_map_vs_pacbio.xmap", header=T)
xmap<-read.table("D_Ko_C29_map_vs_illumina.xmap", header=T)
xmap<-read.table("D_Ko_C29_map_vs_pacbio.xmap", header=T)

output<-matrix(data=NA, ncol=8)

  for (i in 1:length(repeat_contigs)){
    for (j in 1:length(xmap$RefContigID[which(xmap$QryContigID==repeat_contigs[i])])){
      line<-c(xmap$RefContigID[which(xmap$QryContigID==repeat_contigs[i])][j], 
              repeat_contigs[i],
              xmap$QryStartPos[which(xmap$QryContigID==repeat_contigs[i])][j],
              xmap$QryEndPos[which(xmap$QryContigID==repeat_contigs[i])][j],
              xmap$RefStartPos[which(xmap$QryContigID==repeat_contigs[i])][j],
              xmap$RefEndPos[which(xmap$QryContigID==repeat_contigs[i])][j],
              xmap$RefLen[which(xmap$QryContigID==repeat_contigs[i])][j],
              as.character(xmap$Orientation[which(xmap$QryContigID==repeat_contigs[i])][j]))
      print(line)
      output<-rbind(output, line)
    }
  }

output_map_vs_illumina_alignment<-output
output_map_vs_pacbio_alignment<-output


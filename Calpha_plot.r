library(ggplot2)
library(plyr)



#load distances
dist<-read.table("CA_dist_107.dat",header=T)
distCryst<-read.table("CA_dist_crystal.dat",header=T)

#how to subdivide data set
div<- 713
dist_1<-dist[1:div,]
dist_2<-dist[(div+1):(div*2),]
dist_3<-dist[(div*2+1):(div*3),]
dist_4<-dist[(div*3+1):(div*4),]

rownames(dist_1)<-seq(length=(nrow(dist_1)))
rownames(dist_2)<-seq(length=(nrow(dist_2)))
rownames(dist_3)<-seq(length=(nrow(dist_3)))
rownames(dist_4)<-seq(length=(nrow(dist_4)))

distCryst_1<-distCryst[1:div,]
distCryst_2<-distCryst[(div+1):(div*2),]
distCryst_3<-distCryst[(div*2+1):(div*3),]
distCryst_4<-distCryst[(div*3+1):(div*4),]

rownames(distCryst_1)<-seq(length=(nrow(distCryst_1)))
rownames(distCryst_2)<-seq(length=(nrow(distCryst_2)))
rownames(distCryst_3)<-seq(length=(nrow(distCryst_3)))
rownames(distCryst_4)<-seq(length=(nrow(distCryst_4)))

#get pore residues
dist_1p <-dist_1[c(426:432,452:470),]
dist_2p <-dist_2[c(426:432,452:470),]
dist_3p <-dist_3[c(426:432,452:470),]
dist_4p <-dist_4[c(426:432,452:470),]


dist_allp <-data.frame(dist_1p,dist_2p["CA_dist_to_center"],dist_3p["CA_dist_to_center"],dist_4p["CA_dist_to_center"])
temp<- dist_allp[,-c(1,2)]
dist_allp$CA_dist_to_center_sd<-apply(temp,1,sd)
dist_allp$CA_dist_to_center_se <-apply(temp,1,sd)/sqrt((length(temp)))


dist_allp$mean_CA_dist_center <-rowMeans(dist_allp[,-c(1,2)])
#smallest residues appear at the top of graph and cloest distances on the ned
#factored at highest level
dist_allp$Residue_number <- factor(dist_allp$Residue_number, levels=rev(dist_allp$Residue_number))

distCryst_1p <-distCryst_1[c(426:432,452:470),]
distCryst_2p <-distCryst_2[c(426:432,452:470),]
distCryst_3p <-distCryst_3[c(426:432,452:470),]
distCryst_4p <-distCryst_4[c(426:432,452:470),]

distCryst_allp <-data.frame(distCryst_1p,distCryst_2p["CA_dist_to_center"],distCryst_3p["CA_dist_to_center"],distCryst_4p["CA_dist_to_center"])
temp<- distCryst_allp[,-c(1,2)]
distCryst_allp$CA_dist_to_center_crystal_sd<-apply(temp,1,sd)
distCryst_allp$CA_dist_to_center_crystal_se <-apply(temp,1,sd)/sqrt((length(temp)))

distCryst_allp$mean_CA_dist_center_crystal <-rowMeans(distCryst_allp[,-c(1,2)])


distCryst_allp$Residue_number <- factor(distCryst_allp$Residue_number, levels=rev(distCryst_allp$Residue_number))

dist_comb_p <-data.frame(distCryst_allp[,c("Residue_name","Residue_number","mean_CA_dist_center_crystal","CA_dist_to_center_crystal_sd","CA_dist_to_center_crystal_se")],
                         dist_allp[,c("mean_CA_dist_center","CA_dist_to_center_sd","CA_dist_to_center_se")])


#make graph

dist_pore_graph <- ggplot(dist_comb_p) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  geom_point(aes(Residue_number,mean_CA_dist_center)) +
  geom_line(aes(as.numeric(Residue_number),mean_CA_dist_center)) +
  geom_point(aes(Residue_number,mean_CA_dist_center_crystal), color="red") +
  geom_line(aes(as.numeric(Residue_number),mean_CA_dist_center_crystal), color="red") +
  scale_y_reverse()+
  coord_flip()

dist_pore_graph_single <- ggplot(dist_allp) +
  theme_bw() +
  #  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  geom_point(aes(Residue_number,mean_CA_dist_center)) +
#  geom_smooth(aes(as.numeric(Residue_number),mean_CA_dist_center),span=.078,se=F)+
  geom_line(aes(as.numeric(Residue_number),mean_CA_dist_center)) +
  scale_y_reverse()+
  coord_flip()

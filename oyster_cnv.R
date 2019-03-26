# Script used for analysis of CNV in oyster genomic data

#install.packages("Rcpp")
# install.packages("scales")
# install.packages("tidyverse")
# install.packages("stringr")
library(tidyverse)
library(RColorBrewer)
library(stringr)
library(ggplot2)
#library(plyr)
#library(dplyr)
library(scales)

#gene data - sizes but not more info
oysterdup <- read.table("/Users/tejashree/Documents/Projects/cnv/delly/oysterduplicate_sort.bed",stringsAsFactors = FALSE)
oysterdup$l <- oysterdup$V3 - oysterdup$V2
ggplot(oysterdup, aes(l))+geom_histogram(binwidth = 60)+ylim(c(0,100))+
  xlim(c(0,10000))

#all vcf data for each individual for each duplication
oysterdup2 <- read.table("/Users/tejashree/Documents/Projects/cnv/delly/germline_nohead_dup.vcf",stringsAsFactors = FALSE)
header <- strsplit("CHROM POS ID      REF     ALT     QUAL    FILTER  
                   INFO    FORMAT  CL_1    CL_2    CL_3    CL_4    CL_5    CL_6    CLP_1   CLP_2   
                   CLP_3   CLP_4   CLP_5   CLP_6   CS_1    CS_2    CS_3    CS_5    CS_6    CS_7    
                   DEBY_1  DEBY_2  DEBY_3  DEBY_4  DEBY_5  DEBY_6  HC_1    HC_3    HC_4    HC_5    
                   HC_6    HC_7    HCVA_1 HCVA_2 HCVA_3 HCVA_4 HCVA_5 HCVA_6 HG_HG0F2       
                   HG_HG2F1        HG_HG2M5
                   HI_1    HI_2    HI_3    HI_4    HI_5    HI_6    LM_1_pool       LM_3    LM_4    
                   LM_7    LM_8    LOLA_1  LOLA_2  LOLA_3  LOLA_4  LOLA_5  LOLA_6  NEH_1   NEH_2   
                   NEH_3   NEH_4   NEH_5   NEH_6   NG_NH0H4        NG_NH2F6        NG_NH2F8        
                   NG_NH2M1        OBOYS2_1        OBOYS2_2        OBOYS2_3        OBOYS2_4        
                   OBOYS2_5        OBOYS2_6        SL_1    SL_2    SL_3    SL_4    SL_5    SL_6
                   SM_10   SM_11   SM_12   SM_7    SM_8    SM_9    UMFS_1  UMFS_2  UMFS_3  UMFS_4  
                   UMFS_5  UMFS_6", "\\s+")[[1]]
colnames(oysterdup2)<-header
oysterdup3 <-dplyr::filter(oysterdup2,FILTER=="PASS")
oysterdup3$end <- str_split(oysterdup3$INFO, ';') %>%
  map_chr(5) %>% str_split('=') %>% map_chr(2) %>% as.integer()
oysterdup3$length <- oysterdup3$end - oysterdup3$POS

#funtion to pull out number of each genotype from a col in the vcf for a sample
gett <- function(bedout_col){
  sep_out <- str_split( bedout_col, ':') %>% as.data.frame()
  g <- sep_out[1,] %>% unname() %>% t()
  g2<-g[,1]
  table(g2)
}
# gett(oysterdup3$SM_8)

#run above function on all samples (columns)
gtypes <- select(oysterdup3,CL_1:UMFS_6) %>% map_dfr(gett)
gtypes2 <-as.data.frame(t(gtypes))
colnames(gtypes2) <- c("./.", "0/0",  "0/1",  "1/1")
gtypes2$pop <- map_chr(str_split(rownames(gtypes2),"_"),1) #add pop info (no sample num)
gtypes3 <- gather(gtypes2,genotype,number,-pop) #make long

ggplot(gtypes3,aes(genotype,number,color=pop))+geom_boxplot()+
  ylim(c(0,3000)) #plot numbers but this doesn't account for variable quality - use proportions below
ggplot(gtypes3,aes(genotype,number,color=pop))+geom_boxplot()

#proportions
gtypes_p <- gtypes2 %>% mutate(sum=rowSums(select(gtypes2,("0/0":"1/1")))) %>%
  mutate(p0 = gtypes2$"0/0"/sum) %>% mutate(p01 = gtypes2$"0/1"/sum) %>%
  mutate(p1 = gtypes2$"1/1"/sum) %>% select(pop,p0,p01,p1)
gtypesp2 <- gather(gtypes_p,genotype,number,-pop)
gtypesp2$pop <- as.factor(gtypesp2$pop)
levels(gtypesp2$pop) <- c("OBOYS2","UMFS","NEH","DEBY","LOLA" ,
                          "HI","SM","HC","HCVA",  "CS", "CLP",
                          "SL","CL","LM",
                          "HG","NG")  #reorder pops
ggplot(gtypesp2,aes(genotype,number,color=pop))+geom_boxplot()#+
#  geom_point(position=position_jitterdodge(dodge.width=0.9))

anova0 <- lm(p1 ~ pop, data = gtypes_p)
anova(anova0)       
summary(anova0)

mean_gtype_prop <- gtypes_p %>% gather(gtype,prop,-pop) %>% group_by(pop,gtype) %>%
  summarise(mean_prop = mean(prop))

#funtion to pull out genotype from a col in the vcf for a sample
getg <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(1)
}

gtypes_only <- map_dfr(select(oysterdup3,CL_1:UMFS_6),getg)
gtypes_only$ID <- oysterdup3$ID
gtypes_long <- gather(gtypes_only,key=sample,value=gtype,-ID)
gtypes_long$pop <- str_split(gtypes_long$sample,'_') %>% map(1) %>% as.character()
gtypes_long$pop <- as.vector(gtypes_long$pop)
gtypes_long$num_alts <- str_split(gtypes_long$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
#adding dups in all individuals of same pop to give pop count
pop_num_alts <- gtypes_long %>% filter(!is.na(num_alts)) %>%
  group_by(pop,ID) %>% summarize(num_alts = sum(num_alts)) #try creating a new variable

pop_num_alts <- left_join(pop_num_alts,select(oysterdup3,ID,length) )
pop_num_alts_present <- filter(pop_num_alts,num_alts >0)
ggplot(pop_num_alts_present, aes(pop,length)) +geom_violin() + ylim(c(0,2500))

meanl <- group_by(pop_num_alts_present,pop) %>% summarize(mean_len = mean(length),sd = sd(length))
ggplot(meanl,aes(pop,mean_len))+geom_point()+
  geom_errorbar(aes(ymin=mean_len+sd,ymax=mean_len-sd))
#populations all have the same length / distribution of duplications

##duplications by positions in all populations all chromosomes
gtypes_pos <- map_dfr(select(oysterdup3,CL_1:UMFS_6),getg)
gtypes_pos$POS <- oysterdup3$POS
gtypes_pos_long <- gather(gtypes_pos,key=sample,value=gtype,-POS)
gtypes_pos_long$pop <- str_split(gtypes_pos_long$sample,'_') %>% map(1) %>% as.character()
gtypes_pos_long$pop <- as.vector(gtypes_pos_long$pop)
gtypes_pos_long$num_alts <- str_split(gtypes_pos_long$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
#adding dups in all individuals of same pop to give pop count
pop_num_pos_alts <- gtypes_pos_long %>% filter(!is.na(num_alts)) %>%
  group_by(pop,POS) %>%
  summarize(num_alts = sum(num_alts))
pop_num_pos_alts_present <- filter(pop_num_pos_alts,num_alts >0)
ggplot(pop_num_pos_alts_present,aes(POS,num_alts,color=pop))+geom_point()
#pulling out only CHROM & POS
chrom_pos <- oysterdup3 %>% select(CHROM, POS)
#adding CHROM column to above frame 
pop_num_pos_alts_present_chrom <- left_join(pop_num_pos_alts_present, chrom_pos, by = "POS") 
#duplications by positions per chromosome
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035780.1")
#xlims were obtained from the scaffold sizes for each chromosome. info present in the shared genome folder or NCBI. 
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr1", 
                                                               x = "Position", y = "Number") + xlim(1,65668440)+ scale_x_continuous(labels = comma)
pos_chr2 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035781.1")
ggplot(pos_chr2,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr2", 
                                                               x = "Position", y = "Number") + xlim(1,61752955)
pos_chr3 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035782.1")
ggplot(pos_chr3,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr3", 
                                                               x = "Position", y = "Number") + xlim(1,77061148)
pos_chr4 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035783.1")
ggplot(pos_chr4,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr4", 
                                                               x = "Position", y = "Number")+ xlim(1,59691872)
pos_chr5 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035784.1")
ggplot(pos_chr5,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr5", 
                                                               x = "Position", y = "Number")+ xlim(1,98698416)
pos_chr6 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035785.1")
ggplot(pos_chr6,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr6", 
                                                               x = "Position", y = "Number")+ xlim(1,51258098)
pos_chr7 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035786.1")
ggplot(pos_chr7,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr7", 
                                                               x = "Position", y = "Number")+ xlim(1,57830854)
pos_chr8 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035787.1")
ggplot(pos_chr8,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr8", 
                                                               x = "Position", y = "Number")+ xlim(1,75944018)
pos_chr9 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035788.1")
ggplot(pos_chr9,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr9", 
                                                               x = "Position", y = "Number")+ xlim(1,104168038)
pos_chr10 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035789.1")
ggplot(pos_chr10,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr10", 
                                                               x = "Position", y = "Number")+ xlim(1,32650045)

#plot them all together subsetting 4 at a time
#Need to work on how to plot them two at a time in a grid.
# chr <- unique(pop_num_pos_alts_present_chrom$CHROM) %>% sort()
# chr_subset <- chr[1:4]
# pos_chr_plots <- vector("list",length(chr_subset))
# for (i in chr_subset){
#   pos_chr_plots[[i]] <- pop_num_pos_alts_present_chrom %>% filter(CHROM == i) %>% ggplot(aes(POS,num_alts,color=pop)) +
#     geom_point()+labs( x = "Position", y = "Number")+ theme(legend.position = "none")
# }
# cowplot::plot_grid(plotlist=pos_chr_plots,ncol=2,nrow=2, labels = chr_subset)

#try adding a new column called type with inbred,wild and selected values for each population and while plotting use group=type
pop <- c("CL","CLP","CS","DEBY","HC","HCVA","HG","HI","LM","LOLA","NEH","NG","OBOYS2","SL","SM","UMFS")
type <- c("W","W","W","S","W","W","I","W","W","S","S","I","S","W","W","S")
pop_type <- data.frame(pop,type)
colnames(pop_type) <- c("pop","type")
pop_type[] <- lapply(pop_type, as.character)
pop_alts_per_chrom <- left_join(pop_alts_per_chrom,pop_type,by='pop')
#Frequency of dups per chromosome for all populations
pop_alts_per_chrom <- pop_num_pos_alts_present_chrom %>% group_by(pop,CHROM) %>% 
  summarize(num_alts = sum(num_alts))
ggplot(pop_alts_per_chrom, aes(x=CHROM,y=num_alts, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Frequency of CNVs")
# Proportion of dups per chromosome for all populations 
pop_alts_per_chrom %>% group_by(pop) %>% mutate(prop_alts = num_alts/sum(num_alts)) %>% 
  ggplot(aes(x=CHROM,y=prop_alts, color=pop, group=type)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Proportion of CNVs per chromosome")
#The plot shows grouped by results but the key is still the same so hard to understand. Need to order the key as well. 
#If want to go back to the plot without ordering remove group=
#Chr5 has most cnvs and chr10 has least in all populations. To compare the totals per population for chr5 and 10. 
pop_alts_per_chrom %>% filter(CHROM=='NC_035784.1') %>% ggplot(aes(x=pop, y=num_alts)) + 
  geom_bar(stat = "identity", fill='dark gray')+ ggtitle("Duplication freq per population on chr5") +
  xlab("Populations") + ylab("Duplication frequency")
pop_alts_per_chrom %>% filter(CHROM=='NC_035789.1') %>% ggplot(aes(x=pop, y=num_alts)) + 
  geom_bar(stat = "identity", fill='dark gray')+ ggtitle("Duplication freq per population on chr10") +
  xlab("Populations") + ylab("Duplication frequency")

# get a de-duplicated list of locus id's
ids <- unique(pop_num_alts_present$ID)

# for each id, get a binary indicator of whether it is in a pop
#and bind to one dataframe
pops <- unique(pop_num_alts_present$pop)
binaries <- pops %>% 
  map_dfc(~ ifelse(ids %in% filter(pop_num_alts_present, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles

# set column names
names(binaries) <- pops

# have a look at the data
head(binaries)  

# plot the sets with UpSetR
library(UpSetR)
upset(binaries, nsets = length(pops), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")

#The bars on the bottom left show total number of duplicate loci for that population
#The bars at the top show the count of the intersections denoted in the dot matrix below them.
#So for columns in the matrix with only 1 dot, the bar above it shows the count of unique 
#duplicated loci. 
#This would be the outer area of a Venn diagram that has not intersected with anything else. 
#Columns with 2 or more dots show the count of shared duplicates 
#(intersecting sections of a Venn diagram).

##Getting the genotype and copy number for all pop
#funtion to pull out copy num from a col in the vcf for a sample
getcn <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(8)
}
cn_only <- map_dfr(select(oysterdup3,CL_1:UMFS_6),getcn)
cn_only$ID <- oysterdup3$ID
cn_only$POS <- oysterdup3$POS
cn_only$CHROM <- oysterdup3$CHROM
#pulling out gtype and cn for each pop side by side to visually compare
gtypes_cn <- left_join(cn_only, gtypes_only, by = 'ID') 
gtypes_cn <- gtypes_cn[,order(colnames(gtypes_cn))]
#long table with both cn and gtype info
cn_long <- gather(cn_only,key=sample,value=cn,-ID, -POS, -CHROM)
cn_long$pop <- str_split(cn_long$sample,'_') %>% map(1) %>% as.character()
cn_long$pop <- as.vector(cn_long$pop)
#removing 0/0 and ./. genotypes since the cn is not real for those
#cn_gtypes_long <- left_join(cn_long, gtypes_long) %>% filter(gtype != "0/0" & gtype != "./.") %>% select(ID, CHROM, POS, pop, sample, cn, gtype, num_alts) 
cn_gtypes_long <- left_join(cn_long, gtypes_long)
#Converting cn value to 0 for genotypes 0/0 and ./. because they are assumed homologous to reference. 
cn_gtypes_long <- within(cn_gtypes_long, cn[gtype == '0/0'] <- 0)
cn_gtypes_long <- within(cn_gtypes_long, cn[gtype == './.'] <- 0)
cn_gtypes_long$cn <- as.numeric(as.character(cn_gtypes_long$cn))
#cn stats
min(cn_gtypes_long[,5], na.rm=T) #-1
max(cn_gtypes_long[,5], na.rm=T) #20025
hist(cn_gtypes_long$cn)
#cn on chr1
cn_gtypes_long_chr1 <- filter(cn_gtypes_long, CHROM == "NC_035780.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr1$cn <- as.numeric(as.character(cn_gtypes_long_chr1$cn))
cn_chr1_hmap <- ggplot(data = cn_gtypes_long_chr1, mapping = aes(x = POS,y = sample,color = -cn)) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+scale_color_viridis_c()
cn_chr1_hmap
#Chr5
cn_gtypes_long_chr5 <- filter(cn_gtypes_long, CHROM == "NC_035784.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr5$cn <- as.numeric(as.character(cn_gtypes_long_chr5$cn))
cn_chr5_hmap <- ggplot(data = cn_gtypes_long_chr5, mapping = aes(x = POS,y = sample,color = -cn)) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+scale_color_viridis_c()
cn_chr5_hmap
#Chr9
#Chr5
cn_gtypes_long_chr9 <- filter(cn_gtypes_long, CHROM == "NC_035788.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr9$cn <- as.numeric(as.character(cn_gtypes_long_chr9$cn))
cn_chr9_hmap <- ggplot(data = cn_gtypes_long_chr9, mapping = aes(x = POS,y = sample,color = -cn)) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+scale_color_viridis_c()
cn_chr9_hmap

#trials
cn_chr1_hmap2 <- ggplot(data = cn_gtypes_long_chr1, mapping = aes(x = POS,y = sample,color = cn)) + 
  geom_point(aes(cn)) + xlab(label = "Position") +  
  scale_colour_gradient(name = "Copy Number", low = "#FFFF00", high = "#000080")
cn_chr1_hmap2
tmp2 <- filter(cn_gtypes_long_chr1, sample == "CL_1") 
tmp3 <- ggplot(data= tmp2, mapping = aes(x = POS,y = sample,color = cn), stat = "identity") + geom_tile() + xlab(label = "Position")
tmp3
ggplot(tmp2, aes(POS, sample)) +
  geom_raster(aes(fill = cn), interpolate = TRUE)
tmp2 <- filter(tmp2, cn != 0) 
  ggplot(tmp2, aes(POS, sample)) + geom_jitter(aes(color = cn))
# tmp <- data.matrix(cn_gtypes_long_chr1, rownames.force = NA)
# heatmap(tmp)

  # ggplot(cn_gtypes_long_chr1, aes(POS,cn,fill=cn))+ geom_point() +
  # labs(title = "Chr1", x = "Position", y = "Copy Number") +
  # xlim(1,65668440)+ scale_x_continuous(labels = comma)

###### ANALYSIS POST ANNOTATION #######
# Annotating dups 
# Read in annotations from ref genome
ref_annot <- read.table("/Users/tejashree/Documents/Projects/cnv/annot/ref_annot", 
                        sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE)
colnames(ref_annot) <- c("LOC", "annot")
# Read in bed file of dups mapped to LOCs from ref genome
map_dup <- read.table("/Users/tejashree/Documents/Projects/cnv/annot/Oyster_Dup_gene", sep="\t" , stringsAsFactors = FALSE)
colnames(map_dup) <- c("ID", "LOC")
dup_annot <- left_join(map_dup, ref_annot, by = "LOC") 
dup_annot %>% group_by(annot) %>% summarize(count=n()) # There are 593 with no annotation
#The proportion shows that highest dups mapped on chr5 and lowest on chr10. Need to see what are the annot for those dups.
#chrom_pos_id <- 
dup_annot_chr5 <- oysterdup3 %>% select(CHROM,POS,ID) %>% left_join(dup_annot, by='ID') %>% 
  filter(CHROM=="NC_035784.1") %>% filter(!is.na(LOC)) %>% filter(!is.na(annot)) %>% 
  filter(!grepl("uncharacterized",annot))
#lectins
dup_annot_chr10 <- oysterdup3 %>% select(CHROM,POS,ID) %>% left_join(dup_annot, by='ID') %>% 
  filter(CHROM=="NC_035789.1") %>% filter(!is.na(LOC)) %>% filter(!is.na(annot)) %>% 
  filter(!grepl("uncharacterized",annot))
#toll-like receptors


##Mapping GO_IDs, EC_num, LOCs and DUPs
#Protein_ids (XP_IDs) and LOCs were pulled out from ref gff3 file using awk (detailed in log)
ref_annot_prot <- read.table("/Users/tejashree/Documents/Projects/cnv/annot/ref_annot_prot", 
                        sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE)
colnames(ref_annot_prot) <- c("LOC", "Sequence_name")
#Join this with DUP_IDs (map_dup has DUP_ID and LOC)
dup_loc_xp <- left_join(map_dup, ref_annot_prot, by="LOC") 
#Join that with XP_sequences_Cvirginica_GCF_002022765.2_GO.tab by XP 
ref_annot_go_kegg <- read.table("/Users/tejashree/Documents/Projects/oyster/exp_data/Spring2016/genome/XP_sequences_Cvirginica_GCF_002022765.2_GO.tab", 
                             sep="\t" , quote="", fill=FALSE, stringsAsFactors = FALSE, header = TRUE)
colnames(ref_annot_go_kegg) <- c("Sequence_name","Sequence_length","Sequence_description","GO_ID","Enzyme_code","Enzyme_name")
left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% head()
#What % LOCs are mapped to kegg
ref_annot_go_kegg %>% filter(!is.na(Enzyme_name)) %>% filter(Enzyme_name != "") %>% nrow() #10.4 % (100*6273)/60213
ref_annot_go_kegg %>% filter(!is.na(GO_ID)) %>% filter(GO_ID != "") %>% nrow() #58.2% (100*35081)/30213

#Extract DUP_IDs, GO_IDs
dup_go <- left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% select(ID, GO_ID) %>% unique() 
#  separate(GO_ID, sep = ";", into = paste("V", 1:13, sep = "_")) 
#What % dups mapped to GO terms
dup_go %>% filter(!is.na(GO_ID)) %>% filter(GO_ID != "") %>% nrow() # 62% (8300*100)/13387
dup_go %>% filter(!is.na(GO_ID)) %>% filter(GO_ID != "") %>% 
  write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/dup_go.txt", append = FALSE, sep = " ",quote = FALSE,
                                                                         row.names = F, col.names = FALSE)

#separate the GO_IDs and get count for each
go_vector <- as.data.frame(table(unlist(strsplit(as.character(dup_go$GO_ID), ";"))))
go_vector_sorted <-  go_vector[order(go_vector$Freq, decreasing=TRUE),] 
write.table(go_vector_sorted, "/Users/tejashree/Documents/Projects/cnv/delly/go_vector_sorted.txt", append = FALSE, sep = " ",quote = FALSE,
            row.names = F, col.names = TRUE)
write.table(go_vector_sorted$Var1, "/Users/tejashree/Documents/Projects/cnv/delly/go_only.txt", append = FALSE, sep = ",",quote = FALSE,
            row.names = F, col.names = FALSE)
#Highest freq: molecular funtion, ion binding, cellular component, signal transduction, cellular protein modification process
#Extract DUP_IDs, EC#s :pathway mapping 
dup_kegg <- left_join(dup_loc_xp, ref_annot_go_kegg, by="Sequence_name") %>% select(ID, Enzyme_name) %>% unique() 
#What % dups mapped to an EC number via kegg
dup_kegg %>% filter(!is.na(Enzyme_name)) %>% filter(Enzyme_name != "") %>% nrow() # 11% (1420*100)/13387
#separate the enzyme names and get count for each
kegg_vector <- as.data.frame(table(unlist(strsplit(as.character(dup_kegg$Enzyme_name), ";"))))
#Highest:Nucleoside-triphosphate phosphatase,Acting on peptide bonds (peptidases),Protein-serine/threonine phosphatase,
#Protein-tyrosine-phosphatase,Adenosinetriphosphatase

##DO dups exist in expanded gene families?
##Pull out dups that are annotated as the gene members that we know belong to expanded families
dplyr::filter(dup_annot, grepl('interferon-induced protein 44', annot)) %>% select('ID') %>% unique() %>% tally() #15
dplyr::filter(dup_annot, grepl('interferon-induced protein 44', annot)) %>% select('LOC') %>% unique() %>% tally() #17
ifi44_ID <- dplyr::filter(dup_annot, grepl('interferon-induced protein 44', annot)) %>% select('ID') %>% unique()
gtypes_long_pres <- gtypes_long %>% mutate(pres = ifelse(gtypes_long$num_alts > 0, 'yes', 'no'))
ifi44 <- left_join(ifi44_ID, gtypes_long_pres ,by = "ID") %>% filter(pres == 'yes')
#plotting freq of the dup IDs that map to IFI44 across all INDIVIDUALS not separated by populations.
ifi44 %>% group_by(ID) %>% tally() %>% ggplot(aes(x=ID, y=n)) + geom_col(fill='darkblue') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ifi44_sub <- select(ifi44,ID,sample,pop)
ifi44_cn <- select(cn_gtypes_long,ID,sample,pop,cn)
ifi44_cn$cn <- as.numeric(as.character(ifi44_cn$cn))
left_join(ifi44_sub,ifi44_cn) %>% ggplot(aes(pop,cn)) + geom_col()


dplyr::filter(dup_annot, grepl('GTPase IMAP family member', annot)) %>% select('ID') %>% unique() %>% tally() #23 multiple members 4,7,8
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 4', annot)) %>% select('ID') %>% unique() %>% tally() #21
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 7', annot)) %>% select('ID') %>% unique() %>% tally() #9
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 8', annot)) %>% select('ID') %>% unique() %>% tally() #5
#The individual tally doesnt add up because for some DUPs they are mapped to multiple LOCs
dplyr::filter(dup_annot, grepl('scavenger receptor', annot)) %>% select('ID') %>% unique() %>% tally() #37 (multiple types/classes)
#now find out how many of these dups are present in each population



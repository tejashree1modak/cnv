# Script used for analysis of CNV in oyster genomic data

#install.packages("Rcpp")
# install.packages("scales")
# install.packages("tidyverse")
# install.packages("stringr")
library(tidyverse)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(dplyr)

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
  group_by(pop,ID) %>% 
  summarize(num_alts = sum(num_alts))

pop_num_alts <- left_join(pop_num_alts,select(oysterdup3,ID,length) )
pop_num_alts_present <- filter(pop_num_alts,num_alts >0)
ggplot(pop_num_alts_present, aes(pop,length)) +geom_violin() + ylim(c(0,2500))

meanl <- group_by(pop_num_alts_present,pop) %>% 
  summarize(mean_len = mean(length),
            sd = sd(length))
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
pop_num_pos_alts_present_chrom <- left_join(pop_num__pos_alts_present, chrom_pos, by = "POS") 
#duplications by positions per chromosome
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035780.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr1", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035781.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr2", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035782.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr3", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035783.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr4", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035784.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr5", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035785.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr6", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035786.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr7", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035787.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr8", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035788.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr9", 
                                                               x = "Position", y = "Number")
pos_chr1 <- filter(pop_num_pos_alts_present_chrom, CHROM == "NC_035789.1")
ggplot(pos_chr1,aes(POS,num_alts,color=pop))+geom_point()+labs(title = "Chr10", 
                                                               x = "Position", y = "Number")

#plot them all together subsetting 4 at a time
chr <- unique(pop_num_pos_alts_present_chrom$CHROM) %>% sort()
chr_subset <- chr[1:4]
pos_chr_plots <- vector("list",length(chr_subset))
for (i in chr_subset){
  pos_chr_plots[[i]] <- pop_num_pos_alts_present_chrom %>% filter(CHROM == i) %>% ggplot(aes(POS,num_alts,color=pop)) +
    geom_point()+labs( x = "Position", y = "Number")+ theme(legend.position = "none")
}
cowplot::plot_grid(plotlist=pos_chr_plots,ncol=2,nrow=2, labels = chr_subset)


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
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 1), order.by = "freq")

#The bars on the bottom left show total number of duplicate loci for that population
#The bars at the top show the count of the intersections denoted in the dot matrix below them.
#So for columns in the matrix with only 1 dot, the bar above it shows the count of unique 
#duplicated loci. 
#This would be the outer area of a Venn diagram that has not intersected with anything else. 
#Columns with 2 or more dots show the count of shared duplicates 
#(intersecting sections of a Venn diagram).

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


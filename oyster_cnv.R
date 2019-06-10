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

#Labels, colors and shapes for consistent plotting 
labels = c("CL"="LA-HSCL","CLP"="CB-LSCP","CS"="DB-HSCS","DEBY"="CB-DEBY","HC"="DB-LSHC",
           "HCVA"="CB-HSHC","HG"="DB-NEHG","HI"="ME-HSHI","LM"="TX-HSLM","LOLA"="CB-LOLA",
           "NEH"="DB-NEHD","NG"="DB-NEHS","OBOYS2"="LA-OBOY","SL"="LA-LSSL","SM"="ME-LSSM","UMFS"="ME-UMFS")
sample_lables = c("CL_1"="LA-HSCL_1", "CL_2"="LA-HSCL_2","CL_3"="LA-HSCL_3","CL_4"="LA-HSCL_4","CL_5"="LA-HSCL_5","CL_6"="LA-HSCL_6",
                  "CLP_1"="CB-LSCP_1","CLP_2"="CB-LSCP_2","CLP_3"="CB-LSCP_3","CLP_4"="CB-LSCP_4","CLP_5"="CB-LSCP_5","CLP_6"="CB-LSCP_6",
                  "CS_1"="DB-HSCS_1","CS_2"="DB-HSCS_2","CS_3"="DB-HSCS_3","CS_5"="DB-HSCS_5","CS_6"="DB-HSCS_6","CS_7"="DB-HSCS_7",
                  "DEBY_1"="CB-DEBY_1","DEBY_2"="CB-DEBY_2","DEBY_3"="CB-DEBY_3","DEBY_4"="CB-DEBY_4","DEBY_5"="CB-DEBY_5","DEBY_6"="CB-DEBY_6",
                  "HC_1"="DB-LSHC_1","HC_3"="DB-LSHC_3","HC_4"="DB-LSHC_4","HC_5"="DB-LSHC_5","HC_6"="DB-LSHC_6","HC_7"="DB-LSHC_7",
                  "HCVA_1"="CB-HSHC_1","HCVA_2"="CB-HSHC_2","HCVA_3"="CB-HSHC_3","HCVA_4"="CB-HSHC_4","HCVA_5"="CB-HSHC_5","HCVA_6"="CB-HSHC_6",
                  "HG_HG0F2"="DB-NEHG_1",  "HG_HG2F1"="DB-NEHG_2","HG_HG2M5"="DB-NEHG_3",
                  "HI_1"="ME-HSHI_1","HI_2"="ME-HSHI_2","HI_3"="ME-HSHI_3","HI_4"="ME-HSHI_4","HI_5"="ME-HSHI_5", "HI_6"="ME-HSHI_6",
                  "LM_1_pool"="TX-HSLM_1", "LM_3"="TX-HSLM_3","LM_4"="TX-HSLM_4","LM_7"="TX-HSLM_7","LM_8"="TX-HSLM_8",
                  "LOLA_1"="CB-LOLA_1","LOLA_2"="CB-LOLA_2","LOLA_3"="CB-LOLA_3","LOLA_4"="CB-LOLA_4","LOLA_5"="CB-LOLA_5","LOLA_6"="CB-LOLA_6",
                  "NEH_1"="DB-NEHD_1","NEH_2"="DB-NEHD_2","NEH_3"="DB-NEHD_3","NEH_4"="DB-NEHD_4","NEH_5"="DB-NEHD_5","NEH_6"="DB-NEHD_6",
                  "NG_NH0H4"="DB-NEHS_1", "NG_NH2F6"="DB-NEHS_2","NG_NH2F8"="DB-NEHS_3","NG_NH2M1"="DB-NEHS_4", 
                   "OBOYS2_1"="LA-OBOY_1","OBOYS2_2"="LA-OBOY_2","OBOYS2_3"="LA-OBOY_3","OBOYS2_4"="LA-OBOY_4","OBOYS2_5"="LA-OBOY_5","OBOYS2_6"="LA-OBOY_6",
                  "SL_1"="LA-LSSL_1", "SL_2"="LA-LSSL_2", "SL_3"="LA-LSSL_3", "SL_4"="LA-LSSL_4", "SL_5"="LA-LSSL_5","SL_6"="LA-LSSL_6",
                  "SM_10"="ME-LSSM_1","SM_11"="ME-LSSM_2","SM_12"="ME-LSSM_3","SM_7"="ME-LSSM_4","SM_8"="ME-LSSM_5","SM_9"="ME-LSSM_6",
                  "UMFS_1"="ME-UMFS_1","UMFS_2"="ME-UMFS_2","UMFS_3"="ME-UMFS_3","UMFS_4"="ME-UMFS_4","UMFS_5"="ME-UMFS_5","UMFS_6"="ME-UMFS_6")
values = c("CL"="#a63603","CLP"="#a6d854","CS"="#08519c","DEBY"="#006d2c","HC"="skyblue",
           "HCVA"="#31a354","HG"="#6a51a3","HI"="#fbb4b9","LM"="black","LOLA"="#bae4b3",
           "NEH"="#9e9ac8","NG"="#dadaeb","OBOYS2"="#fd8d3c","SL"="#fdae6b","SM"="#7a0177","UMFS"="#f768a1")
shapes = c("CL"=17,"CLP"=17,"CS"=17,"DEBY"=16,"HC"=17,
           "HCVA"=17,"HG"=13,"HI"=17,"LM"=17,"LOLA"=16,
           "NEH"=16,"NG"=13,"OBOYS2"=16,"SL"=17,"SM"=17,"UMFS"=16)

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
oysterdup3$length <- oysterdup3$end - oysterdup3$POS #Smallest is 160bp and largest is 999,122bp
filter(oysterdup3, length > 1000) %>% nrow() #6551 So if we filter the dups based on len 
                                                    #and only keep >1kb we will only have 6551 instead of 13387

#function to pull out number of each genotype from a col in the vcf for a sample
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
# pop <- c("CL","CLP","CS","DEBY","HC","HCVA","HG","HI","LM","LOLA","NEH","NG","OBOYS2","SL","SM","UMFS")
# type <- c("W","W","W","S","W","W","I","W","W","S","S","I","S","W","W","S")
# pop_type <- data.frame(pop,type)
# colnames(pop_type) <- c("pop","type")
# pop_type[] <- lapply(pop_type, as.character)
# pop_alts_per_chrom <- left_join(pop_alts_per_chrom,pop_type,by='pop')
#Frequency of dups per chromosome for all populations
pop_alts_per_chrom <- pop_num_pos_alts_present_chrom %>% group_by(pop,CHROM) %>% 
  summarize(num_alts = sum(num_alts))
ggplot(pop_alts_per_chrom, aes(x=CHROM,y=num_alts, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Frequency of CNVs")
# normalized by chromosome size
chrom_len <- data.frame(CHROM=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                 "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
           start=c(1,1,1,1,1,1,1,1,1,1), 
           end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045))
chrom_len$len <- chrom_len$end - chrom_len$start
pop_alts_per_chrom_len <- left_join(pop_alts_per_chrom, chrom_len, by = "CHROM")
ggplot(pop_alts_per_chrom_len, aes(x=CHROM,y=(num_alts/len), color=pop)) + geom_bar(stat = "identity", fill="white") + 
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

### Upset plot of duplications in populations ###
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
#function to pull out copy num from a col in the vcf for a sample
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
filter(cn_gtypes_long_chr1, cn > 500) %>% select(cn)%>% hist()
cn_chr1_hmap <- ggplot(data = cn_gtypes_long_chr1, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 1") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr1_hmap + geom_point(data=cn_gtypes_long_chr1, aes(x=65668439/2, y=1), col="red",pch=24, cex = 3)
cn_chr1_hmap + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)
#Chr2
cn_gtypes_long_chr2 <- filter(cn_gtypes_long, CHROM == "NC_035781.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr2$cn <- as.numeric(as.character(cn_gtypes_long_chr2$cn))
cn_chr2_hmap <- ggplot(data = cn_gtypes_long_chr2, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 2") + scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr2_hmap + geom_point(data=cn_gtypes_long_chr2, aes(x=30876477, y=1), col="red",pch=24, cex = 3)
cn_chr2_hmap + geom_vline(xintercept = (61752954/2), color = "red", size=0.3)
#Chr3
cn_gtypes_long_chr3 <- filter(cn_gtypes_long, CHROM == "NC_035782.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr3$cn <- as.numeric(as.character(cn_gtypes_long_chr3$cn))
cn_chr3_hmap <- ggplot(data = cn_gtypes_long_chr3, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 3") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr3_hmap + geom_point(data=cn_gtypes_long_chr3, aes(x=(77061147/2), y=1), col="red",pch=24, cex = 3)
cn_chr3_hmap + geom_vline(xintercept = (77061147/2), color = "red", size=0.3)
#Chr4
cn_gtypes_long_chr4 <- filter(cn_gtypes_long, CHROM == "NC_035783.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr4$cn <- as.numeric(as.character(cn_gtypes_long_chr4$cn))
cn_chr4_hmap <- ggplot(data = cn_gtypes_long_chr4, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 4") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr4_hmap + geom_point(data=cn_gtypes_long_chr4, aes(x=(59691871/2), y=1), col="red",pch=24, cex = 3)
cn_chr4_hmap + geom_vline(xintercept = (59691871/2), color = "red", size=0.3)
#Chr5
cn_gtypes_long_chr5 <- filter(cn_gtypes_long, CHROM == "NC_035784.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr5$cn <- as.numeric(as.character(cn_gtypes_long_chr5$cn))
cn_chr5_hmap <- ggplot(data = cn_gtypes_long_chr5, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 5") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr5_hmap + geom_point(data=cn_gtypes_long_chr5, aes(x=(98698415/2), y=1), col="red",pch=24, cex = 3)
cn_chr5_hmap + geom_vline(xintercept = (98698415/2), color = "red", size=0.3)
#Chr6
cn_gtypes_long_chr6 <- filter(cn_gtypes_long, CHROM == "NC_035785.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr6$cn <- as.numeric(as.character(cn_gtypes_long_chr6$cn))
cn_chr6_hmap <- ggplot(data = cn_gtypes_long_chr6, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 6") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr6_hmap + geom_point(data=cn_gtypes_long_chr6, aes(x=(51258097/2), y=1), col="red",pch=24, cex = 3)
cn_chr6_hmap + geom_vline(xintercept = (51258097/2), color = "red", size=0.3)
#Chr7
cn_gtypes_long_chr7 <- filter(cn_gtypes_long, CHROM == "NC_035786.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr7$cn <- as.numeric(as.character(cn_gtypes_long_chr7$cn))
cn_chr7_hmap <- ggplot(data = cn_gtypes_long_chr7, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 7") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr7_hmap + geom_point(data=cn_gtypes_long_chr7, aes(x=(57830853/2), y=1), col="red",pch=24, cex = 3)
cn_chr7_hmap + geom_vline(xintercept = (57830853/2), color = "red", size=0.3)
#Chr8
cn_gtypes_long_chr8 <- filter(cn_gtypes_long, CHROM == "NC_035787.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr8cn <- as.numeric(as.character(cn_gtypes_long_chr8$cn))
cn_chr8_hmap <- ggplot(data = cn_gtypes_long_chr8, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 8") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr8_hmap + geom_point(data=cn_gtypes_long_chr8, aes(x=(75944017/2), y=1), col="red",pch=24, cex = 3)
cn_chr8_hmap + geom_vline(xintercept = (75944017/2), color = "red", size=0.3)
#Chr9
cn_gtypes_long_chr9 <- filter(cn_gtypes_long, CHROM == "NC_035788.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr9$cn <- as.numeric(as.character(cn_gtypes_long_chr9$cn))
cn_chr9_hmap <- ggplot(data = cn_gtypes_long_chr9, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 9") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr9_hmap + geom_point(data=cn_gtypes_long_chr9, aes(x=(104168037/2), y=1), col="red",pch=24, cex = 3)
cn_chr9_hmap + geom_vline(xintercept = (104168037/2), color = "red", size=0.3)
#Chr10
cn_gtypes_long_chr10 <- filter(cn_gtypes_long, CHROM == "NC_035789.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr10$cn <- as.numeric(as.character(cn_gtypes_long_chr10$cn))
cn_chr10_hmap <- ggplot(data = cn_gtypes_long_chr10, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 10") + scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
# cn_chr10_hmap + geom_point(data=cn_gtypes_long_chr10, aes(x=32650044/2, y=1), col="red",pch=24, cex = 3)
cn_chr10_hmap + geom_vline(xintercept = (32650044/2), color = "red", size=0.3)

#### Filteration #####
#Get dups common to all populations since they are likely artifacts 
common_dups <- pop_num_alts_present %>% group_by(ID) %>% tally(sort = TRUE) %>% head(961) %>% select(ID)
#Get dups common to all populations but present in all samples of each population
sample_num_alts <- gtypes_long %>% filter(!is.na(num_alts)) %>% filter(num_alts >0)
#counting the number of total samples
apply(gtypes_long, 2, function(x) length(unique(x))) #there are 90 samples in total
#   ID   sample    gtype      pop num_alts 
#13387       90        4       16        4
sample_num_alts %>% group_by(ID) %>% summarize(count=n()) %>% View() #No dups with count=90, 44 dups with count=89 
#to verify that there are 961 dups present in all populations
pop_num_alts_present %>% group_by(ID) %>% summarize(count=n()) %>% View()
#Criteria for filteration of common dups: 
#Out of the 961 dups that are present in all populations how many are present in >90% samples (i.e. in 81 samples) 
common_filter_dups <- 
  semi_join(sample_num_alts,common_dups, by="ID") %>% group_by(ID) %>% summarize(count=n()) %>% filter(count > 81) %>% select("ID")#298 

##Repeat Masker##
#Repeatmasker output was obtained from the ftp server where C.virginica genome files are at
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/765/GCA_002022765.4_C_virginica-3.0
#The file is GCA_002022765.4_C_virginica-3.0_rm.out.gz
#Repeat locations in the C.vir genome as shown by repeat masker out file
#only start and stop pulled out of output
repeats <- read.table("/Users/tejashree/Documents/Projects/cnv/genome/repeat_masker/Cvir_genome_repeats.txt", 
                      sep="\t" , header = TRUE, skip = 1, stringsAsFactors = FALSE)
repeats$len <- repeats$end - repeats$begin
#A bedfile of repeats was made (Cvir_genome_repeats.bed)
#The repeat output has different chromosome numbers so the file was edited to match names (eg. NC_035789.1) (Cvir_genome_repeats_mod.bed)
#To get overlaps of dups and repeats bedtools.sh was used.
#Output file was cleaned up to remove CHROM column that was present twice 
#(since I used -wo tag for bedtools intersect it gives location of both dup and repeat so CHROM appeared twice in the output)
#Due to presence of overlapping repeats, bedtools intersect was rerun with merged repeats.  
# dup_repeat_overlap <- read.table("/Users/tejashree/Documents/Projects/cnv/genome/repeat_masker/dup_repeat_overlap_len_mod.bed", 
#                       sep="\t" , stringsAsFactors = FALSE)
# colnames(dup_repeat_overlap) <- c("CHROM", "POS","end","ID","R_POS","R_end","R_ID","l")
dup_repeat_overlap <- read.table("/Users/tejashree/Documents/Projects/cnv/genome/repeat_masker/dup_repeat_merged_overlap_mod.bed", 
                                 sep="\t" , stringsAsFactors = FALSE)
colnames(dup_repeat_overlap) <- c("CHROM", "POS","end","ID","R_POS","R_end","R_ID","l")
#Number of repeats mapped to each duplicate
dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% tally() %>% View()
#Total len of repeats included in dups
overlap_total_len <- dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% summarise(sum(l))
colnames(overlap_total_len) <- c("ID","total_len")
#tried to plot but too many entries for a simple plot 
# overlap_total_len%>% ggplot(aes(x=ID, y=total_len)) + geom_col(fill='darkblue') + ylab(label = "Total length of repeats in dups") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#Formatted data to pull out only dups and join repeat len by "," to get see in a gimpse but too many entries joined so not very readable
aggregate(l~ID,data=dup_repeat_overlap,paste,collapse=",") %>% View()
#Get % of dup len covered by repeats for each dup that overlaps with repeats
percent_overlap <- oysterdup3 %>% select(ID,length) %>% left_join(overlap_total_len,by = 'ID') %>% na.omit() 
percent_overlap$percent <- (percent_overlap$total_len/percent_overlap$length)*100
#dups with >10% repeat coverage
percent_overlap %>% filter(percent > 10) %>% nrow() #filter out 1778 dups
repeat_filter_dups <- percent_overlap %>% filter(percent > 10) %>% select("ID")

#There are 28 dups present in both common_filter_dups and repeat_filter_dups 
anti_join(common_filter_dups,repeat_filter_dups) %>% nrow()
#Combining list of dups to be filtered because they are shared among >90% samples or have high repeat coverage.
#Also removing dups present in both common_filter_dups and repeat_filter_dups using distinct()
filter_dups <- rbind(common_filter_dups, repeat_filter_dups) %>% distinct() #2048
# Remove these dups from dups bedfile.
cvir_dup_bed <- oysterdup %>% select(V1:V4)
colnames(cvir_dup_bed) <- c("CHROM","start","stop","ID")
anti_join(cvir_dup_bed,filter_dups) %>% group_by(ID) %>% summarize(count=n()) %>% nrow() #11339 
anti_join(cvir_dup_bed,filter_dups) %>%
write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/cvir_filtered_dups.bed", append = FALSE, sep = "\t",quote = FALSE,
            row.names = F, col.names = FALSE)
#Dups shared by all pops post filtration
pop_shared_fil <- pop_num_alts_present_fil %>% group_by(ID) %>% tally(sort = TRUE) %>% head(608) %>% select(ID)
semi_join(cvir_dup_bed, pop_shared_fil) %>% group_by(ID) %>% summarize(count=n()) %>% nrow()
semi_join(cvir_dup_bed, pop_shared_fil) %>%
  write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/cvir_filtered_fully_shared_dups.bed", append = FALSE, sep = "\t",quote = FALSE,
              row.names = F, col.names = FALSE)

##Repeat analysis post filteration##
#Total length of dups/length of genome
#In order to get total length of genome covered by duplications 
#I merged the duplications to avoid counting overlapping dups multiple times
#This was done using bedtools merge on bluewaves
dups_fil_merged <- read.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/bluewaves/cvir_filtered_dups_merged.bed", 
                              sep="\t" , stringsAsFactors = FALSE)
colnames(dups_fil_merged) <- c("CHROM", "POS","end","count","POS_collapse","end_collapse")
#Length of merged dups
dups_fil_merged$len <- (dups_fil_merged$end - dups_fil_merged$POS) + 1
#Total number of bases of all duplications
sum(dups_fil_merged$len) #118732426
#% of genome covered by duplications
(118732426/684000000)*100 #17%  
#Zebrafish is 14.6% with 192,460,331 bp (Brown et al., 2012)
#Humans ~4% (Conrad et al., 2010)

#Filter original datasets for further analysis:  oysterdup3
oysterdup3_fil <- anti_join(oysterdup3, filter_dups)
pop_num_alts_present_fil <- anti_join(pop_num_alts_present,filter_dups)
#Count if number of dups is the same
length(unique(pop_num_alts_present_fil[["ID"]])) #11339
#Verify number of dups present in all pops
#anti_join(common_dups,common_filter_dups) %>% anti_join(repeat_filter_dups) %>% nrow() 
#both common and repeat filter criteria eliminate 353 dups from 961 that were common in all pops
pop_num_alts_present_fil %>% group_by(ID) %>% summarize(count=n()) %>% filter(count==16) %>% nrow() #608
## Upset plot of duplications in populations POST FILTERATION ##
    # get a de-duplicated list of locus id's
ids <- unique(pop_num_alts_present_fil$ID)
# for each id, get a binary indicator of whether it is in a pop and bind to one dataframe
pops <- unique(pop_num_alts_present_fil$pop)
binaries <- pops %>% 
  map_dfc(~ ifelse(ids %in% filter(pop_num_alts_present_fil, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles
# set column names
names(binaries) <- pops
# have a look at the data
head(binaries)  
# plot the sets with UpSetR
library(UpSetR)
upset(binaries, nsets = length(pops), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")
# Dups per population POST FILTRATION #
pop_sum_fil <- as.data.frame(colSums(binaries))
pop_sum_fil <- data.frame(pop = names(binaries),total_dups=colSums(binaries))
pop_sum_fil$prop <- pop_sum_fil$total_dups/11339  #number of filtered dups are 11339 
ggplot(pop_sum_fil, aes(x=pop,y=prop, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Populations", y="Proportion of total duplications per population", title ="Post filteration") + scale_color_manual(values=values,labels=labels)

##Comparison of dups within samples for inbred pop ##
anti_join(gtypes_long, filter_dups) %>% nrow()
# For HG
gtypes_long_HG <- anti_join(gtypes_long, filter_dups) %>% filter(num_alts > 0) %>% 
  filter(pop == 'HG')
ids_HG <- unique(gtypes_long_HG$ID)
samples <- unique(gtypes_long_HG$sample)
binaries_HG <- samples %>% 
  map_dfc(~ ifelse(ids_HG %in% filter(gtypes_long_HG, sample == .x)$ID, 1, 0) %>% 
            as.data.frame)
names(binaries_HG) <- samples
library(UpSetR)
upset(binaries_HG, nsets = length(samples), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")
colSums(binaries_HG)
#48% dups are shared between atleast 2 samples
# HG_HG0F2 HG_HG2F1 HG_HG2M5 
# 914     1010      776 
# For NG
gtypes_long_NG <- anti_join(gtypes_long, filter_dups) %>% filter(num_alts > 0) %>% 
  filter(pop == 'NG')
ids_NG <- unique(gtypes_long_NG$ID)
samples_NG <- unique(gtypes_long_NG$sample)
binaries_NG <- samples_NG %>% 
  map_dfc(~ ifelse(ids_NG %in% filter(gtypes_long_NG, sample == .x)$ID, 1, 0) %>% 
            as.data.frame)
names(binaries_NG) <- samples_NG
#library(UpSetR)
upset(binaries_NG, nsets = length(samples_NG), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")
# 75% dups are shared between atleast 2 samples
colSums(binaries_NG)
# NG_NH0H4 NG_NH2F6 NG_NH2F8 NG_NH2M1 
# 2342     2095     2137     2242 

## Frequency of duplications per chromosome POST FILTERATION ##
gtypes_pos_fil <- map_dfr(select(oysterdup3_fil,CL_1:UMFS_6),getg)
gtypes_pos_fil$POS <- oysterdup3_fil$POS
gtypes_pos_long_fil <- gather(gtypes_pos_fil,key=sample,value=gtype,-POS)
gtypes_pos_long_fil$pop <- str_split(gtypes_pos_long_fil$sample,'_') %>% map(1) %>% as.character()
gtypes_pos_long_fil$pop <- as.vector(gtypes_pos_long_fil$pop)
gtypes_pos_long_fil$num_alts <- str_split(gtypes_pos_long_fil$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
#adding dups in all individuals of same pop to give pop count
pop_num_pos_alts_fil <- gtypes_pos_long_fil %>% filter(!is.na(num_alts)) %>%
  group_by(pop,POS) %>%
  summarize(num_alts = sum(num_alts))
pop_num_pos_alts_present_fil <- filter(pop_num_pos_alts_fil,num_alts >0)
chrom_pos_fil <- oysterdup3_fil %>% select(CHROM, POS) 
pop_num_pos_alts_present_chrom_fil <- left_join(pop_num_pos_alts_present_fil, chrom_pos, by = "POS")
pop_alts_per_chrom_fil <- pop_num_pos_alts_present_chrom_fil %>% group_by(pop,CHROM) %>% 
  summarize(num_alts = sum(num_alts)) 
#defining order of pop levels for plotting 
pop_alts_per_chrom_fil$pop <- factor (as.character(pop_alts_per_chrom_fil$pop), 
                                  levels=c("HI","SM","CS","HC","HCVA","CLP","CL","SL","LM","UMFS","NEH","HG","NG","DEBY","LOLA","OBOYS2"))
ggplot(pop_alts_per_chrom_fil, aes(x=CHROM,y=num_alts, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Frequency of CNVs", title ="Post filteration") + scale_color_manual(values=values,labels=labels)
# normalized by chromosome size
chrom_len <- data.frame(CHROM=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                                "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
                        start=c(1,1,1,1,1,1,1,1,1,1), 
                        end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045))
chrom_len$len <- chrom_len$end - chrom_len$start
pop_alts_per_chrom_len_fil <- left_join(pop_alts_per_chrom_fil, chrom_len, by = "CHROM")
pop_alts_per_chrom_len_fil$pop <- factor (as.character(pop_alts_per_chrom_len_fil$pop), 
                                      levels=c("HI","SM","CS","HC","HCVA","CLP","CL","SL","LM","UMFS","NEH","HG","NG","DEBY","LOLA","OBOYS2"))
ggplot(pop_alts_per_chrom_len_fil, aes(x=CHROM,y=(num_alts/len), color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Frequency of CNVs",title = "Post filteration")+ scale_color_manual(values=values,labels=labels)

##Copy number analysis POST FILTERATION ##
cn_gtypes_long_fil <- anti_join(cn_gtypes_long, filter_dups)
#what does copy nymber look like for the filtered dups themselves 
semi_join(cn_gtypes_long, filter_dups) %>% View()
#check number of dups to verify correct filteration
length(unique(cn_gtypes_long_fil[["ID"]])) #11339
#cn stats POST FILTERATION unchanged after filteration# 
min(cn_gtypes_long_fil[,5], na.rm=T) #-1
max(cn_gtypes_long_fil[,5], na.rm=T) #20025
hist(cn_gtypes_long_fil$cn)
#change the sample names and order the levels 
#cn on chr1 POST FILTERATION #
cn_gtypes_long_chr1_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035780.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr1_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr1_fil$cn))
#put in sample names as you want them to appear in the plot and in the order you want#
cn_gtypes_long_chr1_fil$sample <- factor (as.character(cn_gtypes_long_chr1_fil$sample), 
                                      levels=c("CL_1", "CL_2","CL_3","CL_4","CL_5","CL_6",
                                               "CLP_1","CLP_2","CLP_3","CLP_4","CLP_5","CLP_6",
                                               "CS_1","CS_2","CS_3","CS_5","CS_6","CS_7",
                                               "DEBY_1","DEBY_2","DEBY_3","DEBY_4","DEBY_5","DEBY_6",
                                               "HC_1","HC_3","HC_4","HC_5","HC_6","HC_7",
                                               "HCVA_1","HCVA_2","HCVA_3","HCVA_4","HCVA_5","HCVA_6",
                                               "HG_HG0F2","HG_HG2F1","HG_HG2M5",
                                               "HI_1","HI_2","HI_3","HI_4","HI_5", "HI_6",
                                               "LM_1_pool", "LM_3","LM_4","LM_7","LM_8",
                                               "LOLA_1","LOLA_2","LOLA_3","LOLA_4","LOLA_5","LOLA_6",
                                               "NEH_1","NEH_2","NEH_3","NEH_4","NEH_5","NEH_6",
                                               "NG_NH0H4", "NG_NH2F6","NG_NH2F8","NG_NH2M1", 
                                               "OBOYS2_1","OBOYS2_2","OBOYS2_3","OBOYS2_4","OBOYS2_5","OBOYS2_6",
                                               "SL_1", "SL_2", "SL_3", "SL_4", "SL_5","SL_6",
                                               "SM_10","SM_11","SM_12","SM_7","SM_8","SM_9",
                                               "UMFS_1","UMFS_2","UMFS_3","UMFS_4","UMFS_5","UMFS_6"))
cn_chr1_hmap_fil <- ggplot(data = cn_gtypes_long_chr1_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 1") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10)) #+ scale_fill_manual(labels=sample_lables)
cn_chr1_hmap_fil
# cn_chr1_hmap + geom_point(data=cn_gtypes_long_chr1, aes(x=65668439/2, y=1), col="red",pch=24, cex = 3)
cn_chr1_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3) + scale_fill_manual(labels=sample_lables)

cn_gtypes_long_chr2_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035781.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr2_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr2_fil$cn))
cn_chr2_hmap_fil <- ggplot(data = cn_gtypes_long_chr2_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 2") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr2_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr3_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035782.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr3_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr3_fil$cn))
cn_chr3_hmap_fil <- ggplot(data = cn_gtypes_long_chr3_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 3") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr3_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr4_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035783.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr4_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr4_fil$cn))
cn_chr4_hmap_fil <- ggplot(data = cn_gtypes_long_chr4_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 4") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr4_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr5_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035784.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr5_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr5_fil$cn))
cn_chr5_hmap_fil <- ggplot(data = cn_gtypes_long_chr5_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 5") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr5_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr6_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035785.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr6_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr6_fil$cn))
cn_chr6_hmap_fil <- ggplot(data = cn_gtypes_long_chr6_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 6") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr6_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr7_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035786.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr7_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr7_fil$cn))
cn_chr7_hmap_fil <- ggplot(data = cn_gtypes_long_chr7_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 7") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr7_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr8_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035787.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr8_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr8_fil$cn))
cn_chr8_hmap_fil <- ggplot(data = cn_gtypes_long_chr8_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 8") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr8_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr9_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035788.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr9_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr9_fil$cn))
cn_chr9_hmap_fil <- ggplot(data = cn_gtypes_long_chr9_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 9") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr9_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

cn_gtypes_long_chr10_fil <- filter(cn_gtypes_long_fil, CHROM == "NC_035789.1") %>% select(POS, sample, cn) 
cn_gtypes_long_chr10_fil$cn <- as.numeric(as.character(cn_gtypes_long_chr10_fil$cn))
cn_chr10_hmap_fil <- ggplot(data = cn_gtypes_long_chr10_fil, mapping = aes(x = POS,y = sample,color = log(cn))) + 
  geom_point(aes(cex=cn/100)) + xlab(label = "Position")+ggtitle(label = "Chr 10") +scale_color_viridis_c(direction = -1, na.value = "#f6f7a4",limits = c(0, 10))
cn_chr10_hmap_fil + geom_vline(xintercept = (65668439/2), color = "red", size=0.3)

##Dups mapped to different genome features ###
#Overlap between filtered dups and ALL genome features was obtained by making a BED file 
#from the GFF3 file for all features (just pulling out the 4 relevant columns)
# and using bedtools intersect to get overlaps with dups. 
dup_genome_overlap <- read.table("/Users/tejashree/Documents/Projects/cnv//scripts/output_files/bluewaves/dup_genome_feat_overlap_mod.bed", 
                              sep="\t" , stringsAsFactors = FALSE)
colnames(dup_fam_overlap) <- c("CHROM", "POS","end","ID","G_POS","G_end","Feature_name","l")
#Instead of all genome features I pulled out exons and genes only 
#then after intersect with dups just do wc -l on the output file to get
#count of dups completely within exons and within genes


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
ifi44 %>% group_by(ID) %>% tally() %>% ggplot(aes(x=ID, y=n)) + geom_col(fill='darkblue') + ylab(label = "Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ifi44_sub <- select(ifi44,ID,sample,pop)
cn <- select(cn_gtypes_long,ID,sample,pop,cn)
cn$cn <- as.numeric(as.character(cn$cn))
left_join(ifi44_sub,cn) %>% ggplot(aes(pop,cn)) + geom_col()
left_join(ifi44_sub,cn) %>% ggplot(aes(cn,sample, color =ID)) + geom_jitter() + xlim(c(0,15))#removing outlier 
left_join(ifi44_sub,cn) %>% ggplot(aes(pop,cn,color=ID)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Populations", y="Copy Number")
left_join(ifi44_sub,cn) %>% ggplot(aes(cn, sample,color=ID)) + geom_bar(stat = "identity", fill="white") + 
  labs(y="Samples", x="Copy Number")
#Inter and intra population copy number variation
left_join(ifi44_sub,cn) %>% ggplot(aes(cn,pop, color =pop)) + facet_wrap(~ID) + geom_jitter() + xlim(c(0,15))#removing outlier 
left_join(ifi44_sub,cn) %>% ggplot(aes(cn,ID, color =ID)) + facet_wrap(~pop) + geom_jitter() + xlim(c(0,15))#removing outlier 
#Getting the sequences of all DUP_IDs mapped as IFI44
#Getting ready a bedfile for bedtools to get the fasta
left_join(ifi44_ID,oysterdup3) %>% select(CHROM, POS, end, ID) %>%
  write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/ifi44_dup.bed", append = FALSE, sep = "\t",quote = FALSE,
              row.names = F, col.names = FALSE)
#The DUPs are on 3 diff chromosomes 1,8 and 9. The ref genome has 33 genes (unique LOCs annotated as IFI44) and are present on the same 3 chromosomes
#Dups and exons mapped on ref genome using script density_plot.R shows that some dups are overlapping, some span the exon some dont, 
#dups are present in 3 of the 4 regions where IFI44 exons are mapped. 

#GIMAP genes
dplyr::filter(dup_annot, grepl('GTPase IMAP family member', annot)) %>% select('ID') %>% unique() %>% tally() #23 multiple members 4,7,8
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 4', annot)) %>% select('ID') %>% unique() %>% tally() #21
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 7', annot)) %>% select('ID') %>% unique() %>% tally() #9
dplyr::filter(dup_annot, grepl('GTPase IMAP family member 8', annot)) %>% select('ID') %>% unique() %>% tally() #5
#The individual tally doesnt add up because for some DUPs they are mapped to multiple LOCs
dplyr::filter(dup_annot, grepl('scavenger receptor', annot)) %>% select('ID') %>% unique() %>% tally() #37 (multiple types/classes)
#now find out how many of these dups are present in each population
#performing same analysis on GIMAP dups as IFI44
gimap_ID <- dplyr::filter(dup_annot, grepl('GTPase IMAP family member', annot)) %>% select('ID') %>% unique()
gimap <- left_join(gimap_ID, gtypes_long_pres ,by = "ID") %>% filter(pres == 'yes')
#plotting freq of the dup IDs that map to GIMAP across all INDIVIDUALS not separated by populations.
gimap %>% group_by(ID) %>% tally() %>% ggplot(aes(x=ID, y=n)) + geom_col(fill='darkblue') + ylab(label = "Number of samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gimap_sub <- select(gimap,ID,sample,pop)
left_join(gimap_sub,cn) %>% ggplot(aes(pop,cn)) + geom_col()
#Inter and intra population copy number variation
left_join(gimap_sub,cn) %>% ggplot(aes(cn,pop, color =pop)) + facet_wrap(~ID) + geom_jitter()  
left_join(gimap_sub,cn) %>% ggplot(aes(cn,ID, color =ID)) + facet_wrap(~pop) + geom_jitter() 
labels = c("CL"="LA-HSCL","CLP"="CB-LSCP","CS"="DB-HSCS","DEBY"="CB-DEBY","HC"="DB-LSHC",
           "HCVA"="CB-HSHC","HG"="DB-NEHG","HI"="ME-HSHI","LM"="TX-HSLM","LOLA"="CB-LOLA",
           "NEH"="DB-NEHD","NG"="DB-NEHS","OBOYS2"="LA-OBOY","SL"="LA-LSSL","SM"="ME-LSSM","UMFS"="ME-UMFS")
values = c("CL"="#a63603","CLP"="#a6d854","CS"="#08519c","DEBY"="#006d2c","HC"="skyblue",
           "HCVA"="#31a354","HG"="#6a51a3","HI"="#fbb4b9","LM"="black","LOLA"="#bae4b3",
           "NEH"="#9e9ac8","NG"="#dadaeb","OBOYS2"="#fd8d3c","SL"="#fdae6b","SM"="#7a0177","UMFS"="#f768a1")
shapes = c("CL"=17,"CLP"=17,"CS"=17,"DEBY"=16,"HC"=17,
           "HCVA"=17,"HG"=13,"HI"=17,"LM"=17,"LOLA"=16,
           "NEH"=16,"NG"=13,"OBOYS2"=16,"SL"=17,"SM"=17,"UMFS"=16)
#Replotting for consistent pop names,colors and shapes.
left_join(gimap_sub,cn) %>% ggplot(aes(cn,pop, color =pop, shape=pop, label=pop)) + facet_wrap(~ID) + geom_jitter() + 
  scale_color_manual(values=values,labels=labels) + 
  scale_shape_manual(values=shapes,labels=labels) + scale_y_discrete(labels=labels)
#Getting the sequences of all DUP_IDs mapped as GIMAP
#Getting ready a bedfile for bedtools to get the fasta
left_join(gimap_ID,oysterdup3) %>% select(CHROM, POS, end, ID) %>% 
  write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/gimap_dup.bed", append = FALSE, sep = "\t",quote = FALSE,
              row.names = F, col.names = FALSE)
#GIMAP dups present on 4 chromosomes:2,6,7,8
#To get fasta of the dups scp gimap_dup.bed to bluewaves and use 
#bedtools getfasta -fi /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/cvir_edited.fa -bed gimap_dup.bed -fo gimap_dup.fa -name
#scp file back to cnv_output_files
#cat /Users/tejashree/Documents/Projects/cnv/annot/gimap_exon.bed | cut -f1 | uniq -c will show the distribution of exons on each chromosome
#gimap genes are present on 7 different chromosomes! #2,4,5,6,7,8,9
#The ref genome has 55 genes (unique LOCs annotated as GTPase IMAP family member)
#In order to get to the phylogeny of GIMAP genes we need to first get their sequences
# I get the LOCs mapped to each GIMAP gene from ref_annot
gimap_genes <- dplyr::filter(ref_annot, grepl('GTPase IMAP family member', annot)) %>% select('LOC') %>% unique()
oyster_genes <- read.table("/Users/tejashree/Documents/Projects/cnv/annot/Oyster_gene.bed",stringsAsFactors = FALSE)
colnames(oyster_genes)  <- c("CHROM", "POS","end", "LOC")
gimap_genes_bed <- semi_join(oyster_genes,gimap_genes)
gimap_genes_bed2 <- semi_join(oyster_genes,gimap_genes) 
gimap_genes_bed2$len <- gimap_genes_bed2$end - gimap_genes_bed2$POS
gimap_genes_bed %>% write.table("/Users/tejashree/Documents/Projects/cnv/scripts/output_files/oyster_cnv/gimap_genes.bed", append = FALSE, sep = "\t",quote = FALSE,
            row.names = F, col.names = FALSE)
#scp file to bluewaves tejashree/oyster_cnv/ and used bedtools to get fasta gimap_genes.fa 
## POST FILTERATION## 
# 2 dups mapped to gimap need to be filtered due to filteration criteria. 
#DUP00223590, DUP01190157

## pulling out dups mapped to histone genes
dplyr::filter(dup_annot, grepl('histone', annot)) %>% left_join(cn_gtypes_long,by="ID") %>% View()


### Duplications and expanded family mapping ###
# Intersect between filtered dups and expanded families (as defined by CAFE analysis) were obtained using bedtools
#Skipping those with 0 overlap
dup_fam_overlap <- read.table("/Users/tejashree/Documents/Projects/cnv//scripts/output_files/bluewaves/dup_cv_expanded_unique_overlap_mod.bed", 
                                 sep="\t" , stringsAsFactors = FALSE, skip = 7)
colnames(dup_fam_overlap) <- c("CHROM", "POS","end","ID","F_POS","F_end","Sequence_name","l")
#Number of dups overlapping expanded families by at least 23bp
#There are 545 duplications overlapping with 669 LOCs and 1254 protein IDs annotated to 880 unique annotations (usually diff isoforms or transcript variants). 
#So there are multiple mappings per duplication.
#Out of these only 274 are annotated as something other than 'uncharacterized'
dup_fam_overlap %>% select("ID","l") %>% group_by(ID) %>% tally() %>% nrow() #545
dup_fam_overlap %>% select("ID","l") %>% group_by(ID) %>% tally() %>% View()
dup_fam_overlap %>% select("Sequence_name") %>% group_by(Sequence_name) %>% tally() %>% nrow() #1254 
dup_fam_overlap %>% select("Sequence_name") %>% group_by(Sequence_name) %>% tally() %>% View()
#Cannot get the % overlap like I did with the repeats 
#because there are multiple mapping of DUP IDs that cannot be aggregated because they map
#with different XP IDs but same LOC IDs
#get annotation for families that overlap with dups
dup_fam_overlap %>% left_join(ref_annot_prot)%>% left_join(ref_annot) %>% View()
dup_fam_overlap %>% left_join(ref_annot_prot)%>% left_join(ref_annot) %>% group_by(LOC) %>% 
  filter(!is.na(LOC)) %>% tally() %>% nrow() #669
dup_fam_overlap %>% left_join(ref_annot_prot)%>% left_join(ref_annot) %>% filter(!is.na(LOC)) %>% 
  filter(!is.na(annot)) %>% 
  filter(!grepl("uncharacterized",annot)) %>% group_by(LOC) %>% tally() %>% nrow() #274
#Interesting finds: complement C1q-like protein 4,complement C1q tumor necrosis factor-related protein 3,complement C1q subcomponent subunit B-like
#interferon-induced protein 44-like,interleukin-17 receptor D,
#mucin-1-like, mucin-2-like, mucin-5AC-like and other isoforms, mucin-3A-like,mucin-17-like, natural resistance-associated macrophage protein 2
# toll-like receptor 3, toll-like receptor 13, lectin BRA-3,hepatic lectin, fucolectin-4-like, scavenger receptor class F member 1-like, protease inhibitors-like,metalloproteinase inhibitor 3-like
#nuclear apoptosis-inducing factor 1-like,death domain-containing protein CRADD-like
dup_fam_overlap %>% left_join(ref_annot_prot)%>% left_join(ref_annot) %>% filter(!is.na(LOC)) %>% filter(!is.na(annot)) %>% 
  filter(!grepl("uncharacterized",annot))%>% select("annot") %>% distinct() %>% View()
#Why dont any dups map to GIMAP members with this analysis?
left_join(ref_annot, ref_annot_prot) %>% filter(grepl('GTPase IMAP family member',annot)) %>% 
  select('Sequence_name') %>% distinct() %>% semi_join(tmp2) #GIMAPs seem to be absent from the CAFE output both shared and unique!
#Once this mystery is solved I can make a table of how many duplications mapped to these expanded genes of interest.
#Also have to think about how much overlap should be used as a cutoff  

###Vst calculations###
# Vpopx is the CN variance for each respective population
# Npopx is the number of individuals in each population
# Ntotal is the total sample size
# Vtotal is total variance
#Above is calculated for each 500bp window on each chromosome






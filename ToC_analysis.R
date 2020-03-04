## Analysis for theories of change paper

##Load required libraries & functions
library(tidyr)
library(stringr)
library(Amelia)
library(reshape2)
library(pscl)
library(aod)
library(car)
library(dplyr)
library(Hmisc)
library(corrplot)
library(polycor)
library(foreign)
library(nnet)
library(ggplot2)
library(stargazer)
library(mlogit)
library(xtable)

##define functions
assignBiomeGroups <- function(input,output){
  rows <- c(1:nrow(input))
  biome_groups <- matrix(nrow=nrow(input),ncol=1)
  rownames(biome_groups) <- rows
  colnames(biome_groups) <- c("biome_group")
  
  for (i in rows){
    biome <- input$Biome.[i]
    if (is.na(biome)) {
      group <- "NA"
    } else if (biome == "M_P" | biome == "M_TSS" | biome == "M_TU" | biome == "M_TRU" | biome == "M_TRC" | biome == "M_TSTSS") {
      group <- "MAR"
    } else if (biome == "FW_LL" | biome == "FW_LRD" | biome == "FW_PF" | biome == "FW_MF" | biome == "FW_TCR" | biome == "FW_TFRW" | biome == "FW_TUR" | biome == "FW_TSTCR" | biome == "FW_TSTFRW" | biome == "FW_TSTUR" | biome == "FW_XFEB" | biome == "FW_OI") {
      group <- "FRW"
    } else if (biome == "T_TSTMBF" | biome == "T_TSTDBF" | biome == "T_TSTCF" | biome == "T_TBMF" | biome == "T_TCF" | biome == "T_BFT" | biome == "T_MFWS") {
      group <- "FOR"
    } else if (biome == "T_TSTGSS" | biome == "T_TGSS" | biome == "T_FGS" | biome == "T_MGS") {
      group <- "GRS"
    } else if (biome == "T_T") {
      group <- "TUN"
    } else if (biome == "T_DXS") {
      group <- "DES"
    } else if (biome == "T_M") {
      group <- "MAN"
    } 
    biome_groups[i,"biome_group"] <- group
  }
  biome_groups <- as.data.frame(biome_groups)
  output <- bind_cols(input,biome_groups)
  output <- filter(output,!is.na(biome_groups))
  return(output) 
}

setwd("~/Documents/Manuscripts/SNAPP_ToC_Paper/ToC/")
data <- readRDS("~/Documents/github/McKinnon_et_al_2016/map_data_final.rds")
responses <- read.csv("Responses.csv",header=TRUE) %>% distinct()

## replace NAs in region
data$region <- as.character(data$region)
data$region[is.na(data$region)] <- "Global"

## create column with major habitat types
data <- assignBiomeGroups(data)

##Fix typos and NAs
data$Pub_year[data$Pub_year == "205"] <- "2015"
data$Pub_year[data$Pub_year == "N/A"] <- NA

##Create dummy variables for major habitat type, region, assessor
mht <- select(data,aid,biome_group) %>% distinct()
for(level in unique(mht$biome_group)){
  mht[paste("mht",level,sep="_")] <- ifelse(mht$biome_group == level, 1, 0)
}
mht <- select(mht,-biome_group)

mht <- mht %>% group_by(aid) %>% summarise_all(funs(sum))

region <- select(data,aid,region) %>% distinct()
for(level in unique(region$region)){
  region[paste("region",level,sep="_")] <- ifelse(region$region == level, 1, 0)
}
region <- select(region,-region)

region <- region %>% group_by(aid) %>% summarise_all(funs(sum))
colnames(region) <- c("aid","region_Africa","region_Asia","region_LatinAmerica","region_Oceania","region_Europe","region_Global")

assessor <- select(data.biblio,aid,Assessor) %>% distinct()
assessor$Assessor[grep("MCM ",assessor$Assessor)] <- "MCM"
assessor$Assessor[grep("mcm",assessor$Assessor)] <- "MCM"
assessor$Assessor[grep("MCB",assessor$Assessor)] <- "MCM"

for(level in unique(assessor$Assessor)){
  assessor[paste("Assessor",level,sep="_")] <- ifelse(assessor$Assessor == level, 1, 0)
}
assessor <- select(assessor,-Assessor)

## Convert publication year to continuous scale starting at 1
data$Pub_year <- as.numeric(data$Pub_year)
min <- min(data$Pub_year)
max <- max(data$Pub_year)
range <- seq(min,max,by=1)
len <- c(1:length(range))
year <- as.data.frame(range)
year <- cbind(year,as.data.frame(len))
colnames(year) <- c("year","year_mod")

pub_year <- select(data,aid,Pub_year)
v <- c(1:nrow(pub_year))
pub_year$year_mod <- c("")
for (i in v){
  y <- pub_year[i,"Pub_year"]
  n <- filter(year, year == y)
  pub_year[i,"year_mod"] <- n[1,"year_mod"]
}

pub_year <- select(pub_year,-Pub_year)
colnames(pub_year) <- c("aid","Pub_year")

pub_year <- distinct(pub_year)

## Calculate number of interventions and outcomes
n <- unique(unlist(data$aid))
number <- matrix(ncol=3,nrow=length(n))
colnames(number) <- c("aid","no_int","no_out")
dat <- dplyr::select(data,aid,Int_type,Outcome) %>% distinct()

for (i in 1:length(n)){
  sub <- filter(dat, aid == n[i])
  number[i,"aid"] <- n[i]
  number[i,"no_int"] <- n_distinct(sub$Int_type)
  number[i,"no_out"] <- n_distinct(sub$Outcome)
}

number <- as.data.frame(number)

## Create dummy variable for interventions and outcomes
number$no_int.f <- cut(number$no_int, breaks=c(0,1,100), include.lowest=TRUE,labels=c("1","GreaterThan_1"))

for(level in unique(number$no_int.f)){
  number[paste("No_int",level,sep="_")] <- ifelse(number$no_int.f == level, 1, 0)
}

number$no_out.f <- cut(number$no_out, breaks=c(0,1,100), include.lowest=TRUE,labels=c("1","GreaterThan_1"))

for(level in unique(number$no_out.f)){
  number[paste("No_out",level,sep="_")] <- ifelse(number$no_out.f == level, 1, 0)
}

number <- select(number,-no_out.f,-no_int.f)

#Code variables for impact evaluations

ie <- select(data,aid,IE) %>% distinct()
ie[ie=="Y"] <- 1
ie[ie=="N"] <- 0

##Set up analysis frames
all_data <- responses %>% left_join(assessor,by="aid") %>% distinct() %>% left_join(region,by="aid") %>% distinct() %>% left_join(pub_year,by="aid") %>% distinct() %>% left_join(number,by="aid") %>% distinct() %>% left_join(mht,by="aid") %>% distinct() %>% left_join(ie,by="aid") %>% distinct() %>% select(-X,-Authors,-Title,-Pub_year.y)

sapply(all_data, function(x) is.factor(x)) ##Check for factors
cols <- c(colnames(mht),colnames(region),colnames(assessor),"No_int_GreaterThan_1","No_int_1","No_out_GreaterThan_1","No_out_1","Concept_mod","Graphical","Robust","IE")
cols <- cols[!grepl("aid",cols)]
all_data[cols] <- lapply(all_data[cols],factor)
all_data$Pub_year.x <- as.numeric(all_data$Pub_year.x)

all_data$Response2 <- as.character(all_data$Response)
all_data$Response2[all_data$Response2 == "None"] <- 1
all_data$Response2[all_data$Response2 == "Non-graphical"] <- 2
all_data$Response2[all_data$Response2 == "Graphical non-robust"] <- 3
all_data$Response2[all_data$Response2 == "Graphical robust"] <- 0
all_data$Response2 <- as.numeric(all_data$Response2)
all_data$Response2 <- factor(all_data$Response2, levels=c(0,1,2,3))

all_data <- select(all_data,-aid)

continuous_vars <- select(all_data,Pub_year.x,no_int,no_out)
categorical_vars <- select(all_data,-Pub_year.x,-no_int,-no_out,-Concept_mod,-Graphical,-Robust,-Response,-Response2)
categorical_vars[] <- lapply(categorical_vars, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

all_vars <- select(all_data,-Concept_mod,-Graphical,-Robust,-Response,-No_int_GreaterThan_1,-No_out_GreaterThan_1,-No_int_1, -No_out_1,-Response2)
all_vars[] <- lapply(all_vars, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

##Check missing data and number of unique values
sapply(all_data,function(x) sum(is.na(x)))
sapply(all_data, function(x) length(unique(x)))

##Fix missing data
pdf("Missing data.pdf",height=8.5,width=11)
missmap(all_data, main = "Missing values vs. observed")
dev.off()

##Correlation matrices
continuous_cor <- cor(continuous_vars,method="pearson")
continuous_cor_p <- rcorr(as.matrix(continuous_vars))
corrplot(continuous_cor_p$r, type="upper",order="hclust",p.mat=continuous_cor_p$P, sig.level=0.05, insig="blank")

categorical_cor <- polychor(categorical_vars)
categorical_cor_p <- rcorr(as.matrix(categorical_vars))
heatmap(x=as.matrix(categorical_cor_p), col=c("blue","white","red"),symm=TRUE)

test <- cor.ci(categorical_vars,n.iter=10,poly=TRUE)

all_cor <- hetcor(as.data.frame(all_vars))
all_cor_p <- rcorr(as.matrix(all_vars))
pdf("Correlation_matrix_all.pdf",height=8.5,width=11)
corrplot(all_cor_p$r, type="upper", order="hclust", p.mat=all_cor_p$P, sig.level=0.05, insig="blank")
dev.off()
all_cor <- as.data.frame(all_cor)


##Analysis - binomial logistic regression model
##Response variable - conceptual model or not
resp_1 <- select(all_data,-Graphical,-Robust,-Response,-Response2,-No_int_GreaterThan_1,-No_out_GreaterThan_1,-No_int_1, -No_out_1)

##Response variable - non model, graphical non-robust model, non-graphical model, robust graphical model
resp_2 <- select(all_data,-Graphical,-Robust,-Response2,-Concept_mod,-No_int_GreaterThan_1,-No_out_GreaterThan_1,-No_int_1, -No_out_1)

resp_2 <- select(all_data,-Graphical,-Robust,-Response2,-Concept_mod,-No_int_GreaterThan_1,-No_out_GreaterThan_1,-No_int_1, -No_out_1,-Assessor_SHC, -Assessor_MCM,-Assessor_JE, -Assessor_ES,-Assessor_SZS)


##Response variable - robust vs. not robust ...
resp_3 <- filter(all_data, Graphical == 1)
resp_3 <- select(resp_3, -Graphical, -Concept_mod, -No_int_GreaterThan_1,-No_out_GreaterThan_1,-No_int_1, -No_out_1,-Response,-Response2)

##Check missing data and number of unique values
sapply(resp_1,function(x) sum(is.na(x)))
sapply(resp_1, function(x) length(unique(x)))

resp_1 <- select(resp_1,-Assessor_ES,-Assessor_SZS,-region_Europe,-mht_TUN)

#Run binomial logistic regression model
##Set 1
model <- glm(Concept_mod ~.,family=binomial,data=resp_1,maxit=5000) #run bionomial logit
summary(model)
anova(model,test="Chisq") #Analyze deviance
pR2(model) ##calculate pR2 model fit
OR_CI <- exp(cbind(OR = coef(model), confint(model)))
model.rrr <- exp(coef(model))
model.rrr2 <- coef(model)

##Calculate p-value of model with predictors versus null model
with(model, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE)) 

stargazer(model, type="html",coef=list(model.rrr),p.auto=FALSE,out="binom_oddsratio_4_1_2018.htm")
stargazer(model, type="html",coef=list(model.rrr2),p.auto=FALSE,out="binom_coef_4_1_2018.htm")

##Set 3 
##t-test
t.test(resp_3$Pub_year.x~resp_3$Robust)
t.test(resp_3$no_int~resp_3$Robust)
t.test(resp_3$no_out~resp_3$Robust)

##chi-square
library(MASS)
tbl = table(resp_3$IE,resp_3$Robust)
chisq.test(tbl)

##Summary statistic graphs
years <- c(min(all_data$Pub_year.x):max(all_data$Pub_year.x))
year_count <- matrix(nrow=26,ncol=2)
rownames(year_count) <-years
colnames(year_count) <- c("Robust","Not robust")
for (i in years){
  sub <- filter(resp_3, Pub_year.x == i)
  r <- count(sub, Robust)
  year_count[i, "Robust"] <- as.numeric(r[2,2])
  year_count[i, "Not robust"] <- as.numeric(r[1,2])
}
year_count[is.na(year_count)] <- 0
Pub_year <- as.data.frame(rownames(year_count))
year_count <- as.data.frame(year_count)
year_count <- bind_cols(year_count,Pub_year)
colnames(year_count) <- c("Robust","Not_robust","Publication_year")
year_count2 <- melt(year_count, id.vars="Publication_year",measure.vars=c("Robust","Not_robust"))

p1 <- ggplot(year_count2,aes(x=Publication_year, y=value, fill=variable)) + 
      geom_bar(stat="identity", position="dodge") +
      labs(x="PUBLICATION YEAR", y="NUMBER OF ARTICLES") +
  scale_x_discrete(breaks=c(1:26),labels=c(1990:2015),limits=c(1:26)) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  theme(legend.title=element_blank())

int_out_count <- select(resp_3,no_int,no_out,Robust)
int_out_count2 <- melt(int_out_count)

p2 <- ggplot(int_out_count2, aes(x=variable, y=value, fill=Robust)) +
  geom_boxplot() +
  theme(axis.title.x = element_blank(), legend.title=element_blank(), legend.position="none") +
  ylab("COUNT") +
  scale_x_discrete(breaks=c("no_int","no_out"),labels=c("INTERVENTIONS","OUTCOMES")) +
  scale_fill_discrete(breaks=c(0,1),labels=c("Non-robust","Robust")) +
  expand_limits(y=0)

mht_count <- select(resp_3,mht_GRS:mht_NA,Robust)
mht_count[,1:7] <- lapply(mht_count[,1:7], function(x) as.numeric(as.character(x)))
mht_count2 <- mht_count %>% group_by(Robust) %>% summarise_each(funs(sum))
mht_count2 <- melt(mht_count2, id.vars="Robust",measure.vars=c("mht_GRS","mht_FRW","mht_FOR","mht_MAR","mht_DES","mht_MAN","mht_NA"))

p3 <- ggplot(mht_count2, aes(x=variable, y=value, fill=Robust)) +
  geom_bar(stat="identity",position="dodge") +
  theme(axis.title.x = element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle=45,hjust=1),legend.position="none") +
  ylab("NUMBER OF ARTICLES") +
  scale_x_discrete(breaks=c("mht_GRS","mht_FRW","mht_FOR","mht_MAR","mht_MAN","mht_DES","mht_NA"),labels=c("Grasslands","Freshwater","Forests","Marine","Mangroves","Deserts","Not specified")) +
  scale_fill_discrete(breaks=c(0,1),labels=c("Non-robust","Robust"))

reg_count <- select(resp_3,region_Africa:region_Global,Robust)
reg_count[,1:5] <- lapply(reg_count[,1:5], function(x) as.numeric(as.character(x)))
reg_count2 <- reg_count %>% group_by(Robust) %>% summarise_each(funs(sum))
reg_count2 <- melt(reg_count2, id.vars="Robust",measure.vars=c("region_Africa","region_Asia","region_LatinAmerica","region_Oceania","region_Global"))

p4 <- ggplot(reg_count2, aes(x=variable, y=value, fill=Robust)) +
  geom_bar(stat="identity", position="dodge") +
  theme(axis.title.x = element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle=45,hjust=1), legend.position="none") +
  ylab("NUMBER OF ARTICLES") +
  scale_x_discrete(breaks=c("region_Africa","region_Asia","region_LatinAmerica","region_Oceania","region_Global"),labels=c("Africa","Asia","Latin America","Oceania","Global")) +
  scale_fill_discrete(breaks=c(0,1),labels=c("Non-robust","Robust"))

ie_count <- select(resp_3,IE,Robust)
ie_count[,1] <- lapply(ie_count[,1], function(x) as.numeric(as.character(x)))
ie_count2 <- ie_count %>% group_by(Robust) %>% count(IE)
ie_count2$IE[ie_count2$IE == 1] <- "IE"
ie_count2$IE[ie_count2$IE == 0] <- "not_IE"

p5 <- ggplot(ie_count2, aes(x=IE, y=n, fill=Robust)) +
  geom_bar(stat="identity", position="dodge") +
  theme(axis.title.x = element_blank(), legend.title=element_blank(), axis.text.x = element_text(angle=45,hjust=1), legend.position="none") +
  ylab("NUMBER OF ARTICLES") +
  scale_x_discrete(breaks=c("IE","not_IE"),labels=c("Impact evaluation","Not impact evaluation"))
  scale_fill_discrete(breaks=c(0,1),labels=c("Non-robust","Robust"))

library(gridExtra)
library(cowplot)

  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")

pdf("Robust_v_not_Robust.pdf", height=8.5,width=11)

grid.arrange(p1,p2,p3,p4,p5, ncol=3,nrow=2,
             layout_matrix = cbind(c(1,3),c(1,4),c(2,5)))

dev.off()

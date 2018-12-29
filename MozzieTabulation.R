library(ggplot2)
library(readr)

setwd(wd <- "~/Projects/Malaria collab/Exported results/")   # Working directory name
filepath <- paste(wd, "Mozzie mosquito all plates readonly.csv", sep="")    # CSV file path
dat <- read.csv(filepath, header=TRUE, stringsAsFactors=FALSE)   # Raw data # read_csv?
dat <- dat[-which(dat$Sample.Name==""), ]    # Remove rows where sample name is blank
dat[dat=="Undetermined"] <- NA
dat$HbtubCT1 <- as.numeric(dat$HbtubCT1)
dat$HbtubCT2 <- as.numeric(dat$HbtubCT2)
dat$pfr364CT1 <- as.numeric(dat$pfr364CT1)
dat$pfr364CT2 <- as.numeric(dat$pfr364CT2)
dat[dat <= 0] <- NA   # Change non-positive values to NAs
dat$HbtubMeanCT <- rowMeans(dat[, c("HbtubCT1","HbtubCT2")], na.rm=TRUE)
dat$pfr364MeanCT <- rowMeans(dat[, c("pfr364CT1","pfr364CT2")], na.rm=TRUE)
dat[is.na(dat)] <- NA   # Change NaNs to NAs

hbtub_dat <- dat[which(dat$HbtubCT1 > 0 | dat$HbtubCT2 > 0), ]    # Data with positive Hbtub CT
pfr_dat <- dat[which(dat$pfr364CT1 > 0 | dat$pfr364CT2 > 0), ]    # Data with positive pfr364 CT
ids <- names(table(dat$Identifier))    # All unique mosquito identifiers
both_ids <- dat$Identifier[duplicated(dat$Identifier)]    # IDs that have both H and A

tab <- matrix(rep(NA, 90), nrow=10, ncol=9,
        dimnames=list(c("in.H.villageK","in.H.villageM","in.H.villageS","in.H.all",
                        "in.A.villageK","in.A.villageM","in.A.villageS","in.A.all",
                        "in.either","in.both"),
                      c("Total","Hbtub.Positive","Hbtub.Prevalence","Hbtub.Mean","Hbtub.SD",
                        "pfr364.Positive","pfr364.Prevalence","pfr364.Mean","pfr364.SD")))    # Tabulated data

tabulate <- function(type, region) {
  in_type <- paste("in", type, sep=".")
  if(type %in% c("H","A")) {
    if(missing(region)) {
      in_type <- paste(in_type, "all", sep=".")
      tab[in_type, "Total"] <<- sum(dat$Specimen.Type==type)
      hbtub_type <- hbtub_dat$Specimen.Type==type
      pfr_type <- pfr_dat$Specimen.Type==type
    } else {
      in_type <- paste(in_type, region, sep=".village")
      tab[in_type, "Total"] <<- sum(dat$Specimen.Type==type & dat$Region==region)
      hbtub_type <- hbtub_dat$Specimen.Type==type & hbtub_dat$Region==region
      pfr_type <- pfr_dat$Specimen.Type==type & pfr_dat$Region==region
    }
  } else if(type=="either") {    # Mosquitoes with either H or A or both
    tab[in_type, "Total"] <<- length(ids)
    hbtub_type <- ids %in% hbtub_dat$Identifier
    pfr_type <- ids %in% pfr_dat$Identifier
  } else if(type=="both") {    # Mosquitoes with both H and A
    tab[in_type, "Total"] <<- length(both_ids)
    hbtub_type <- both_ids %in% hbtub_dat$Identifier
    pfr_type <- both_ids %in% pfr_dat$Identifier
  }
  tab[in_type, "Hbtub.Positive"] <<- sum(hbtub_type)
  tab[in_type, "Hbtub.Prevalence"] <<- tab[in_type, "Hbtub.Positive"] / tab[in_type, "Total"]
  tab[in_type, "Hbtub.Mean"] <<- mean(hbtub_dat$HbtubMeanCT[which(hbtub_type)], na.rm=TRUE)
  tab[in_type, "Hbtub.SD"]   <<-   sd(hbtub_dat$HbtubMeanCT[which(hbtub_type)], na.rm=TRUE)
  tab[in_type, "pfr364.Positive"] <<- sum(pfr_type)
  tab[in_type, "pfr364.Prevalence"] <<- tab[in_type, "pfr364.Positive"] / tab[in_type, "Total"]
  tab[in_type, "pfr364.Mean"] <<- mean(pfr_dat$pfr364MeanCT[which(pfr_type)], na.rm=TRUE)
  tab[in_type, "pfr364.SD"]   <<-   sd(pfr_dat$pfr364MeanCT[which(pfr_type)], na.rm=TRUE)
}

for(type in c("H", "A")) {
  tabulate(type, "K")
  tabulate(type, "M")
  tabulate(type, "S")
  tabulate(type)
}
tabulate("either")
tabulate("both")

pfr_ct_dat <- as.data.frame(sapply(dat[, c("pfr364CT1","pfr364CT2")], as.numeric))    # Only pfr364 CTs
pfr_ct_dat <- pfr_ct_dat[complete.cases(pfr_ct_dat), ]    # Remove NAs
pfr_ct_dat <- pfr_ct_dat[which(pfr_ct_dat$pfr364CT1 > 1e-6 & pfr_ct_dat$pfr364CT2 > 1e-6), ]    # Remove rows where CT == 0
print(ggplot(pfr_ct_dat, aes(x=pfr364CT1, y=pfr364CT2))) + 
  geom_point() + 
  geom_smooth(method=lm)    # Compare the replicate pfr364 CT values for each mosquito sample

write.csv(tab, file = "tabulation.csv")    # Export tabulated data
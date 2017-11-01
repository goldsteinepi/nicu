#################
# NICU spatial analysis
# Citation: Goldstein ND, Tuttle D, Tabb LP, Paul DA, Eppes SC. Spatial and Environmental Correlates of Organism Colonization and Infection in the Neonatal Intensive Care Unit. Manuscript in preparation.
# 10/13/16 -- Neal Goldstein
#################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(lme4) #nonlinear mixed effects (glmer)
library(blme) #bayesian nonlinear mixed effects (bglmer)
library(maptools) #mapping
library(RColorBrewer) #color palette
library(spdep) #spatial modeling
library(boot) #bootstrapping
library(RODBC) #connect to SQL Server

#returns the MOR from a glmer model, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
bootMOR = function(model)
{
  return(exp(0.95*getME(model,"theta")))
}


### READ DATA ###

#ET cultures: 2006-2015 intubated infants with at least one ETT culture performed
load("et_data_sepsis.RData")
ET_cultures = NICU
rm(Vent4plus,Vent6plus,Vent8plus)

#NICU data (NOTE: must use the same source data as for et_data)
load("NICU.2016-11-18.RData")

#read NICU map, see http://www.goldsteinepi.com/blog/creatingyourownspatialanalysismapsasshapefileswithqgis
NICU_map = readShapePoly("Shapefiles/NICU", proj4string=CRS("+proj=tmerc +lat_0=38 +lon_0=-75.41666666666667 +k=0.999995 +x_0=200000.0001016002 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs"))

#MRN and DOB are added from external data but redacted from this code


### SUBSET and CREATE COHORT ###

#no birthweight means not admitted to NICU
NICU = NICU[!is.na(NICU$Birthweight), ]

#limit data to 2006 to 2015 corresponding to SafeySurveillor data
NICU = NICU[NICU$Admission_year>=2006 & NICU$Admission_year<=2015, ]

#limit to single admissions only
NICU = NICU[NICU$Admission_n==1, ]


### ADD DATA and RECODE ###

#late onset sepsis
NICU$Sepsis_late = ifelse(!is.na(NICU$Sepsis_onset) & NICU$Sepsis_onset==1, 1, 0)

#MRSA colonization
set.seed(777)
NICU$MRSA_colonization_09 = ifelse((NICU$MRSA_colonization==0) & (runif(nrow(NICU), min=0, max=1)<rep((0.009 - sum(NICU$MRSA_colonization)/nrow(NICU)),nrow(NICU))), 1, NICU$MRSA_colonization)
NICU$MRSA_colonization_15 = ifelse((NICU$MRSA_colonization==0) & (runif(nrow(NICU), min=0, max=1)<rep((0.015 - sum(NICU$MRSA_colonization)/nrow(NICU)),nrow(NICU))), 1, NICU$MRSA_colonization)
NICU$MRSA_colonization_22 = ifelse((NICU$MRSA_colonization==0) & (runif(nrow(NICU), min=0, max=1)<rep((0.022 - sum(NICU$MRSA_colonization)/nrow(NICU)),nrow(NICU))), 1, NICU$MRSA_colonization)

#mean center
NICU$Admission_year_centered = scale(NICU$Admission_year, center=T, scale=F)
NICU$Gestational_age_centered = scale(NICU$Gestational_age, center=T, scale=F)
NICU$Birthweight_centered = scale(NICU$Birthweight, center=T, scale=F)
NICU$LOS_centered = scale(NICU$Birthweight, center=T, scale=F)
NICU$Census_average_centered = scale(NICU$Census_average, center=T, scale=F)

#bed location to mapping location
NICU$Location_map = ifelse(NICU$Location=="201A" | NICU$Location=="201B" | NICU$Location=="201C", "201", NICU$Location)
NICU$Location_map = ifelse(NICU$Location=="202A" | NICU$Location=="202B" | NICU$Location=="202C", "202", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="203A" | NICU$Location=="203B" | NICU$Location=="203C", "203", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="204A" | NICU$Location=="204B" | NICU$Location=="204C", "204", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="205A" | NICU$Location=="205B" | NICU$Location=="205C", "205", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="206A" | NICU$Location=="206B" | NICU$Location=="206C", "206", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="207A" | NICU$Location=="207B" | NICU$Location=="207C", "207", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="208C", NA, NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="209A" | NICU$Location=="209B" | NICU$Location=="209C", "209", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="211A" | NICU$Location=="211B" | NICU$Location=="211C", "211", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="212A" | NICU$Location=="212B" | NICU$Location=="212C", "212", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="213A" | NICU$Location=="213B" | NICU$Location=="213C", "213", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="214A" | NICU$Location=="214B" | NICU$Location=="214C", "214", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="215A" | NICU$Location=="215B" | NICU$Location=="215C", "215", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="216C", NA, NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="217A" | NICU$Location=="217B" | NICU$Location=="217C", "217", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="218A" | NICU$Location=="218B" | NICU$Location=="218C", "218", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="219A" | NICU$Location=="219B" | NICU$Location=="219C", "219", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="220A" | NICU$Location=="220B" | NICU$Location=="220C", "220", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="221A" | NICU$Location=="221B" | NICU$Location=="221C", "221", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="222A" | NICU$Location=="222B" | NICU$Location=="222C", "222", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="223A" | NICU$Location=="223B" | NICU$Location=="223C", "223", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="224A" | NICU$Location=="224B" | NICU$Location=="224C", "224", NICU$Location_map)
NICU$Location_map = ifelse(NICU$Location=="225A" | NICU$Location=="225B" | NICU$Location=="225C", "225", NICU$Location_map)

#bed location to map location, culture results
ET_cultures$Location_map = ifelse(ET_cultures$Location=="201A" | ET_cultures$Location=="201B" | ET_cultures$Location=="201C", "201", ET_cultures$Location)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="202A" | ET_cultures$Location=="202B" | ET_cultures$Location=="202C", "202", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="203A" | ET_cultures$Location=="203B" | ET_cultures$Location=="203C", "203", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="204A" | ET_cultures$Location=="204B" | ET_cultures$Location=="204C", "204", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="205A" | ET_cultures$Location=="205B" | ET_cultures$Location=="205C", "205", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="206A" | ET_cultures$Location=="206B" | ET_cultures$Location=="206C", "206", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="207A" | ET_cultures$Location=="207B" | ET_cultures$Location=="207C", "207", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="208C", NA, ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="209A" | ET_cultures$Location=="209B" | ET_cultures$Location=="209C", "209", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="211A" | ET_cultures$Location=="211B" | ET_cultures$Location=="211C", "211", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="212A" | ET_cultures$Location=="212B" | ET_cultures$Location=="212C", "212", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="213A" | ET_cultures$Location=="213B" | ET_cultures$Location=="213C", "213", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="214A" | ET_cultures$Location=="214B" | ET_cultures$Location=="214C", "214", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="215A" | ET_cultures$Location=="215B" | ET_cultures$Location=="215C", "215", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="216C", NA, ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="217A" | ET_cultures$Location=="217B" | ET_cultures$Location=="217C", "217", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="218A" | ET_cultures$Location=="218B" | ET_cultures$Location=="218C", "218", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="219A" | ET_cultures$Location=="219B" | ET_cultures$Location=="219C", "219", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="220A" | ET_cultures$Location=="220B" | ET_cultures$Location=="220C", "220", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="221A" | ET_cultures$Location=="221B" | ET_cultures$Location=="221C", "221", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="222A" | ET_cultures$Location=="222B" | ET_cultures$Location=="222C", "222", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="223A" | ET_cultures$Location=="223B" | ET_cultures$Location=="223C", "223", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="224A" | ET_cultures$Location=="224B" | ET_cultures$Location=="224C", "224", ET_cultures$Location_map)
ET_cultures$Location_map = ifelse(ET_cultures$Location=="225A" | ET_cultures$Location=="225B" | ET_cultures$Location=="225C", "225", ET_cultures$Location_map)

#add outcome and covariate data
NICU$Vent_start_day = NA
NICU$Vent_culture = NA
NICU$Pseudomonas_colonization = NA
NICU$Klebsiella_colonization = NA
NICU$Location_map_single = NA
NICU$Location_equipment = NA
NICU$Gloving = NA
NICU$Location_ceiling_tile = NA
NICU$Location_carpet = NA
NICU$Location_refrigerator = NA
NICU$Infection_policy = NA

#add MRN to NICU data
NICU$MRN = MRNlist$MRN[which(MRNlist$ID %in% NICU$ID)]

#add DOB to NICU data
NICU$DOB = as.Date(DOBlist$DateOfBirth[which(DOBlist$ID %in% NICU$MRN)])

for (i in 1:nrow(NICU))
{
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #vent cultures performed
  
  if (NICU$Vent[i]==1) {
    #day of week vent started
    NICU$Vent_start_day[i] = weekdays(NICU$DOB[i]+NICU$Vent_start[i])
    
    #check if intubated on a Monday (when respiratory cultures are taken)
    if (NICU$Vent_start_day[i]=="Sunday" && NICU$Vent_length[i]>=1) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Monday" && NICU$Vent_length[i]>=7) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Tuesday" && NICU$Vent_length[i]>=6) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Wednesday" && NICU$Vent_length[i]>=5) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Thursday" && NICU$Vent_length[i]>=4) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Friday" && NICU$Vent_length[i]>=3) {
      NICU$Vent_culture[i] = 1
    } else if (NICU$Vent_start_day[i]=="Saturday" && NICU$Vent_length[i]>=2) {
      NICU$Vent_culture[i] = 1
    } else {
      NICU$Vent_culture[i] = 0
    }
  }
  
  #pseudomonas
  NICU$Pseudomonas_colonization[i] = ifelse(length(grep("Pseudomonas", ET_cultures$Vent_culture_organisms[ET_cultures$ID==NICU$ID[i]]))>0, 1, ifelse(NICU$Vent_culture[i]==1, 0, NA))
  
  #klebseilla
  NICU$Klebsiella_colonization[i] = ifelse(length(grep("Klebsiella", ET_cultures$Vent_culture_organisms[ET_cultures$ID==NICU$ID[i]]))>0, 1, ifelse(NICU$Vent_culture[i]==1, 0, NA))
  
  #single patient rooms
  NICU$Location_map_single[i] = ifelse(NICU$Location_map[i]=="208A" | NICU$Location_map[i]=="208B" | NICU$Location_map[i]=="210A" | NICU$Location_map[i]=="210B" | NICU$Location_map[i]=="217" | NICU$Location_map[i]=="216A" | NICU$Location_map[i]=="216B" | NICU$Location_map[i]=="220", 1, 0)
  
  #patient rooms where more equipment can fit
  NICU$Location_equipment[i] = ifelse(NICU$Location_map[i]=="208A" | NICU$Location_map[i]=="209" | NICU$Location_map[i]=="217" | NICU$Location_map[i]=="216A" | NICU$Location_map[i]=="216B" | NICU$Location_map[i]=="220", 1, 0)
  
  #universal gloving: based on dates from Deb Tuttle
  NICU$Gloving[i] = ifelse((NICU$Date_admission[i]>=as.Date("2008-08-18") & NICU$Date_discharge_initial[i]<=as.Date("2010-02-28")) | (NICU$Date_admission[i]>=as.Date("2012-10-01") & NICU$Date_discharge_initial[i]<=as.Date("2013-10-15")), 1, 0)
  
  #leaking ceiling tile/pseudomonas
  NICU$Location_ceiling_tile[i] = ifelse(NICU$Date_admission[i]>=as.Date("2006-05-01") & NICU$Date_admission[i]<=as.Date("2006-05-10") & NICU$Location_map[i]=="207", 1, 0)
  
  #carpet in quiet room
  NICU$Location_carpet[i] = ifelse(NICU$Date_admission[i]<=as.Date("2007-09-15") & (NICU$Location_map[i]=="222" | NICU$Location_map[i]=="220" | NICU$Location_map[i]=="218" | NICU$Location_map[i]=="223" | NICU$Location_map[i]=="221" | NICU$Location_map[i]=="219"), 1, 0)
  
  #refrigerator
  NICU$Location_refrigerator[i] = ifelse(NICU$Location_map[i]=="208B" | NICU$Location_map[i]=="216B" | NICU$Location_map[i]=="220", 1, 0)
  
  #unit specific infection policy
  NICU$Infection_policy[i] = ifelse(NICU$Date_admission[i]>=as.Date("2013-03-01"), 1, 0)
}
rm(i,MRNlist,DOBlist)
NICU$MRN = NULL
NICU$DOB = NULL
rownames(NICU) = NULL

#indicators and outcomes in the ecological data
map_data=data.frame("bed"=as.character(slot(NICU_map, "data")$bed), "MRSA_colonization"=NA, "RSV"=NA, "Sepsis_late"=NA, "Pseudomonas_colonization"=NA, "Klebsiella_colonization"=NA, "Location_equipment"=NA, stringsAsFactors=F)

for (i in 1:nrow(map_data))
{
  map_data$MRSA_colonization[i] = sum(NICU$MRSA_colonization[NICU$Location_map==map_data$bed[i]], na.rm=T)
  map_data$RSV[i] = sum(NICU$RSV[NICU$Location_map==map_data$bed[i]], na.rm=T)
  map_data$Sepsis_late[i] = sum(NICU$Sepsis_late[NICU$Location_map==map_data$bed[i]], na.rm=T)
  map_data$Pseudomonas_colonization[i] = length(grep("Pseudomonas", ET_cultures$Vent_culture_organisms[ET_cultures$Location_map==map_data$bed[i]]))
  map_data$Klebsiella_colonization[i] = length(grep("Klebsiella", ET_cultures$Vent_culture_organisms[ET_cultures$Location_map==map_data$bed[i]]))
  map_data$Location_equipment[i] = ifelse(map_data$bed[i]=="208A" | map_data$bed[i]=="209" | map_data$bed[i]=="217" | map_data$bed[i]=="216A" | map_data$bed[i]=="216B" | map_data$bed[i]=="220", 1, 0)
}
rm(i,ET_cultures)

#merge the indicator dataframe to the shapefile, note we do not sort data to preserve the order of the polygons 
#NICU_map@data = merge(x=NICU_map@data, y=map_data, by="bed", all.x=T, sort=F)
NICU_map = merge(x=NICU_map, y=map_data, by="bed", all.x=T)
rm(map_data)

#save data
save.image("spatial.RData")
#load("spatial.RData")


### ANALYSIS INDIVIDUAL LEVEL VARS ###

CrossTable(NICU$MRSA_colonization)
#CrossTable(NICU$RSV)
CrossTable(NICU$Sepsis_late)

CrossTable(NICU$Vent)
CrossTable(NICU$Vent_culture)
CrossTable(NICU$Pseudomonas_colonization)
CrossTable(NICU$Klebsiella_colonization)

describe(NICU$Gestational_age); IQR(NICU$Gestational_age,na.rm=T)
describe(NICU$Birthweight)
CrossTable(NICU$Central_line)
CrossTable(NICU$Antibiotic)
CrossTable(NICU$PROM)
CrossTable(NICU$Chorio_clinical)
CrossTable(NICU$Mom_antibiotic)
CrossTable(NICU$Outborn)
describe(NICU$Census_average)

describe(NICU$Gestational_age[NICU$Vent_culture==1]); IQR(NICU$Gestational_age[NICU$Vent_culture==1],na.rm=T)
describe(NICU$Birthweight[NICU$Vent_culture==1])
CrossTable(NICU$Central_line[NICU$Vent_culture==1])
CrossTable(NICU$Antibiotic[NICU$Vent_culture==1])
CrossTable(NICU$PROM[NICU$Vent_culture==1])
CrossTable(NICU$Chorio_clinical[NICU$Vent_culture==1])
CrossTable(NICU$Mom_antibiotic[NICU$Vent_culture==1])
CrossTable(NICU$Outborn[NICU$Vent_culture==1])
describe(NICU$Census_average[NICU$Vent_culture==1])
describe(NICU$Vent_length[NICU$Vent_culture==1]); IQR(NICU$Vent_length[NICU$Vent_culture==1],na.rm=T)

describeBy(NICU$Admission_year, NICU$MRSA_colonization); t.test(NICU$Admission_year ~ NICU$MRSA_colonization)
describeBy(NICU$Gestational_age, NICU$MRSA_colonization); t.test(NICU$Gestational_age ~ NICU$MRSA_colonization)
describeBy(NICU$Birthweight, NICU$MRSA_colonization); wilcox.test(NICU$Birthweight ~ NICU$MRSA_colonization)
IQR(NICU$Birthweight); IQR(NICU$Birthweight[NICU$MRSA_colonization==0]); IQR(NICU$Birthweight[NICU$MRSA_colonization==1])
CrossTable(NICU$Central_line, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Care_team, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(NICU$Vent_length, NICU$MRSA_colonization); wilcox.test(NICU$Vent_length ~ NICU$MRSA_colonization)
#IQR(NICU$Vent_length); IQR(NICU$Vent_length[NICU$MRSA_colonization==0]); IQR(NICU$Vent_length[NICU$MRSA_colonization==1])
CrossTable(NICU$Antibiotic, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Mom_antibiotic, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Census_average, NICU$MRSA_colonization); wilcox.test(NICU$Census_average ~ NICU$MRSA_colonization)

# CrossTable(NICU$Admission_year, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# describeBy(NICU$Gestational_age, NICU$RSV); t.test(NICU$Gestational_age ~ NICU$RSV)
# describeBy(NICU$Birthweight, NICU$RSV); wilcox.test(NICU$Birthweight ~ NICU$RSV)
# IQR(NICU$Birthweight); IQR(NICU$Birthweight[NICU$RSV==0]); IQR(NICU$Birthweight[NICU$RSV==1])
# CrossTable(NICU$Central_line, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Care_team, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# #describeBy(NICU$Vent_length, NICU$RSV); wilcox.test(NICU$Vent_length ~ NICU$RSV)
# #IQR(NICU$Vent_length); IQR(NICU$Vent_length[NICU$RSV==0]); IQR(NICU$Vent_length[NICU$RSV==1])
# CrossTable(NICU$Antibiotic, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$PROM, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Chorio_clinical, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Mom_antibiotic, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

CrossTable(NICU$Admission_year, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Gestational_age, NICU$Sepsis_late); t.test(NICU$Gestational_age ~ NICU$Sepsis_late)
describeBy(NICU$Birthweight, NICU$Sepsis_late); wilcox.test(NICU$Birthweight ~ NICU$Sepsis_late)
IQR(NICU$Birthweight); IQR(NICU$Birthweight[NICU$Sepsis_late==0]); IQR(NICU$Birthweight[NICU$Sepsis_late==1])
CrossTable(NICU$Central_line, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Care_team, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(NICU$Vent_length, NICU$Sepsis_late); wilcox.test(NICU$Vent_length ~ NICU$Sepsis_late)
#IQR(NICU$Vent_length); IQR(NICU$Vent_length[NICU$Sepsis_late==0]); IQR(NICU$Vent_length[NICU$Sepsis_late==1])
CrossTable(NICU$Antibiotic, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Mom_antibiotic, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Census_average, NICU$Sepsis_late); wilcox.test(NICU$Census_average ~ NICU$Sepsis_late)

CrossTable(NICU$Admission_year, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Gestational_age, NICU$Pseudomonas_colonization); t.test(NICU$Gestational_age ~ NICU$Pseudomonas_colonization)
describeBy(NICU$Birthweight, NICU$Pseudomonas_colonization); wilcox.test(NICU$Birthweight ~ NICU$Pseudomonas_colonization)
IQR(NICU$Birthweight); IQR(NICU$Birthweight[NICU$Pseudomonas_colonization==0]); IQR(NICU$Birthweight[NICU$Pseudomonas_colonization==1])
CrossTable(NICU$Central_line, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Care_team, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(NICU$Vent_length, NICU$Pseudomonas_colonization); wilcox.test(NICU$Vent_length ~ NICU$Pseudomonas_colonization)
#IQR(NICU$Vent_length); IQR(NICU$Vent_length[NICU$Pseudomonas_colonization==0]); IQR(NICU$Vent_length[NICU$Pseudomonas_colonization==1])
CrossTable(NICU$Antibiotic, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Mom_antibiotic, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Census_average, NICU$Pseudomonas_colonization); wilcox.test(NICU$Census_average ~ NICU$Pseudomonas_colonization)

CrossTable(NICU$Admission_year, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Gestational_age, NICU$Klebsiella_colonization); t.test(NICU$Gestational_age ~ NICU$Klebsiella_colonization)
describeBy(NICU$Birthweight, NICU$Klebsiella_colonization); wilcox.test(NICU$Birthweight ~ NICU$Klebsiella_colonization)
IQR(NICU$Birthweight); IQR(NICU$Birthweight[NICU$Klebsiella_colonization==0]); IQR(NICU$Birthweight[NICU$Klebsiella_colonization==1])
CrossTable(NICU$Central_line, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Care_team, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#describeBy(NICU$Vent_length, NICU$Klebsiella_colonization); wilcox.test(NICU$Vent_length ~ NICU$Klebsiella_colonization)
#IQR(NICU$Vent_length); IQR(NICU$Vent_length[NICU$Klebsiella_colonization==0]); IQR(NICU$Vent_length[NICU$Klebsiella_colonization==1])
CrossTable(NICU$Antibiotic, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Mom_antibiotic, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Census_average, NICU$Klebsiella_colonization); wilcox.test(NICU$Census_average ~ NICU$Klebsiella_colonization)


### ANALYSIS ENVIRONMENTAL VARS ###

CrossTable(NICU$Location_map_single)
CrossTable(NICU$Location_equipment)
CrossTable(NICU$Location_refrigerator)

CrossTable(NICU$Location_map_single[NICU$Vent_culture==1])
CrossTable(NICU$Location_equipment[NICU$Vent_culture==1])
CrossTable(NICU$Location_refrigerator[NICU$Vent_culture==1])

CrossTable(NICU$Location_map_single, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(MRSA_colonization ~ (1 | Location_map) + Location_map_single, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_equipment, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(MRSA_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_refrigerator, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(MRSA_colonization ~ (1 | Location_map) + Location_refrigerator, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))

# CrossTable(NICU$Location_map_single, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(RSV ~ (1 | Location_map) + Location_map_single, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
# CrossTable(NICU$Location_equipment, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(RSV ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
# CrossTable(NICU$Location_refrigerator, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(RSV ~ (1 | Location_map) + Location_refrigerator, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))

CrossTable(NICU$Location_map_single, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Sepsis_late ~ (1 | Location_map) + Location_map_single, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_equipment, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Sepsis_late ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_refrigerator, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Sepsis_late ~ (1 | Location_map) + Location_refrigerator, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))

CrossTable(NICU$Location_map_single, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Pseudomonas_colonization ~ (1 | Location_map) + Location_map_single, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_equipment, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Pseudomonas_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_refrigerator, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Pseudomonas_colonization ~ (1 | Location_map) + Location_refrigerator, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))

CrossTable(NICU$Location_map_single, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Klebsiella_colonization ~ (1 | Location_map) + Location_map_single, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_equipment, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Klebsiella_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))
CrossTable(NICU$Location_refrigerator, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=F); summary(glmer(Klebsiella_colonization ~ (1 | Location_map) + Location_refrigerator, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa")))

#check correlation among group level vars
#CrossTable(NICU$Location_map_single, NICU$Location_equipment, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#CrossTable(NICU$Location_map_single, NICU$Location_refrigerator, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#CrossTable(NICU$Location_equipment, NICU$Location_refrigerator, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### ANALYSIS TIME DEPENDENT VARS ###

CrossTable(NICU$Location_ceiling_tile, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Location_carpet, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Gloving, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Infection_policy, NICU$MRSA_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

# CrossTable(NICU$Location_ceiling_tile, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Location_carpet, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Gloving, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
# CrossTable(NICU$Infection_policy, NICU$RSV, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

CrossTable(NICU$Location_ceiling_tile, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Location_carpet, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Gloving, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Infection_policy, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

CrossTable(NICU$Location_ceiling_tile, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Location_carpet, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Gloving, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Infection_policy, NICU$Pseudomonas_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

CrossTable(NICU$Location_ceiling_tile, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Location_carpet, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Gloving, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Infection_policy, NICU$Klebsiella_colonization, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)


### MIXED EFFECTS REGRESSION MODELING ###

#note if glmer fails to converge, try a different optimizer: https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html; http://www.ats.ucla.edu/stat/r/dae/melogit.htm
model = glmer(MRSA_colonization ~ (1 | Location_map), family=binomial(), data=NICU)
model = glmer(MRSA_colonization ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(MRSA_colonization ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="nloptwrap"))
summary(model)

#compare log likelihood values from models; as they are similar to the warning message, can treat message as a false positive
logLik(model)

## EMPTY MODELS

model = bglmer(MRSA_colonization ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
#model = glmer(RSV ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Sepsis_late ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Pseudomonas_colonization ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Klebsiella_colonization ~ (1 | Location_map), family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

summary(model)

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)

## ENVIRONMENTAL MODELS

#based on crude correlations, Location_equipment was the most important predictor
model = glmer(MRSA_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
#model = glmer(RSV ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Sepsis_late ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Pseudomonas_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Klebsiella_colonization ~ (1 | Location_map) + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

#OR and CI estimates for fixed effects
summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)

## INDIVIDUAL MODELS

model = glmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
#model = glmer(RSV ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Sepsis_late ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Pseudomonas_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Klebsiella_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

#OR and CI estimates for fixed effects
summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

## INDIVIDUAL+ENVIRONMENTAL MODELS

model = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
#model = bglmer(RSV ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Sepsis_late ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Pseudomonas_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(Klebsiella_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

#ETT colonization models have random effect variances estimated as zero, see this discussion: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-models-random-effect-variances-estimated-as-zero-or-correlations-estimated-as---1 
#using Bayesian linear mixed effect models to estimate random variance

#sensitivity analyses
model = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU[NICU$Admission_year>=2013,], control=glmerControl(optimizer="bobyqa"))
model = bglmer(MRSA_colonization_09 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(MRSA_colonization_15 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(MRSA_colonization_22 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU[NICU$Outborn==0,], control=glmerControl(optimizer="bobyqa"))
model = bglmer(Sepsis_late ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU[NICU$Outborn==0,], control=glmerControl(optimizer="bobyqa"))
model = bglmer(Pseudomonas_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU[NICU$Outborn==0,], control=glmerControl(optimizer="bobyqa"))
model = bglmer(Klebsiella_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU[NICU$Outborn==0,], control=glmerControl(optimizer="bobyqa"))

#OR and CI estimates for fixed effects
summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
set.seed(777)
NICU_boot = NICU
#NICU_boot = NICU[NICU$Outborn==0,]
#NICU_boot = NICU[NICU$Admission_year>=2013,]
#NICU_boot = NICU[!is.na(NICU$Vent_culture) & NICU$Vent_culture==1, ]
MOR = NA
for (i in 1:1000)
{
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  bootindex = sample(nrow(NICU_boot), nrow(NICU_boot), replace=T)
  bootdata = NICU_boot[bootindex,]
  #model = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  #model = bglmer(MRSA_colonization_09 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  model = bglmer(MRSA_colonization_15 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  #model = bglmer(MRSA_colonization_22 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  #model = bglmer(Sepsis_late ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  #model = bglmer(Pseudomonas_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  #model = bglmer(Klebsiella_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=bootdata, control=glmerControl(optimizer="bobyqa"))
  MOR = c(MOR, exp(0.95*getME(model,"theta")))
}
rm(i,NICU_boot,bootindex,bootdata)
MOR = MOR[-1]
quantile(MOR, probs=c(0.025, 0.5, 0.975))

#spatial autocorrelation
moran.mc(ranef(model)[[1]][[1]], listw=sa.wt, nsim=1000)

## INDIVIDUAL+ENVIRONMENTAL+TIME DEPENDENT MODELS

#note these models do not improve fit
model = glmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Infection_policy, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
#model = glmer(RSV ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Infection_policy, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Sepsis_late ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Infection_policy, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Pseudomonas_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Location_equipment + Infection_policy, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model = glmer(Klebsiella_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + log(Vent_length) + PROM + Chorio_clinical + Location_equipment + Infection_policy, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

#OR and CI estimates for fixed effects
summary(model)
exp(fixef(model))
exp(confint.merMod(model, method="Wald"))

#area level variance
getME(model,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model,"theta"))

#bootstrap precision around MOR
boot_model = bootMer(model, bootMOR, nsim=1000, parallel="multicore", ncpus=4)
boot.ci(boot_model, type="norm", index=1)


### SPATIAL MODELING ###

#output to tif for publication 
tiff("Figure 1d.tif",height=6,width=6,units='in',res=1200) 

#choropleth plots by 5-levels of shading
spplot(NICU_map, "MRSA_colonization", cuts=4, sp.layout=list("sp.text", coordinates(NICU_map), NICU_map$bed), col.regions=brewer.pal(5, "Reds"))
#spplot(NICU_map, "RSV", cuts=4, sp.layout=list("sp.text", coordinates(NICU_map), NICU_map$bed), col.regions=brewer.pal(5, "Reds"))
spplot(NICU_map, "Sepsis_late", cuts=4, sp.layout=list("sp.text", coordinates(NICU_map), NICU_map$bed), col.regions=brewer.pal(5, "Reds"))
spplot(NICU_map, "Pseudomonas_colonization", cuts=4, sp.layout=list("sp.text", coordinates(NICU_map), NICU_map$bed), col.regions=brewer.pal(5, "Reds"))
spplot(NICU_map, "Klebsiella_colonization", cuts=4, sp.layout=list("sp.text", coordinates(NICU_map), NICU_map$bed), col.regions=brewer.pal(5, "Reds"))

#close file 
dev.off() 

#define neighbors based on shared boundary or distance of 10 apart
sa.nb = poly2nb(NICU_map, queen=T, snap=10)
summary(sa.nb)

#adjust neighbors based on typical nursing assignments; see NICU_map@data for slot/bed list
sa.nb[[1]] = which(slot(NICU_map, "data")$bed %in% c("202","203"))
sa.nb[[2]] = which(slot(NICU_map, "data")$bed %in% c("201","205","202","204"))
sa.nb[[3]] = which(slot(NICU_map, "data")$bed %in% c("203","207","204","206"))
sa.nb[[4]] = which(slot(NICU_map, "data")$bed %in% c("205","208B","206","208A"))
sa.nb[[5]] = which(slot(NICU_map, "data")$bed %in% c("208B","208A"))
sa.nb[[6]] = which(slot(NICU_map, "data")$bed %in% c("207","209","208A"))
sa.nb[[7]] = which(slot(NICU_map, "data")$bed %in% c("209","208B","207","206"))
sa.nb[[8]] = which(slot(NICU_map, "data")$bed %in% c("208A","204","207","205"))
sa.nb[[9]] = which(slot(NICU_map, "data")$bed %in% c("202","206","205","203"))
sa.nb[[10]] = which(slot(NICU_map, "data")$bed %in% c("201","203","204"))
sa.nb[[11]] = which(slot(NICU_map, "data")$bed %in% c("210A","214","213","215"))
sa.nb[[12]] = which(slot(NICU_map, "data")$bed %in% c("212","215","217","216A","216B"))
sa.nb[[13]] = which(slot(NICU_map, "data")$bed %in% c("210B","210A","213A"))
sa.nb[[14]] = which(slot(NICU_map, "data")$bed %in% c("211","215","210A","212"))
sa.nb[[15]] = which(slot(NICU_map, "data")$bed %in% c("213","217","212","214"))
sa.nb[[16]] = which(slot(NICU_map, "data")$bed %in% c("215","214","216A","216B"))
sa.nb[[17]] = which(slot(NICU_map, "data")$bed %in% c("214","217","216A"))
sa.nb[[18]] = which(slot(NICU_map, "data")$bed %in% c("214","217","216B"))
sa.nb[[19]] = which(slot(NICU_map, "data")$bed %in% c("218","220","221"))
sa.nb[[20]] = which(slot(NICU_map, "data")$bed %in% c("223","219","218","220","222"))
sa.nb[[21]] = which(slot(NICU_map, "data")$bed %in% c("222","220","221"))
sa.nb[[22]] = which(slot(NICU_map, "data")$bed %in% c("220","221","223"))
sa.nb[[23]] = which(slot(NICU_map, "data")$bed %in% c("222","218","223","221","219"))
sa.nb[[24]] = which(slot(NICU_map, "data")$bed %in% c("220","221","219"))
sa.nb[[25]] = which(slot(NICU_map, "data")$bed %in% c("224"))
sa.nb[[26]] = which(slot(NICU_map, "data")$bed %in% c("225"))
sa.nb[[27]] = which(slot(NICU_map, "data")$bed %in% c("210B","212","211","213"))
sa.nb[[28]] = which(slot(NICU_map, "data")$bed %in% c("210A","211"))

#create a representation of the binary weights
sa.wt = nb2listw(neighbours=sa.nb, style="B")
summary(sa.wt)

#plot neighbor network
par(xpd=T, mar=c(1,2,1,2)) 
plot(NICU_map)
plot(sa.nb, coordinates(NICU_map), add=T, col="red")

#moran's spatial autocorrelation
moran.mc(NICU_map$MRSA_colonization, listw=sa.wt, nsim=1000)
#moran.mc(NICU_map$RSV, listw=sa.wt, nsim=1000)
moran.mc(NICU_map$Sepsis_late, listw=sa.wt, nsim=1000)
moran.mc(NICU_map$Pseudomonas_colonization, listw=sa.wt, nsim=1000)
moran.mc(NICU_map$Klebsiella_colonization, listw=sa.wt, nsim=1000)

#spatial correlogram, shows autocorrelation by order of neighbor (1=neighbor, 2=neighbor's neighbor, etc)
plot(sp.correlogram(neighbours=sa.nb,var=NICU_map$MRSA_colonization,order=4,method="I",style="B",zero.policy=T), main="")
#plot(sp.correlogram(neighbours=sa.nb,var=NICU_map$RSV,order=4,method="I",style="B",zero.policy=T), main="Spatial correlogram of RSV infection")
plot(sp.correlogram(neighbours=sa.nb,var=NICU_map$Sepsis_late,order=4,method="I",style="B",zero.policy=T), main="Spatial correlogram of Lab-confirmed late onset sepsis")
plot(sp.correlogram(neighbours=sa.nb,var=NICU_map$Pseudomonas_colonization,order=4,method="I",style="B",zero.policy=T), main="Spatial correlogram of Pseudomonas colonization")
plot(sp.correlogram(neighbours=sa.nb,var=NICU_map$Klebsiella_colonization,order=4,method="I",style="B",zero.policy=T), main="Spatial correlogram of Klebsiella colonization")

#spatial regression for outcomes that have spatial autocorrelation (ie MRSA)
summary(lagsarlm(MRSA_colonization ~ Location_equipment, data=NICU_map@data, listw=sa.wt, type="lag"))
summary(errorsarlm(MRSA_colonization ~ Location_equipment, data=NICU_map@data, listw=sa.wt))
summary(lagsarlm(MRSA_colonization ~ Location_equipment, data=NICU_map@data, listw=sa.wt, type="mixed"))

#select eigenvectors for multilevel modeling, see: https://www.ncbi.nlm.nih.gov/pubmed/24571639 
spatial_eigen = SpatialFiltering(MRSA_colonization ~ 1, ~ Location_equipment - 1, data=NICU_map@data, nb=sa.nb, style="W", ExactEV=T)
ncol(fitted(spatial_eigen)) #4 eigenvectors selected

#add to NICU dataset
NICU$Eigenvector1 = NA
NICU$Eigenvector2 = NA
NICU$Eigenvector3 = NA
NICU$Eigenvector4 = NA

for (i in 1:nrow(NICU))
{
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  NICU$Eigenvector1[i] = ifelse(!is.na(NICU$Location_map[i]), fitted(spatial_eigen)[which(slot(NICU_map, "data")$bed==NICU$Location_map[i]), 1], NA)
  NICU$Eigenvector2[i] = ifelse(!is.na(NICU$Location_map[i]), fitted(spatial_eigen)[which(slot(NICU_map, "data")$bed==NICU$Location_map[i]), 2], NA)
  NICU$Eigenvector3[i] = ifelse(!is.na(NICU$Location_map[i]), fitted(spatial_eigen)[which(slot(NICU_map, "data")$bed==NICU$Location_map[i]), 3], NA)
  NICU$Eigenvector4[i] = ifelse(!is.na(NICU$Location_map[i]), fitted(spatial_eigen)[which(slot(NICU_map, "data")$bed==NICU$Location_map[i]), 4], NA)
}
rm(i,spatial_eigen)

#mixed effects spatial model, with eigenvectors as random effects
model_nospatial = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model_spatial = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Census_average_centered + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))

#sensitivity analyses
model_spatial = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU[NICU$Admission_year>=2013,], control=glmerControl(optimizer="bobyqa"))
model_spatial = bglmer(MRSA_colonization_09 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model_spatial = bglmer(MRSA_colonization_15 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model_spatial = bglmer(MRSA_colonization_22 ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU, control=glmerControl(optimizer="bobyqa"))
model_spatial = bglmer(MRSA_colonization ~ (1 | Location_map) + Admission_year_centered + Gestational_age_centered + Central_line + Antibiotic + PROM + Chorio_clinical + Location_equipment + Eigenvector1 + Eigenvector2 + Eigenvector3 + Eigenvector4, family=binomial(), data=NICU[NICU$Outborn==0,], control=glmerControl(optimizer="bobyqa"))

#compare models
summary(model_nospatial)
summary(model_spatial)

#OR and CI estimates for fixed effects
summary(model_spatial)
exp(fixef(model_spatial))
exp(confint.merMod(model_spatial, method="Wald"))

#area level variance
getME(model_spatial,"theta")

#compute median OR, see: http://www.ncbi.nlm.nih.gov/pubmed/16537344
exp(0.95*getME(model_spatial,"theta"))

#moran's spatial autocorrelation using random intercepts from glmm models
moran.mc(ranef(model_nospatial)[[1]][[1]], listw=sa.wt, nsim=1000)
moran.mc(ranef(model_spatial)[[1]][[1]], listw=sa.wt, nsim=1000)

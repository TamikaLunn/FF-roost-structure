## Title: Code to structure roosting data for GAM model analysis on flying-fox roost structure in SE QLD and NE NSW
## Manuscript: Lunn et al (20201) Counterintuitive scaling between population abundance and local density: implications for modelling transmission of infectious diseases in bat populations
## Author: Tamika Lunn, Griffith University
## Version: Submission, created 12th November 2021

## V1-7 - Manuscript preparation
## VSubmission - Code transferred from V7 "Bat quantitative analysis_FF#2_V7-revision.R"

rm(list=ls())

##############################################################################################
##---------------------------------Load data & set functions--------------------------------##
##############################################################################################

### Load helper files and any additional packages
source ("Code/00_functions_VSubmission.R")
library(mgcv) #for GAMs
library(binom)

##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##--------------------------------------------------------------- Read in data ---------------------------------------------------------------##
## Tree survey data:
treesurvey <- read.csv("Data/Raw/tree survey data.csv")

## Extra centroid data
centroids <- read.csv("Data/Raw/Centroids.csv") %>%
  mutate(roost.centroid.N = ifelse(site.accession=="DLIS010",NA,roost.centroid.N)) %>% ## Remove May 2019 Lismore area point (DLIS010)
  mutate(roost.centroid.E = ifelse(site.accession=="DLIS010",NA,roost.centroid.E)) 
centroids$session <- as.factor(centroids$session)

## Bat structure data:
treebat <- read.csv("Data/Raw/spatial-bat-structure-data.csv", row.names=1)
treebat$subplot <- as.factor (treebat$subplot)
treebat$session <- as.factor (treebat$session)
head(treebat)[,c("site.code", "tree.accession", "site.accession", "BFF.index", "GHFF.index", "LRFF.index", "subplot", "crowngroup")]

## Pixel density data:
##with zero values:
kernelbat <- read.csv("Data/Raw/pixel-density-data.csv", row.names=1)
kernelbat$subplot <- as.factor (kernelbat$subplot)
kernelbat <- arrange(kernelbat, site.accession, subplot) #arrange makes customising kables easier
##without zero values:
kernelbatnz <- read.csv("Data/Raw/pixel-density-data-nonzero.csv", row.names=1)
kernelbatnz$subplot <- as.factor (kernelbatnz$subplot)
kernelbatnz <- arrange(kernelbatnz, site.accession, subplot) #arrange makes customising kables easier

## Roost-level data:
roostbat<-read.csv("Data/Raw/ALL-roost-use-data.csv", row.names=1) %>% filter(!is.na(site.code)) %>%
  mutate(roost.area = ifelse(site.accession=="DLIS010",NA,roost.area)) %>%
  mutate(roost.perimeter = ifelse(site.accession=="DLIS010",NA,roost.perimeter)) ## Remove May 2019 Lismore area point (DLIS010)
roostbat$session <- as.factor (roostbat$session)
head(roostbat)[,c("site.code", "session", "site.accession", "pop.estimate.BFF", "pop.estimate.GHFF", "pop.estimate.LRFF", "pop.estimate.total")]
temp.tree <- treebat %>%
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #removes missing trees and rare occasions where trees were missed. This is important to do so that proportion is calculated correctly
  plyr::ddply(c("site.code", "session", "site.accession"), summarise,
              total.trees = sum(!is.na(tree.accession)))
roostbat <- dplyr::left_join(roostbat, temp.tree, by=c("site.code", "session", "site.accession"))

## Additional plot data:
plotdata <- read.csv("Data/Raw/plot-data.csv", row.names=1) %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(subplot = as.factor(subplot)) %>%
  mutate(roost.type = as.factor(roost.type)) %>%
  mutate(plotID = as.factor(plotID))

## Tree tessellation data:
tree.tessellation <- read.csv("Data/Raw/Voronoi Polygons/tree.tessellation.csv", row.names = 1) %>%
  mutate(tree.accession = as.factor(tree.accession))  %>%
  mutate(site.code = as.factor(site.code))  %>%
  mutate(subplot = as.factor(subplot))  %>%
  mutate(crown = as.factor(crown))  %>%
  mutate(value_dirichlet = ifelse(crown=="I",Liqr(value_dirichlet), ## Set area of intermediate trees as the first quantile of the canopy trees (5.8 m2)
                                  ifelse(crown=="E",199, 
                                         ifelse(crown=="OG",199,
                                                value_dirichlet)))) %>%
  mutate(value_dirichlet = ifelse(value_dirichlet>199,199,value_dirichlet)) ## set maximum value, based on mean crown projection area for eucalypts in NSW, reported in Verma (2014) An allometric model for estimating DBH of isolated and clustered Eucalyptus trees from measurements of crown projection area

nrow(tree.tessellation[which(tree.tessellation$value_dirichlet==199),]) ##218 trees were over-rided with the max crown value
nrow(tree.tessellation)
##-------------------------------------------------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------- Format data ---------------------------------------------------------------##
## There are two main data types for bat data: 'index-based' and 'count-based'. Index-based data were taken for all trees (and include index and weighted index values), whereas count was taken for a subset of trees (and include count, min and max height).
## 'pixel-based' data are kernel density estimates, calculated from the spatial point pattern of tree distributions within plots weighted by the 'index-based' data per tree. 
## There are two possible summary levels: plot-level and roost-level. Note - in the published paper, 'plot' is referred to as 'subplot', and 'site' as 'roost'.
## There are two possible calculations for roost-level summaries: with tree as the lowest replicate ('bytree'), and plot as the lowest replicate ('byplot')
## Summary data frames are returned in either a wide format (species in one column, then separate columns for each summary value) OR long format (species in one column, then one column for 'value' and one column for 'measure'). May also be an 'extra-wide' format which reflects the way the data was entered (separate columns for site*measure combinations)
## Wide format for when named columns are to be selected for use. Long format for when values need to be filtered from a single column

### Naming convention for functions:
###      Functions: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format'
###      Output: 'data type' '.' 'LOWEST REPLICATE' '.' 'returned df format' '_' 'summary level'
###      data type = index | count
###      LOWEST REPLICATE = TREE | PLOT
###      returned df format = wide |Ewide | long
###      summary level = site | plot

###############################################
##### Plot-level (bat structure) overview #####
###############################################

### Calculate data summaries:
byplot <- c("site.code", "session", "site.accession", "subplot", "rep") ## Group by plot
bytree <- c("site.code", "session", "site.accession", "tree.accession", "id", "crowngroup", "subplot", "rep") ## Group by tree
list.count <- treebat  %>% #get list of accessions with tree count 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session", "subplot", "rep"), summarise,
        tree.count = sum(!is.na(tree.accession)))  
list.kernel <- treebat  %>% #get list of accessions with tree count
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session", "subplot"), summarise,
        tree.count = sum(!is.na(tree.accession)))  

## Calculate summaries & transpose dataframes into different formats
index.TREE.Ewide_plot <- index.TREE.Ewide(treebat, byplot) ## Dataframe with separate column per measure*species
index.TREE.wide_plot <- index.TREE.wide(index.TREE.Ewide_plot, byplot, list.count) ## Dataframe with separate column per measure
index.TREE.long_plot <- index.TREE.long(index.TREE.Ewide_plot, byplot) 
heights.subset <- calculate.heightdiff(treebat, bytree)  %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% 
  mutate(vert.density = count/(height.diff*value_dirichlet)) %>% 
  mutate(vert.density = ifelse(vert.density==Inf,NA,vert.density))
count.TREE.wide_plot <- count.TREE.wide(heights.subset, byplot, list.count) 
count.TREE.long_plot <- count.TREE.long(count.TREE.wide_plot, byplot)
kernel.TREE.wide_plot<- kernel.TREE.wide(kernelbat, list.kernel)
kernel.TREE.long_plot<- kernel.TREE.long(kernelbat, list.kernel)
kernelnz.TREE.wide_plot<- kernel.TREE.wide(kernelbatnz, list.kernel)
kernelnz.TREE.long_plot<- kernel.TREE.long(kernelbatnz, list.kernel)

###############################################
##### Site-level (bat structure) overview #####
###############################################

### Calculate data summaries. 
## Note - There are two calculations for site summaries: with tree as the lowest replicate, and plot as the lowest replicate. When plot is the lowest replicate, averages are calculated across occupied plots only. Mean in this case is the mean of the plot means. Standard deviation here is the standard deviation between plots. 

##---------------------------------------- Plot as lowest replicate ---------------------------------------- 
byplot <- c("site.code", "session", "site.accession", "subplot", "rep")
bysite <- c("site.code", "session", "site.accession")

## Create list of accessions and measures for merging:
list.count_site <- treebat  %>% #get list of accessions with tree count, doesn't matter which dataframe 
  filter(!is.na(BFF.index)&!is.na(GHFF.index)&!is.na(LRFF.index)) %>% #remove rare cases where trees missed in survey OR were removed
  ddply(c("site.accession", "site.code", "session"), summarise,
        tree.count = sum(!is.na(tree.accession)))  
list.index_site <- index.TREE.long_plot  %>%
  distinct(site.accession, site.code, session, measure) %>%
  mutate(mean.value = as.numeric(0)) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session")) ##Left join to get tree count per site
list.count_site <-index.TREE.long_plot  %>%
  distinct(site.accession, site.code, session) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session"))
list.pixel_site <- kernel.TREE.long_plot  %>%
  distinct(site.accession, site.code, session, measure) %>%
  mutate(mean.value = as.numeric(0)) %>%
  left_join(list.count_site, by=c("site.accession", "site.code", "session"))

##identify occupied plots:
occ.all.output_plot <- occ.all(treebat, byplot)
occ.BFF.output_plot <- occ.BFF(treebat, byplot)
occ.GHFF.output_plot <- occ.GHFF(treebat, byplot)
occ.LRFF.output_plot <- occ.LRFF(treebat, byplot)

##calculate:
index.PLOT.long_BFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "BFF", occ.BFF.output_plot, list.index_site) #by plot, because you want to select occupied plots 
index.PLOT.long_GHFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "GHFF", occ.GHFF.output_plot, list.index_site)
index.PLOT.long_LRFF.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "LRFF", occ.LRFF.output_plot, list.index_site)
index.PLOT.long_all.site <- index.PLOT.long.sp(index.TREE.long_plot, bysite, "all", occ.all.output_plot, list.index_site)
index.PLOT.long_site <- index.PLOT.long.merge(index.PLOT.long_BFF.site, index.PLOT.long_GHFF.site, index.PLOT.long_LRFF.site, index.PLOT.long_all.site)
index.PLOT.wide_site <- index.PLOT.wide(index.PLOT.long_site)
count.PLOT.wide_site <- count.PLOT.wide(heights.subset, list.count_site)
count.PLOT.long_site <- count.PLOT.long(count.PLOT.wide_site, bysite)
kernel.PLOT.long_site <- kernel.PLOT.long(kernel.TREE.long_plot, list.pixel_site, occ.BFF.output_plot, occ.GHFF.output_plot, occ.LRFF.output_plot, occ.all.output_plot)
kernel.PLOT.wide_site <- kernel.PLOT.wide(kernel.PLOT.long_site)
kernelnz.PLOT.long_site <- kernel.PLOT.long(kernelnz.TREE.long_plot, list.pixel_site, occ.BFF.output_plot, occ.GHFF.output_plot, occ.LRFF.output_plot, occ.all.output_plot)
kernelnz.PLOT.wide_site <- kernel.PLOT.wide(kernelnz.PLOT.long_site)

##---------------------------------------- Tree as lowest replicate ---------------------------------------- 
bysite <- c("site.code", "session", "site.accession")
index.TREE.Ewide_site <- index.TREE.Ewide(treebat, bysite)
index.TREE.wide_site <- index.TREE.wide(index.TREE.Ewide_site, bysite, list.count_site) 
heights.subset <- calculate.heightdiff(treebat, bytree) %>% 
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% 
  mutate(vert.density = count/(height.diff*value_dirichlet)) %>% 
  mutate(vert.density = ifelse(vert.density==Inf,NA,vert.density))
count.TREE.wide_site <- count.TREE.wide(heights.subset, bysite, list.count_site) 
count.TREE.long_site <- count.TREE.long(count.TREE.wide_site, bysite) 


#########################################
##### Get data into a single frame: #####
#########################################

## Reformat roost data so that species is a column:
roostbat_format <- roostbat %>% 
  dplyr::rename(index.abundance.all = index.abundance) %>% 
  dplyr::rename(pop.estimate.all = pop.estimate.total) %>% 
  melt(id.vars = c("site.code","session","site.accession", "roost.area", "roost.perimeter", "total.trees"), measure.vars = c("pop.estimate.all", "pop.estimate.BFF", "pop.estimate.GHFF", "index.abundance.all", "pop.estimate.LRFF"),
       variable.name = c("species.value"), value.name="value") %>%
  dplyr::mutate(species = str_extract(species.value, "all|BFF|GHFF|LRFF")) %>%
  dplyr::mutate(valuecat = str_extract(species.value, "pop.estimate|index.abundance")) %>%
  dplyr::select(-c(species.value)) %>%
  spread(valuecat,value) 

## Reformat treebat data so that species is a column:
## Also assign whether tree is core or peripherall occupied
bysite <- c("site.code", "session", "site.accession")
occ.all.output_site <- occ.all(treebat, bysite)
occ.BFF.output_site <- occ.BFF(treebat, bysite)
occ.GHFF.output_site <- occ.GHFF(treebat, bysite)
occ.LRFF.output_site <- occ.LRFF(treebat, bysite)
threshold <- 0.8

treebat_all <- treebat %>% 
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #remove cases where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #create total index weight
  ## Assign core trees for all species:
  filter(site.accession %in% occ.all.output_site$site.accession) %>% #choose surveys when at least 1 bat was present, only 
  mutate(tree.occ.all = ifelse(all.index.weight==0, 0, #Assign binary value to show if tree is occupied
                               ifelse(all.index.weight>0, 1,
                                      NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ.all = sum(tree.occ.all)/sum(!is.na(tree.occ.all))) %>%
  mutate(occupancy.cat.all = ifelse(tree.occ.all>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

treebat_BFF <- treebat %>% 
  filter(!is.na(BFF.index.weight)) %>% #remove cases where trees were missed
  ## Assign core trees for BFF species:
  filter(site.accession %in% occ.BFF.output_site$site.accession) %>% #choose surveys when at least 1 BFF was present, only 
  mutate(tree.occ.BFF = ifelse(BFF.index.weight==0, 0, #Assign binary value to show if tree is occupied
                               ifelse(BFF.index.weight>0, 1,
                                      NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ.BFF = sum(tree.occ.BFF)/sum(!is.na(tree.occ.BFF))) %>%
  mutate(occupancy.cat.BFF = ifelse(tree.occ.BFF>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

treebat_GHFF <- treebat %>% 
  filter(!is.na(GHFF.index.weight)) %>% #remove cases where trees were missed
  ## Assign core trees for GHFF species:
  filter(site.accession %in% occ.GHFF.output_site$site.accession) %>% #choose surveys when at least 1 GHFF was present, only 
  mutate(tree.occ.GHFF = ifelse(GHFF.index.weight==0, 0, #Assign binary value to show if tree is occupied
                                ifelse(GHFF.index.weight>0, 1,
                                       NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ.GHFF = sum(tree.occ.GHFF)/sum(!is.na(tree.occ.GHFF))) %>%
  mutate(occupancy.cat.GHFF = ifelse(tree.occ.GHFF>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

treebat_LRFF <- treebat %>% 
  filter(!is.na(LRFF.index.weight)) %>% #remove cases where trees were missed
  ## Assign core trees for LRFF species:
  filter(site.accession %in% occ.LRFF.output_site$site.accession) %>% #choose surveys when at least 1 LRFF was present, only 
  mutate(tree.occ.LRFF = ifelse(LRFF.index.weight==0, 0, #Assign binary value to show if tree is occupied
                                ifelse(LRFF.index.weight>0, 1,
                                       NA))) %>% #nrow of tree.occ should be 2522 trees (this is how many were tagged in total)
  ddply(c("site.code", "tree.accession"), summarise, #for each tree calculate the number of times it was occupied, across surveys where there was one bat present at roost
        tree.occ.LRFF = sum(tree.occ.LRFF)/sum(!is.na(tree.occ.LRFF))) %>%
  mutate(occupancy.cat.LRFF = ifelse(tree.occ.LRFF>=threshold,"Core","Peripheral")) #set threshold for "core" occupancy

treebat_coretrees <- treebat_all %>%
  full_join(treebat_BFF, by = c("site.code", "tree.accession")) %>%
  full_join(treebat_GHFF, by = c("site.code", "tree.accession")) %>%
  full_join(treebat_LRFF, by = c("site.code", "tree.accession")) %>%
  full_join(treebat, by = c("site.code", "tree.accession"))

treebat_format <- treebat_coretrees %>% 
  filter(!is.na(BFF.index.weight)&!is.na(GHFF.index.weight)&!is.na(LRFF.index.weight)) %>% #remove cases where trees were missed
  mutate(all.index.weight = BFF.index.weight+GHFF.index.weight+LRFF.index.weight) %>% #create total index weight
  melt(id.vars = c("tree.accession","session","site.accession", "subplot", "site.code", "crowngroup"), measure.vars = c("BFF.index", "GHFF.index", "LRFF.index", "BFF.index.weight", "GHFF.index.weight", "LRFF.index.weight", "all.index.weight", "tree.occ.all","tree.occ.BFF", "tree.occ.GHFF", "tree.occ.LRFF", "occupancy.cat.all", "occupancy.cat.BFF", "occupancy.cat.GHFF", "occupancy.cat.LRFF"),
       variable.name = c("species.value"), value.name="value") %>%
  dplyr::mutate(species = str_extract(species.value, "all|BFF|GHFF|LRFF")) %>%
  dplyr::mutate(valuecat = str_extract(species.value, "index.weight|index|tree.occ|occupancy.cat")) %>%
  dplyr::select(-c(species.value)) %>%
  spread(valuecat,value) 
#Note that index will have NA values for species=all, because index values aren't sum-able. Values are given for index.weight which are sum-able

## Format value types:
index.TREE.wide_plot <- index.TREE.wide_plot %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(site.accession = as.factor(site.accession)) %>%
  mutate(occ = as.numeric(as.character(occ))) %>%
  mutate(rep = as.factor(rep)) %>%
  mutate(species = as.factor(species)) 
#session and subplot are factors

kernelnz.TREE.wide_plot <- kernelnz.TREE.wide_plot %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(site.accession = as.factor(site.accession)) %>%
  mutate(species = as.factor(species)) 
#session and subplot are factors

roostbat_format <- roostbat_format %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(site.accession = as.factor(site.accession)) %>%
  mutate(species = as.factor(species)) 
#session is a factor

treebat_format <- treebat_format %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(site.accession = as.factor(site.accession)) %>%
  mutate(tree.accession = as.factor(tree.accession)) %>%
  mutate(species = as.factor(species)) %>%
  mutate(occupancy.cat = as.factor(occupancy.cat)) %>%
  mutate(tree.occ = as.numeric(as.character(tree.occ))) %>%
  mutate(crowngroup = as.factor(crowngroup)) %>%
  mutate(index.weight = as.numeric(as.character(index.weight)))
#session and subplot are factors

heights.subset <- heights.subset %>%
  mutate(site.code = as.factor(site.code)) %>%
  mutate(site.accession = as.factor(site.accession)) %>%
  mutate(tree.accession = as.factor(tree.accession)) %>%
  mutate(rep = as.factor(rep)) %>%
  mutate(crowngroup = as.factor(crowngroup)) %>%
  mutate(species = as.factor(species)) 
#session is a factor


## Merge all:
data <- index.TREE.wide_plot %>% #local within-tree abundance (average per plot): "N" (total plot abundance), "tree.count" (number of available trees)
  mutate(prop.occ = occ/tree.count) %>% #create proportion of occupied trees "prop.occ"
  full_join(kernelnz.TREE.wide_plot, by = c("site.code", "session", "site.accession", "subplot", "species", "tree.count")) %>% #plot level kernel density: "mean_v"
  full_join(plotdata, by = c("site.code", "subplot")) %>% #additional plot data: "meandist" is the mean distance between trees per plot
  full_join(roostbat_format, by = c("site.code", "session", "site.accession", "species")) %>% #roost level data: "roost.area", "pop.estimate", "index.abundance"
  full_join(treebat_format, by = c("site.code", "session", "site.accession", "subplot", "species")) %>% #tree-level data: "all.index.weight", "occupancy.cat"
  full_join(heights.subset, by = c("site.code", "session", "site.accession", "subplot", "species", "tree.accession", "crowngroup", "rep")) %>% #tree-level subset data: "count", "height.diff", "vert.density"  
  ## Rename columns for easy reference:
  dplyr::rename(Plot.Abundance = N) %>% 
  dplyr::rename(Plot.Available.Trees = tree.count) %>% 
  dplyr::mutate(Plot.Density.Trees = Plot.Available.Trees/(20*20)) %>% #divide by 20x20 plot area
  dplyr::rename(Plot.Prop.Trees.Occupied = prop.occ) %>% 
  dplyr::rename(Plot.Kernel.Density = mean_v) %>% 
  dplyr::rename(Roost.Area = roost.area) %>% 
  dplyr::rename(Roost.Abundance = pop.estimate) %>%
  dplyr::rename(Roost.Index.Abundance = index.abundance) %>% 
  dplyr::rename(Roost.Available.Trees = total.trees) %>% 
  dplyr::mutate(Roost.Density.Trees = Roost.Available.Trees/(10*20*20)) %>% #divide by 10 x 20x20 plot area
  dplyr::rename(Tree.Abundance.all = index.weight) %>% 
  dplyr::rename(Tree.Abundance.subset = count) %>% 
  dplyr::rename(Tree.Height.Range.subset = height.diff) %>% 
  dplyr::rename(Tree.Density.subset = vert.density) %>% 
  dplyr::rename(Tree.Preference.all = occupancy.cat) %>%
  dplyr::rename(Tree.Occupancy.all = tree.occ) %>% 
  dplyr::rename(NN.distance = meandist) %>% 
  mutate(plotID = paste(site.code, subplot, sep=".")) %>%
  mutate(plotID = as.factor(plotID)) %>%
  mutate(session = as.numeric(as.character(session))) %>%
  mutate(roost.cat = ifelse(Roost.Abundance>=1 & Roost.Abundance<=499,1,
                            ifelse(Roost.Abundance>=500 & Roost.Abundance<=2499,2,
                                   ifelse(Roost.Abundance>=2500 & Roost.Abundance<=4999,3,
                                          ifelse(Roost.Abundance>=5000 & Roost.Abundance<=9999,4,
                                                 ifelse(Roost.Abundance>=10000 & Roost.Abundance<=15999,5,
                                                        ifelse(Roost.Abundance>=16000 & Roost.Abundance<=49999,6,
                                                               ifelse(Roost.Abundance>=50000,7,0
                                                               )))))))) %>%
  mutate(Roost.Index.Abundance = ifelse(species=="BFF",roost.cat,
                                        ifelse(species=="GHFF",roost.cat,
                                               ifelse(species=="LRFF",roost.cat,
                                                      Roost.Index.Abundance 
                                               ))))


## Data is formatted so that species = 'all', 'BFF', "GHFF' and 'LRFF' for all values. Rows are duplicated per species and so need to be filtered to the appropriate species before running models:
data.all <- data %>%
  filter(species == "all") #extract combined species measure only

data.BFF <- data %>%
  filter(species == "BFF") #extract BFF species measure only

data.GHFF <- data %>%
  filter(species == "GHFF") #extract GHFF species measure only

data.LRFF <- data %>%
  filter(species == "LRFF") #extract LRFF species measure only


################################################################################################
##-----------------------------------------Save data------------------------------------------##
################################################################################################
write.csv(roostbat, "Data/Processed/roostbat.csv")
write.csv(kernel.TREE.wide_plot, "Data/Processed/kernel.TREE.wide_plot.csv")
write.csv(kernelnz.TREE.wide_plot, "Data/Processed/kernelnz.TREE.wide_plot.csv")
write.csv(heights.subset, "Data/Processed/heights.subset.csv")
write.csv(treebat, "Data/Processed/treebat.csv")
write.csv(data, "Data/Processed/merged.data.csv")
write.csv(data.all, "Data/Processed/merged.data_all.csv")
write.csv(data.BFF, "Data/Processed/merged.data_BFF.csv")
write.csv(data.GHFF, "Data/Processed/merged.data_GHFF.csv")
write.csv(data.LRFF, "Data/Processed/merged.data_LRFF.csv")

saveRDS(roostbat, "Data/Processed/roostbat.RDS")
saveRDS(kernel.TREE.wide_plot, "Data/Processed/kernel.TREE.wide_plot.RDS")
saveRDS(kernelnz.TREE.wide_plot, "Data/Processed/kernelnz.TREE.wide_plot.RDS")
saveRDS(heights.subset, "Data/Processed/heights.subset.RDS")
saveRDS(treebat, "Data/Processed/treebat.RDS")
saveRDS(data, "Data/Processed/merged.data.RDS")
saveRDS(data.all, "Data/Processed/merged.data_all.RDS")
saveRDS(data.BFF, "Data/Processed/merged.data_BFF.RDS")
saveRDS(data.GHFF, "Data/Processed/merged.data_GHFF.RDS")
saveRDS(data.LRFF, "Data/Processed/merged.data_LRFF.RDS")

## Note that models will not run from these saved .csv files - formatting (e.g. as factors) will be incorrect when read in again
## Use saved .csv files for figure generation only. For running models, run this entire script and use generated files, or load RDS files

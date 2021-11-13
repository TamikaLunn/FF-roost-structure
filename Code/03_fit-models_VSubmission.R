## Title: Code to run GAM model analysis on flying-fox roost structure in SE QLD and NE NSW
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

#roostbat <- readRDS("Data/Processed/roostbat.RDS")
#kernel.TREE.wide_plot <- readRDS("Data/Processed/kernel.TREE.wide_plot.RDS")
#kernelnz.TREE.wide_plot <- readRDS("Data/Processed/kernelnz.TREE.wide_plot.RDS")
#heights.subset <- readRDS("Data/Processed/heights.subset.RDS")
#treebat <- readRDS("Data/Processed/treebat.RDS")
#data <- readRDS("Data/Processed/merged.data.RDS")
data.all <- readRDS("Data/Processed/merged.data_all.RDS")
data.BFF <- readRDS("Data/Processed/merged.data_BFF.RDS")
data.GHFF <- readRDS("Data/Processed/merged.data_GHFF.RDS")
data.LRFF <- readRDS("Data/Processed/merged.data_LRFF.RDS")
## Note that models will not run from  saved .csv files - formatting (e.g. as factors) will be incorrect when read in again
## Use saved .csv files for figure generation only. For running models, run the entire script from '01_download-and-clean-data_VSubmission.R' and use generated files, or load RDS files

##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## Plot.Abundance is estimated from index abundance values, recorded for all trees **per plot**
## Plot.Available.Trees is the number of midstory, canopy and overstory trees **per plot**
## Plot.Density.Trees is the density of midstory, canopy and overstory trees **per plot**
## Plot.Prop.Trees.Occupied is the proportion of trees occupied **per plot**
## Plot.Kernel.Density is the kernel density estimate of bats **per plot**. Estimated with zero kernel values (i.e. blank space) removed
## Roost.Area is the total area (meters squared) of the roost **per roost**, calculated from the roost perimeter
## Roost.Abundance is the abundance estimate of the roost **per roost**, estimated from direct census counts or taken from council estimates
## Roost.Index.Abundance is an estimate of roost abundance **per roost**. Index categories were as follows: 1 = 1-499 bats; 2 = 500-2,499 bats; 3 = 2,500 - 4,999 bats; 4 = 5,000 - 9,999 bats; 5 = 10,000 - 15,999 bats; 6 = 16,000 - 49,999 bats; and 7 = > 50,000 bat
## Roost.Available.Trees is a count of all tagged midstory, canopy and overstory trees within **the roost**
## Roost.Density.Trees is the density of tagged midstory, canopy and overstory trees **per roost**
## Tree.Abundance.all is an estimate of abundance **per tree**, from index abundance values for **all trees**
## Tree.Abundance.subset is a direct count of abundance **per tree**, from a **subset of trees** (N=6 per plot + zero values)
## Tree.Height.Range.subset is the difference between the highest and lowest bat **per tree**, taken for a **subset of trees** only
## Tree.Density.subset is the density of bats **per tree**, estimated as the total count by the height range, for a **subset of trees** only
## Tree.Preference.all indicates tree preference for roosting: whether a tree is occupied =>80% of surveys (core trees=1) or less (peripheral trees=0) **per tree**
## Tree.Occupancy.all indicates tree preference for roosting: is the proportion of times a tree is occupied across the survey **per tree**. This is calculated for surveys when bats are present, only
## NN.distance is the average distance between trees **per plot** (i.e. nearest neighbors)

################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##--------------------------------------------------- Set facet labels for plotting model plots ---------------------------------------------------

##ordered by site, and adjusted for species (i.e. site removed if species absent)
##All/BFF:
site.labs.ordered <- c("Avondale", "Sunnybank", "Burleigh", "Redcliffe", "Toowoomba","Canungra", "Clunes", "Lismore")
names(site.labs.ordered) <- c("01-DAVO", "02-DSUN", "03-DBUR", "04-DRED", "05-DTOW", "06-DCAN", "07-DCLU", "08-DLIS")

##GHFF:
site.labs.ordered_GHFF <- c("Avondale", "Sunnybank", "Redcliffe", "Toowoomba","Canungra", "Clunes", "Lismore")
names(site.labs.ordered_GHFF) <- c("01-DAVO", "02-DSUN", "04-DRED", "05-DTOW", "06-DCAN", "07-DCLU", "08-DLIS")

##LRFF:
site.labs.ordered_LRFF <- c("Avondale", "Sunnybank", "Redcliffe", "Toowoomba")
names(site.labs.ordered_LRFF) <- c("01-DAVO", "02-DSUN", "04-DRED", "05-DTOW")

spp.labs <- c("All species", "Black flying-fox", "Grey-headed flying-fox", "Little red flying-fox")
names(spp.labs) <- c("all", "BFF", "GHFF", "LRFF")

##--------------------------------------------------------------- Fit and save models -------------------------------------------------------------

####-----------------------------------------------------------------------------------------------------------------------------------------------
####--------------------------------------- Tree-level 3-D density (bats/m3) -----------------------------------------------------------------------------------
####-----------------------------------------------------------------------------------------------------------------------------------------------

####################################################### ALL SPECIES ###############################################################################
data.all <- data.all[which(data.all$Tree.Density.subset>0),] ## error with Poisson and Gamma family when response has zero value: 'non-positive values not allowed'
data.all_noNA <- data.all[!is.na(data.all$Tree.Density.subset),]

## View tree data:
jpeg("Output/Model outputs/Tree-level 3-D density/Tree Density_all.jpg", width=1200, height=800)
data.all %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Density.subset)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,5)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Tree-level 3-D density (bats/m3)", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
gam_full_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Density.subset ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_roost.Rds")
saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_M2_roost.Rds")
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_roost.Rds")

## Read in models:
#gam_full_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_roost.Rds")
#gam_M2_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_M2_roost.Rds")
#gam_M3_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_M3_roost.Rds")
#gam_null_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_roost.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()


## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_roost
modlist[[2]] <- gam_M2_roost
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Index Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Index Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"

roost_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Roost-level explanatory variables")
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/roost_model_comparison.Rds")
#readRDS("roost_model_comparison.Rds")
##Best model is gam_full_roost

modelname <- gam_full_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
var3 <- "Roost Tree Density"
roost_model_output <- save.output.gamm_interactionfull_ROOST(modelname, var1, var2, interact, var3, "Best roost-level explanatory model")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level 3-D density/roost_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/roost_model_output.Rds")

####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Density.subset ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/plot_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"
plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level 3-D density/plot_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Density.subset ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/tree_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/tree_model_comparison.Rds")
##Best model is gam_null_tree

modelname <- gam_null_tree
tree_model_output <- save.output.gamm_null(modelname, "Best tree-level explanatory model")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level 3-D density/tree_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/tree_model_output.Rds")

####################################################### BFF ###############################################################################
data.BFF <- data.BFF[which(data.BFF$Tree.Density.subset>0),] ## returning error with Poisson and Gamma family when response has zero value: 'non-positive values not allowed'
data.BFF_noNA <- data.BFF[!is.na(data.BFF$Tree.Density.subset),]

## View tree data:
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Tree Density_all.jpg", width=1200, height=800)
data.BFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Density.subset)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,5)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Tree-level 3-D density (bats/m3)", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
gam_full_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Density.subset ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_roost.Rds")
saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M2_roost.Rds")
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_roost.Rds")

## Read in models:
gam_full_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_roost.Rds")
gam_M2_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_roost
modlist[[2]] <- gam_M2_roost
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model" 

roost_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Roost-level explanatory variables - Black flying-fox")
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/BFF/roost_model_comparison.Rds")
readRDS("roost_model_comparison.Rds")
##Best model is gam_full_roost

modelname <- gam_full_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
var3 <- "Roost Tree Density"
roost_model_output <- save.output.gamm_interactionfull_ROOST(modelname, var1, var2, interact, var3, "Best roost-level explanatory model")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level 3-D density/BFF/roost_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/BFF/roost_model_output.Rds")


####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Density.subset ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M2_plot.Rds") #Did not converge
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Black flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/BFF/plot_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/BFF/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"
plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level 3-D density/BFF/plot_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/BFF/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Density.subset ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/BFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/BFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference"
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables - Black flying-fox")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/BFF/tree_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/BFF/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_full_tree
var1 <- "Tree preference"
tree_model_output <- save.output.gamm_single(modelname, var1, "Best tree-level explanatory model - Black flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level 3-D density/BFF/tree_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/BFF/tree_model_output.Rds")

####################################################### GHFF ###############################################################################
data.GHFF <- data.GHFF[which(data.GHFF$Tree.Density.subset>0),] ## returning error with Poisson and Gamma family when response has zero value: 'non-positive values not allowed'
data.GHFF_noNA <- data.GHFF[!is.na(data.GHFF$Tree.Density.subset),]

## View tree data:
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Tree Density_all.jpg", width=1200, height=800)
data.GHFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Density.subset)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,5)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Tree-level 3-D density (bats/m3)", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
#gam_full_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Density.subset ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
#saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_roost.Rds")
saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M2_roost.Rds")
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_roost.Rds")

## Read in models:
#gam_full_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_roost.Rds")
gam_M2_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
#jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_full_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_full_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()


## Evaluate and save output:
## Not through function as full model didn't converge:
modlist <- list()
modlist[[1]] <- NA
modlist[[2]] <- gam_M2_roost
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"
title <- "Roost-level explanatory variables - Grey-headed flying-fox"

nmodels <- length(modlist)
modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")

for(i in c(2,3,4)){
  modelcomp_output[i,1] <- noquote(modlist_str[[i]])
  modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
  modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
}
modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
for(i in c(1,2,3)){ ## Changed positions because non-converged model is last
  modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
}
modelcomp_output[4,1] <- noquote(modlist_str[[1]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[4,3] <- paste("-")
modelcomp_output[4,2] <- paste("-")
modelcomp_output[4,4] <- paste("-")

## Create table
roost_model_comparison <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
  kable_styling("striped", full_width = F)
#column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
#scroll_box(height = "500px")
#footnote(general = "* Indicates models that did not converge")
#save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/GHFF/roost_model_comparison.Rds")
##Best model is gam_M2_roost

modelname <- gam_M2_roost
var1 <- "Roost Abundance"
var2 <- "Roost Area"
interact <- "Roost Abundance * Roost Area"
var3 <- "Roost Tree Density"

roost_model_output <- save.output.gamm_interaction_ROOST(modelname, var1, var2, interact, "Best roost-level explanatory model - Grey-headed flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level 3-D density/GHFF/roost_model_output.Rds")


####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Density.subset ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Grey-headed flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/GHFF/plot_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/plot_model_comparison.Rds")
##Best model is gam_M2_plot

modelname <- gam_M2_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
plot_model_output <- save.output.gamm_interaction(modelname, var1, var2, interact, "Best subplot-level explanatory model - Grey-headed flying-fox")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level 3-D density/GHFF/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Density.subset ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/GHFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables - Grey-headed flying-fox")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/GHFF/tree_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_full_tree
var1 <- "Tree preference"
tree_model_output <- save.output.gamm_single(modelname, var1, "Best tree-level explanatory model - Grey-headed flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level 3-D density/GHFF/tree_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/GHFF/tree_model_output.Rds")

####################################################### LRFF ###############################################################################
data.LRFF <- data.LRFF[which(data.LRFF$Tree.Density.subset>0),] ## returning error with Poisson and Gamma family when response has zero value: 'non-positive values not allowed'
data.LRFF_noNA <- data.LRFF[!is.na(data.LRFF$Tree.Density.subset),]

## View tree data:
jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Tree Density_all.jpg", width=1200, height=800)
data.LRFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Density.subset)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,5)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Tree-level 3-D density (bats/m3)", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
#gam_full_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
#gam_M2_roost <- gamm(Tree.Density.subset ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Density.subset ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
#saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_roost.Rds") #did not converge
#saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M2_roost.Rds") #did not converge
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M3_roost.Rds") #did not converge
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_roost.Rds")

## Read in models:
#gam_full_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_roost.Rds")
#gam_M2_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_full_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_full_roost$gam)
#dev.off()

#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_full_roost$gam)
#dev.off()

#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()

## Evaluate and save output:
## Not through function as full model and M2 didn't converge:
modlist <- list()
modlist[[1]] <- NA
modlist[[2]] <- NA
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model" 
title <- "Roost-level explanatory variables - Little red flying-fox"

nmodels <- length(modlist)
modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")

for(i in c(3,4)){
  modelcomp_output[i,1] <- noquote(modlist_str[[i]])
  modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
  modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
}
modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
for(i in c(3,4)){ ## Changed positions because non-converged model is last
  modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
}
modelcomp_output[4,1] <- noquote(modlist_str[[2]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[4,3] <- paste("-")
modelcomp_output[4,2] <- paste("-")
modelcomp_output[4,4] <- paste("-")
modelcomp_output[3,1] <- noquote(modlist_str[[1]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[3,3] <- paste("-")
modelcomp_output[3,2] <- paste("-")
modelcomp_output[3,4] <- paste("-")
modelcomp_output[2,1] <- noquote(modlist_str[[2]]) #For the non-converged model, in 3rd position because it's ranked last
#modelcomp_output[2,3] <- paste("-")
#modelcomp_output[2,2] <- paste("-")
#modelcomp_output[2,4] <- paste(round(321.1-320.9),digits=4)

## Create table
roost_model_comparison <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
  kable_styling("striped", full_width = F)
#column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
#scroll_box(height = "500px")
#footnote(general = "* Indicates models that did not converge")
#save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/LRFF/roost_model_comparison.Rds")

## Save output - equal best
modelname <- gam_M3_roost
var3 <- "Roost Tree Density"
roost_model_output <- save.output.gamm_single_ROOST(modelname, var3, "Best roost-level explanatory model - Little red flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level 3-D density/LRFF/roost_model_output.Rds")

## Save output - equal best
modelname <- gam_null_roost
roost_model_output <- save.output.gamm_null_ROOST(modelname, "Best roost-level explanatory model - Little red flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level 3-D density/LRFF/roost_model_comparison-equal best.Rds")


####################################################################
###################### Sublot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Density.subset ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log))  #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Density.subset ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_plot.Rds") #did not converge
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M2_plot.Rds") #did not converge
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Little red flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/LRFF/plot_model_comparison.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"
plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level 3-D density/LRFF/plot_model_output.Rds")
readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/plot_model_output.Rds")

## OR (as AIC only 0.3 different)
modelname <- gam_M2_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
plot_model_output <- save.output.gamm_interaction(modelname, var1, var2, interact, "Best subplot-level explanatory model")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level 3-D density/LRFF/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
#gam_full_tree <- gamm(Tree.Density.subset ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Density.subset ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF_noNA,  family=Gamma(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
#saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_tree.Rds") #did not converge
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_tree.Rds")

## Read models: 
#gam_full_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/LRFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_full_tree.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_full_tree$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
#jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_full_tree$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level 3-D density/LRFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- NA
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

nmodels <- length(modlist)
modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")

for(i in c(2)){
  modelcomp_output[i,1] <- noquote(modlist_str[[i]])
  modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
  modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
}
modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
for(i in c(1)){ ## Changed positions because non-converged model is last
  modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
}
modelcomp_output[2,1] <- noquote(modlist_str[[1]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[2,3] <- paste("-")
modelcomp_output[2,2] <- paste("-")
modelcomp_output[2,4] <- paste("-")

## Create table
tree_model_comparison <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
  kable_styling("striped", full_width = F)
#column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
#scroll_box(height = "500px")
#footnote(general = "* Indicates models that did not converge")
#save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level 3-D density/LRFF/tree_model_comparison.Rds")

modelname <- gam_null_tree
tree_model_output <- save.output.gamm_null(modelname, "Best tree-level explanatory model - Little red flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level 3-D density/LRFF/tree_model_output.Rds")


####-----------------------------------------------------------------------------------------------------------------------------------------------
####--------------------------------------- Tree-level abundance -----------------------------------------------------------------------------------
####-----------------------------------------------------------------------------------------------------------------------------------------------

####################################################### ALL SPECIES ###############################################################################

jpeg("Output/Model outputs/Tree-level abundance/Tree Abundance_all.jpg", width=1200, height=800)
data.all %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Abundance.all)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,25)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Within-tree abundance", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
gam_full_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Abundance.all ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_full_roost.Rds")
saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_M2_roost.Rds")
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_null_roost.Rds")

## Read in models:
gam_full_roost <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_full_roost.Rds")
gam_M2_roost <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_null_roost.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_roost
modlist[[2]] <- gam_M2_roost
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Index Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Index Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"

roost_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Roost-level explanatory variables")
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level abundance/roost_model_comparison.Rds")
##Best model is gam_M2_roost

modelname <- gam_M2_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
roost_model_output <- save.output.gamm_interaction_ROOST(modelname, var1, var2, interact, "Best roost-level explanatory model")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level abundance/roost_model_output.Rds")


####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Abundance.all ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level abundance/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"

plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level abundance/plot_model_output.Rds")



####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Abundance.all ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.all,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level abundance/MODEL_gam_null_tree.Rds")

## Read models:
gam_full_tree <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level abundance/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_full_tree
var1 <- "Tree Preference"
tree_model_output <- save.output.gamm_single(modelname, var1, "Best tree-level explanatory model")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level abundance/tree_model_output.Rds")


####################################################### BFF #################################################################################
jpeg("Output/Model outputs/Tree-level abundance/BFF/Tree Abundance_all.jpg", width=1200, height=800)
data.BFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Abundance.all)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,25)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Within-tree abundance", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
gam_full_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Abundance.all ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_roost.Rds")
saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M2_roost.Rds")
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_roost.Rds")

## Read in models:
gam_full_roost <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_roost.Rds")
gam_M2_roost <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()


## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_roost
modlist[[2]] <- gam_M2_roost
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Index Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Index Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"

roost_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Roost-level explanatory variables - Black flying-fox")
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level abundance/BFF/roost_model_comparison.Rds")
##Best model is gam_full_roost

modelname <- gam_full_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
var3 <- "Roost Tree Density"
roost_model_output <- save.output.gamm_interactionfull_ROOST(modelname, var1, var2, interact, var3, "Best roost-level explanatory model - Black flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level abundance/BFF/roost_model_output.Rds")


####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
#gam_full_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
#gam_M2_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Abundance.all ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
#saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_plot.Rds")
#saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_plot.Rds")

## Read models: 
#gam_full_plot <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_plot.Rds")
#gam_M2_plot <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Black flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level abundance/BFF/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"

plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model - Black flying-fox")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level abundance/BFF/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Abundance.all ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.BFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level abundance/BFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/BFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables - Black flying-fox")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level abundance/BFF/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_null_tree
tree_model_output <- save.output.gamm_null(modelname, "Best tree-level explanatory model - Black flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level abundance/BFF/tree_model_output.Rds")


####################################################### GHFF #################################################################################
jpeg("Output/Model outputs/Tree-level abundance/GHFF/Tree Abundance_all.jpg", width=1200, height=800)
data.GHFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Abundance.all)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,25)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Within-tree abundance", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
gam_full_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
#gam_M2_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Abundance.all ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_roost.Rds")
#saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M2_roost.Rds") ## gam_M2_roost did not converge
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_roost.Rds")

## Read in models:
gam_full_roost <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_roost.Rds")
#gam_M2_roost <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_roost$gam)
dev.off()

#jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()


## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_roost$gam)
dev.off()

#jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()


## Evaluate and save output:
## Not through function as M2 didn't converge:
modlist <- list()
modlist[[1]] <- gam_full_roost
modlist[[2]] <- NA
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"
title <- "Roost-level explanatory variables - Grey-headed flying-fox"

nmodels <- length(modlist)
modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")

for(i in c(1,3,4)){
  modelcomp_output[i,1] <- noquote(modlist_str[[i]])
  modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
  modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
}
modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
for(i in c(1,2,3)){ ## Changed positions because non-converged model is last
  modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
}
modelcomp_output[4,1] <- noquote(modlist_str[[2]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[4,3] <- paste("-")
modelcomp_output[4,2] <- paste("-")
modelcomp_output[4,4] <- paste("-")

## Create table
roost_model_comparison <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
  kable_styling("striped", full_width = F)
#column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
#scroll_box(height = "500px")
#footnote(general = "* Indicates models that did not converge")
#save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level abundance/GHFF/roost_model_comparison.Rds")
##Best model is gam_full_roost

modelname <- gam_full_roost
var1 <- "Roost Abundance"
var2 <- "Roost Area"
interact <- "Roost Abundance * Roost Area"
var3 <- "Roost Tree Density"

roost_model_output <- save.output.gamm_interactionfull_ROOST(modelname, var1, var2, interact, var3, "Best roost-level explanatory model - Grey-headed flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level abundance/GHFF/roost_model_output.Rds")



####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Abundance.all ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Grey-headed flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level abundance/GHFF/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion trees occupied"
var3 <- "Subplot Tree Density"

plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model - Grey-headed flying-fox")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level abundance/GHFF/plot_model_output.Rds")


####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Abundance.all ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.GHFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level abundance/GHFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/GHFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables - Grey-headed flying-fox")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level abundance/GHFF/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_null_tree
tree_model_output <- save.output.gamm_null(modelname, "Best tree-level explanatory model - Grey-headed flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level abundance/GHFF/tree_model_output.Rds")



####################################################### LRFF #################################################################################
jpeg("Output/Model outputs/Tree-level abundance/LRFF/Tree Abundance_all.jpg", width=1200, height=800)
data.LRFF %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  ggplot(aes(session, Tree.Abundance.all)) + 
  geom_point(col="grey") + geom_smooth(col="black") + 
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  coord_cartesian(ylim=c(0,25)) + 
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) +  
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +  
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Within-tree abundance", x="Survey month")
dev.off()

####################################################################
###################### Roost-level comparison ###################### 
####################################################################

## Fit models:
#gam_full_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
#gam_M2_roost <- gamm(Tree.Abundance.all ~ Roost.Index.Abundance*Roost.Area + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_roost <- gamm(Tree.Abundance.all ~ Roost.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_roost <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re"), correlation = corARMA(form= ~ +1|site.code, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline

## Save models:
#saveRDS(gam_full_roost, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_roost.Rds") ## did not converge
#saveRDS(gam_M2_roost, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M2_roost.Rds") ## did not converge
saveRDS(gam_M3_roost, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M3_roost.Rds")
saveRDS(gam_null_roost, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_roost.Rds")

## Read in models:
#gam_full_roost <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_roost.Rds")
#gam_M2_roost <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M2_roost.Rds")
gam_M3_roost <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M3_roost.Rds")
gam_null_roost <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_roost.Rds")

## Check and save model fits:
#jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_full_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_full_roost$gam)
#dev.off()

#jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#gam.check(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_roost$gam)
dev.off()

## Save random effects plots (check seasonal fit)
#jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_full_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_full_roost$gam)
#dev.off()

#jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_M2_roost.jpg", width=1200, height=800)
#par(mfrow = c(2,2))
#plot(gam_M2_roost$gam)
#dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_M3_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_roost$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_null_roost.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_roost$gam)
dev.off()

## Evaluate and save output:
## Not through function as full model and M2 didn't converge:
modlist <- list()
modlist[[1]] <- NA
modlist[[2]] <- NA
modlist[[3]] <- gam_M3_roost
modlist[[4]] <- gam_null_roost

modlist_str <- list()
modlist_str[[1]] <- "Roost Abundance * Roost Area + Roost Density Trees"
modlist_str[[2]] <- "Roost Abundance * Roost Area"
modlist_str[[3]] <- "Roost Density Trees"
modlist_str[[4]] <- "Null Model"
title <- "Roost-level explanatory variables - Little red flying-fox"

nmodels <- length(modlist)
modelcomp_output <- as.data.frame(array(dim=c(nmodels,4)))
colnames(modelcomp_output) <- c("Model structure", "R2", "AIC", "Delta AIC")

for(i in c(3,4)){
  modelcomp_output[i,1] <- noquote(modlist_str[[i]])
  modelcomp_output[i,3] <- round(AIC(summary(modlist[[i]]$lme)), digits=1)
  modelcomp_output[i,2] <- round(summary(modlist[[i]]$gam)[["r.sq"]], digits=3) #coefficient of determination, R2, gives the percentage variation in y explained by x-variables. The range is 0 to 1 (i.e. 0% to 100% of the variation in y can be explained by the x-variables)
}
modelcomp_output <- modelcomp_output[order(modelcomp_output$AIC),]
for(i in c(1,2)){ ## Changed positions because non-converged model is last
  modelcomp_output[i,4] <- modelcomp_output[i,3] - modelcomp_output[1,3]
}
modelcomp_output[3,1] <- noquote(modlist_str[[1]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[3,3] <- paste("-")
modelcomp_output[3,2] <- paste("-")
modelcomp_output[3,4] <- paste("-")
modelcomp_output[4,1] <- noquote(modlist_str[[2]]) #For the non-converged model, in 3rd position because it's ranked last
modelcomp_output[4,3] <- paste("-")
modelcomp_output[4,2] <- paste("-")
modelcomp_output[4,4] <- paste("-")
## Create table
roost_model_comparison <- kable(modelcomp_output[,1:4], caption = title,  escape = F) %>%
  kable_styling("striped", full_width = F)
#column_spec(c(3:11), width = "4cm")  #note - collumn_spec slow for long tables
#scroll_box(height = "500px")
#footnote(general = "* Indicates models that did not converge")
#save_kable(paste(title, ".png", sep="")) #fails to save as .png but saves a nicer looking html than write_file()
saveRDS(roost_model_comparison, file = "Output/Model outputs/Tree-level abundance/LRFF/roost_model_comparison.Rds")
##Best model is gam_null_roost

modelname <- gam_null_roost
roost_model_output <- save.output.gamm_null_ROOST(modelname, "Best roost-level explanatory model - Little red flying-fox")
saveRDS(roost_model_output, file = "Output/Model outputs/Tree-level abundance/LRFF/roost_model_output.Rds")



####################################################################
###################### Subplot-level comparison ###################### 
####################################################################

## Fit models:
gam_full_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M2_plot <- gamm(Tree.Abundance.all ~ Plot.Abundance*Plot.Prop.Trees.Occupied + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_M3_plot <- gamm(Tree.Abundance.all ~ Plot.Density.Trees + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_plot <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_plot, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_plot.Rds")
saveRDS(gam_M2_plot, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M2_plot.Rds")
saveRDS(gam_M3_plot, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M3_plot.Rds")
saveRDS(gam_null_plot, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_plot.Rds")

## Read models: 
gam_full_plot <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_plot.Rds")
gam_M2_plot <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M2_plot.Rds")
gam_M3_plot <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_M3_plot.Rds")
gam_null_plot <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_plot.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_plot$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_full_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_M2_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M2_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_M3_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_M3_plot$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_null_plot.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_plot$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_plot
modlist[[2]] <- gam_M2_plot
modlist[[3]] <- gam_M3_plot
modlist[[4]] <- gam_null_plot

modlist_str <- list()
modlist_str[[1]] <- "Subplot Abundance * Subplot Prop Trees Occupied + Subplot Density Trees" 
modlist_str[[2]] <- "Subplot Abundance * Subplot Prop Trees Occupied" 
modlist_str[[3]] <- "Subplot Density Trees" 
modlist_str[[4]] <- "Null Model" 

plot_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Subplot-level explanatory variables - Little red flying-fox")
saveRDS(plot_model_comparison, file = "Output/Model outputs/Tree-level abundance/LRFF/plot_model_comparison.Rds")
##Best model is gam_full_plot

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"

plot_model_output <- save.output.gamm_interactionfull(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model - Little red flying-fox")
saveRDS(plot_model_output, file = "Output/Model outputs/Tree-level abundance/LRFF/plot_model_output.Rds")



####################################################################
###################### Tree-level comparison ###################### 
####################################################################

## Fit models:
gam_full_tree <- gamm(Tree.Abundance.all ~ Tree.Preference.all + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
gam_null_tree <- gamm(Tree.Abundance.all ~ 1 + s(session, bs="cc", k=2) + s(site.code, bs="re") + s(subplot, bs="re"), correlation = corARMA(form= ~ +1|site.code/subplot, p=1), method = "REML", data=data.LRFF,  family=poisson(link=log)) #k=argument to choose the number of knots. Specify 2 because we have 2 repeat periods. bs='cc' : cyclic cubic regression spline
## note that AR part of model isn't shown in summary output but it is having an effect on the model - when altered the fitted values change

## Save models:
saveRDS(gam_full_tree, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_tree.Rds")
saveRDS(gam_null_tree, file = "Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_tree.Rds")

## Read models: 
gam_full_tree <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_full_tree.Rds")
gam_null_tree <- readRDS("Output/Model outputs/Tree-level abundance/LRFF/MODEL_gam_null_tree.Rds")

## Check and save model fits:
jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
gam.check(gam_null_tree$gam)
dev.off()

## Save random effects plots (check seasonal fit)
jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_full_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_full_tree$gam)
dev.off()

jpeg("Output/Model outputs/Tree-level abundance/LRFF/Random effects_gam_null_tree.jpg", width=1200, height=800)
par(mfrow = c(2,2))
plot(gam_null_tree$gam)
dev.off()

## Evaluate and save output:
modlist <- list()
modlist[[1]] <- gam_full_tree
modlist[[2]] <- gam_null_tree

modlist_str <- list()
modlist_str[[1]] <- "Tree Preference" 
modlist_str[[2]] <- "Null Model" 

tree_model_comparison <- gamm_modelcomp(modlist, modlist_str, "Tree-level explanatory variables - Little red flying-fox")
saveRDS(tree_model_comparison, file = "Output/Model outputs/Tree-level abundance/LRFF/tree_model_comparison.Rds")
##Best model is gam_full_tree

modelname <- gam_full_tree
var1 <- "Tree Preference" 
tree_model_output <- save.output.gamm_single(modelname, var1, "Best tree-level explanatory model - Little red flying-fox")
saveRDS(tree_model_output, file = "Output/Model outputs/Tree-level abundance/LRFF/tree_model_output.Rds")


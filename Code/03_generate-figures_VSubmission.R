## Title: 
## Author: Tamika Lunn, Griffith University
## Version: V1, created 

## Title: Code to visualize roosting structure of flying-fox roosts in SE QLD and NE NSW
## Manuscript: Lunn et al (20201) Counterintuitive scaling between population abundance and local density: implications for modelling transmission of infectious diseases in bat populations
## Author: Tamika Lunn, Griffith University
## Version: Submission, created 12th November 2021

## V1-4 - Manuscript preparation
## VSubmission - Code transferred from V4 "Figures_FF-2_V4-revision.R"

rm(list=ls())

##############################################################################################
##---------------------------------Load data & set functions--------------------------------##
##############################################################################################

### Load helper files and any additional packages
source ("Code/00_functions_VSubmission.R")
library(anchors) #for replace.value
library(binom)

##############################################################################################
##------------------------------------ Overview of data: -----------------------------------##
##############################################################################################

## See data_README

summary(tree.tessellation) #check there are 2,522 and no missing values, max of 199
area <- mean(tree.tessellation$value_dirichlet) # mean crown area across plots
radius <- sqrt(area/3.14159) # mean tree radius across plots. Note this is of a circle


################################################################################################
##-----------------------------------------Start code-----------------------------------------##
################################################################################################

##################
#### Figure 4 ####
##################

##----- Tree-level 3-D density (bats/m3) -----##

#### All species, roost level, best model only ####
gam_M2_roost <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_M2_roost.Rds")

modelname <- gam_M2_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
roost_model_output_df <- save.output.gamm_interaction_ROOST_df(modelname, var1, var2, interact, "Best roost-level explanatory model")

roost_model_output_df ##data in dataframe appropriate for ggplot

#jpeg("Output/Figures/3-D density_roost.jpg", width=1200, height=800)
a <- ggplot(roost_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#dev.off()


#### All species, subplot level, best model only ####
gam_full_plot <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_full_plot.Rds")

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"
plot_model_output_df <- save.output.gamm_interactionfull_df(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")

#jpeg("Output/Figures/3-D density_subplot.jpg", width=1200, height=800)
b <- ggplot(plot_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#facet_grid(.~Effect)
#dev.off()


#### All species, tree level, best model only ####
gam_null_tree <- readRDS("Output/Model outputs/Tree-level 3-D density/MODEL_gam_null_tree.Rds")

modelname <- gam_null_tree
tree_model_output_df <- save.output.gamm_null_df(modelname, "Best tree-level explanatory model")

#jpeg("Output/Figures/3-D density_tree.jpg", width=1200, height=800)
c <- ggplot(tree_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#facet_grid(.~Effect)
#dev.off()


##----- Tree-level abundance -----##

#### All species, roost level, best model only ####

gam_M2_roost <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_M2_roost.Rds")
modelname <- gam_M2_roost
var1 <- "Roost Index Abundance"
var2 <- "Roost Area"
interact <- "Roost Index Abundance * Roost Area"
roost_model_output_df <- save.output.gamm_interaction_ROOST_df(modelname, var1, var2, interact, "Best roost-level explanatory model")

#jpeg("Output/Figures/abundance_roost.jpg", width=1200, height=800)
d <- ggplot(roost_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#dev.off()

#### All species, subplot level, best model only ####
gam_full_plot <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_full_plot.Rds")

modelname <- gam_full_plot
var1 <- "Subplot Abundance"
var2 <- "Proportion Trees Occupied"
interact <- "Subplot Abundance * Proportion Trees Occupied"
var3 <- "Subplot Tree Density"

plot_model_output_df <- save.output.gamm_interactionfull_df(modelname, var1, var2, interact, var3, "Best subplot-level explanatory model")

#jpeg("Output/Figures/abundance_subplot.jpg", width=1200, height=800)
e <- ggplot(plot_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#facet_grid(.~Effect)
#dev.off()


#### All species, tree level, best model only ####
gam_full_tree <- readRDS("Output/Model outputs/Tree-level abundance/MODEL_gam_full_tree.Rds")
modelname <- gam_full_tree
var1 <- "Tree Preference"
tree_model_output_df <- save.output.gamm_single_df(modelname, var1, "Best tree-level explanatory model")

#jpeg("Output/Figures/abundance_tree.jpg", width=1200, height=800)
f <- ggplot(tree_model_output_df, aes(x=Variable, y=Coef, fill=Effect)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=Coef-se, ymax=Coef+se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() + background_grid(major="x", colour.major = "grey95") + 
  theme(axis.text.x = element_text(size=20, angle = 20, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  labs(y="Model coefficient", x="Variable") +
  geom_hline(yintercept=0, linetype="dashed", color = "black")
#facet_grid(.~Effect)
#dev.off()

## Combine panels into single figure
jpeg("Output/Figures/Figure 4.png", width=1400, height=2000)
ggarrange(a+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)), 
          d+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)), 
          b+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)),
          e+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)),
          c+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)), 
          f+coord_cartesian(ylim=c(-6, 11.5))+ggtitle("")+ theme(legend.position="none", axis.text.x = element_text(size=25, angle = 25), axis.text.y = element_text(size=25), axis.title.x = element_text(size=30), axis.title.y = element_text(size=30)), 
          nrow = 3, ncol=2, 
          labels = c("A", "B", "C", "D", "E", "F"), font.label = list(size = 24)) #26_plot occupancy when bats are present at site - total per subplot and proportion of times occupied_ZOOMED_1400 x 800
dev.off()


##################
#### Figure 3 ####
##################

#### Read in data for figure ####
roostbat <- read.csv("Data/Processed/roostbat.csv", row.names = 1)
kernel.TREE.wide_plot <- read.csv("Data/Processed/kernel.TREE.wide_plot.csv", row.names = 1)
kernelnz.TREE.wide_plot <- read.csv("Data/Processed/kernelnz.TREE.wide_plot.csv", row.names = 1)
heights.subset <- read.csv("Data/Processed/heights.subset.csv", row.names = 1)
treebat <- read.csv("Data/Processed/treebat.csv", row.names = 1)
data <- read.csv("Data/Processed/merged.data.csv", row.names = 1)
data.all <- read.csv("Data/Processed/merged.data_all.csv", row.names = 1) %>% mutate(subplot = as.factor(subplot))
data.BFF <- read.csv("Data/Processed/merged.data_BFF.csv", row.names = 1) %>% mutate(subplot = as.factor(subplot))
data.GHFF <- read.csv("Data/Processed/merged.data_GHFF.csv", row.names = 1) %>% mutate(subplot = as.factor(subplot))
data.LRFF <- read.csv("Data/Processed/merged.data_LRFF.csv", row.names = 1) %>% mutate(subplot = as.factor(subplot))

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

## Add on: calculate tree 2D density:
data.all <- data.all %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% ## Merge all data with tesselation data to get crown area for ALL trees not just the height subset
  rename(value_dirichlet_subset = value_dirichlet.x) %>% 
  rename(value_dirichlet_all = value_dirichlet.y) %>%
  mutate(TWOD_Tree.Density.all = Tree.Abundance.all/value_dirichlet_all) ## Calculate new density measure which is just tree abundance by crown area

data.BFF <- data.BFF %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% ## Merge all data with tesselation data to get crown area for ALL trees not just the height subset
  rename(value_dirichlet_subset = value_dirichlet.x) %>% 
  rename(value_dirichlet_all = value_dirichlet.y) %>%
  mutate(TWOD_Tree.Density.all = Tree.Abundance.all/value_dirichlet_all) ## Calculate new density measure which is just tree abundance by crown area

data.GHFF <- data.GHFF %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% ## Merge all data with tesselation data to get crown area for ALL trees not just the height subset
  rename(value_dirichlet_subset = value_dirichlet.x) %>% 
  rename(value_dirichlet_all = value_dirichlet.y) %>%
  mutate(TWOD_Tree.Density.all = Tree.Abundance.all/value_dirichlet_all) ## Calculate new density measure which is just tree abundance by crown area

data.LRFF <- data.LRFF %>%
  dplyr::left_join(
    tree.tessellation[,c("tree.accession", "site.code", "subplot","value_dirichlet")], 
    by = c("tree.accession", "site.code", "subplot")) %>% ## Merge all data with tesselation data to get crown area for ALL trees not just the height subset
  rename(value_dirichlet_subset = value_dirichlet.x) %>% 
  rename(value_dirichlet_all = value_dirichlet.y) %>%
  mutate(TWOD_Tree.Density.all = Tree.Abundance.all/value_dirichlet_all) ## Calculate new density measure which is just tree abundance by crown area


#### Set facet labels for plotting model plots ####

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

est.labs <- c("Tree-level 3-D density", "Tree-level 2-D density","Subplot-level kernel density","Subplot-level density","Roost-level density")
names(est.labs) <- c("01_Within-tree packing", "02_TWOD_tree-density", "03_Kernel-based subplot density", "04_Count-based subplot density", "05_Roost density")


#################################
########## All species ##########
#################################

## Format the different data frames for plotting:
roost_format <- data.all[,c("site.code", "session", "site.accession","Roost.Abundance", "Roost.Index.Abundance", "Roost.Area", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "05_Roost density") %>%
  filter(!is.na(Roost.Abundance)) %>%
  filter(Roost.Abundance>0) %>%
  mutate(Roost.Density = Roost.Abundance/Roost.Area)

plot_format_count <- data.all[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(Count.Plot.Density = as.numeric(as.character(Plot.Abundance))/(20*20)) %>% #bats per meter, as 20 squared = 400
  mutate(est.method = "04_Count-based subplot density") %>%
  filter(Count.Plot.Density > 0) #choose occupied plots only

plot_format_kernel <- data.all[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "03_Kernel-based subplot density") %>%
  filter(Plot.Kernel.Density > 0) #choose occupied plots only

tree_format <- data.all[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "01_Within-tree packing") %>%
  filter(Tree.Density.subset > 0) #choose occupied trees only

tree_format_2D <- data.all[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "TWOD_Tree.Density.all", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "02_TWOD_tree-density") %>%
  filter(TWOD_Tree.Density.all > 0) #choose trees plots only


## Mean Roost-level density 
mean(roost_format$Roost.Density)
Liqr(roost_format$Roost.Density)
Uiqr(roost_format$Roost.Density)

## Mean Subplot-level density  
mean(plot_format_count$Count.Plot.Density)
Liqr(plot_format_count$Count.Plot.Density)
Uiqr(plot_format_count$Count.Plot.Density)

## Mean Subplot-level kernel density  
mean(plot_format_kernel$Plot.Kernel.Density)
Liqr(plot_format_kernel$Plot.Kernel.Density)
Uiqr(plot_format_kernel$Plot.Kernel.Density)

## Mean Tree-level 3-D density 
mean(tree_format$Tree.Density.subset)
Liqr(tree_format$Tree.Density.subset)
Uiqr(tree_format$Tree.Density.subset)

## Mean Tree-level 2-D density 
mean(tree_format_2D$TWOD_Tree.Density.all)
Liqr(tree_format_2D$TWOD_Tree.Density.all)
Uiqr(tree_format_2D$TWOD_Tree.Density.all)

#### Create comparison plot  ####
comp_1_roost <- ggplot() + 
  ## Roost scale density (horizontal):
  #geom_point(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method), size=1) + 
  geom_smooth(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=0.5) +
  
  scale_fill_manual(values=c("#FDE725FF"), labels = c("Roost-level density"))+
  scale_color_manual(values=c("#FDE725FF"), labels = c("Roost-level density"))+
  theme_bw() +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Bat density (number per m^2)", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) # .~ facet into collumns and # ~. facet into rows

comp_1_Ksubplot <- ggplot() + 
  ## Plot scale density from kernel density estimate:
  #geom_point(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=0.5) +
  
  scale_fill_manual(values=c("#39568CFF"), labels = c("Subplot-level kernel density"))+
  scale_color_manual(values=c("#39568CFF"), labels = c("Subplot-level kernel density"))+
  theme_bw() +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Bat density (number per m^2)", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) # .~ facet into collumns and # ~. facet into rows

comp_1_subplot <- ggplot() + 
  ## Plot scale density (horizontal), from plot count/plot area:
  #geom_point(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=0.5) +
  
  scale_fill_manual(values=c("#73D055FF"), labels = c("Subplot-level density"))+
  scale_color_manual(values=c("#73D055FF"), labels = c("Subplot-level density"))+
  theme_bw() +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Bat density (number per m^2)", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) # .~ facet into collumns and # ~. facet into rows

comp_1_tree <- ggplot() + 
  ## Tree scale density (vertical), from tree count/height range:
  #geom_point(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=0.5) +
  
  scale_fill_manual(values=c("#440154FF"), labels = c("Tree-level 3-D density"))+
  scale_color_manual(values=c("#440154FF"), labels = c("Tree-level 3-D density"))+
  theme_bw() +
  background_grid(major="x", colour.major = "grey95")+
  labs(y="Bat density (number per m^3)", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) + 
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(species~Site, labeller = labeller(Site = site.labs.ordered, species = spp.labs)) # .~ facet into collumns and # ~. facet into rows

### Arrange and save plot:
jpeg("Output/Figures/Figure 3_ALL SPECIES.png", width=1200, height=1800)
ggarrange(comp_1_tree+ggtitle("")+theme(legend.position="none")+theme(axis.text.x = element_blank(), axis.title.x = element_blank())+coord_cartesian(ylim=c(0, 10)),
          comp_1_Ksubplot+ggtitle("")+theme(legend.position="none")+theme(axis.text.x = element_blank(), axis.title.x = element_blank())+coord_cartesian(ylim=c(0, 25)),
          comp_1_subplot+ggtitle("")+theme(legend.position="none")+theme(axis.text.x = element_blank(), axis.title.x = element_blank())+coord_cartesian(ylim=c(0, 2)),
          comp_1_roost+ggtitle("")+theme(legend.position="none")+coord_cartesian(ylim=c(0, 3)), 
          nrow = 4, common.legend = TRUE, heights = c(1,1,1,1.5), labels = c("A", "B", "C", "D"), font.label = list(size = 24)) #20_Comparison of methods for density estimation_1500 x 800
dev.off()


comp_1 <- ggplot() + 
  ## Tree scale density (vertical), from tree count/height range:
  #geom_point(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method), size=1) + 
  geom_smooth(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  scale_fill_manual(values=c("#440154FF", "#39568CFF", "#73D055FF", "coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density", "Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  scale_color_manual(values=c("#440154FF","#39568CFF","#73D055FF","coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density","Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  
  theme_bw() +
  background_grid("none")+
  labs(y="Bat density", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=24, angle = 70, hjust = 1),
        axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        plot.title = element_text(size=26), 
        strip.text = element_text(size = 26)) + 
  coord_cartesian(ylim=c(0,20)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(est.method~Site, labeller = labeller(Site = site.labs.ordered, est.method = est.labs), scales="free_y")  # .~ facet into collumns and # ~. facet into rows

jpeg("Output/Figures/Figure 3_ALL SPECIES_constrained axis.png", width=1400, height=1800)
ggarrange(comp_1+ggtitle("")+theme(legend.position="none"))
dev.off()


#################################
############ Black FF ###########
#################################

## Format the different data frames for plotting:
roost_format <- data.BFF[,c("site.code", "session", "site.accession","Roost.Abundance", "Roost.Index.Abundance", "Roost.Area", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "05_Roost density") %>%
  filter(!is.na(Roost.Abundance)) %>%
  filter(Roost.Abundance>0) %>%
  mutate(Roost.Density = Roost.Abundance/Roost.Area)

plot_format_count <- data.BFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(Count.Plot.Density = as.numeric(as.character(Plot.Abundance))/(20*20)) %>% #bats per meter, as 20 squared = 400
  mutate(est.method = "04_Count-based subplot density") %>%
  filter(Count.Plot.Density > 0) #choose occupied plots only

plot_format_kernel <- data.BFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "03_Kernel-based subplot density") %>%
  filter(Plot.Kernel.Density > 0) #choose occupied plots only

tree_format_2D <- data.BFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "TWOD_Tree.Density.all", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "02_TWOD_tree-density") %>%
  filter(TWOD_Tree.Density.all > 0) #choose trees plots only

tree_format <- data.BFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "01_Within-tree packing") %>%
  filter(Tree.Density.subset > 0) #choose occupied plots only


## Create comparison plot:  
comp_1 <- ggplot() + 
  ## Tree scale density (vertical), from tree count/height range:
  #geom_point(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method), size=1) + 
  geom_smooth(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  scale_fill_manual(values=c("#440154FF", "#39568CFF", "#73D055FF", "coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density", "Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  scale_color_manual(values=c("#440154FF","#39568CFF","#73D055FF","coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density","Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  
  theme_bw() +
  background_grid("none")+
  labs(y="Bat density", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=24, angle = 70, hjust = 1),
        axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        plot.title = element_text(size=26), 
        strip.text = element_text(size = 26)) + 
  coord_cartesian(ylim=c(0,20)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(est.method~Site, labeller = labeller(Site = site.labs.ordered, est.method = est.labs), scales="free_y")  # .~ facet into collumns and # ~. facet into rows

jpeg("Output/Figures/Figure 3_BFF_constrained axis.png", width=1400, height=1800)
ggarrange(comp_1+ggtitle("")+theme(legend.position="none"))
dev.off()


#######################################
############ Grey-headed FF ###########
#######################################

## Format the different data frames for plotting:
roost_format <- data.GHFF[,c("site.code", "session", "site.accession","Roost.Abundance", "Roost.Index.Abundance", "Roost.Area", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "05_Roost density") %>%
  filter(!is.na(Roost.Abundance)) %>%
  filter(Roost.Abundance>0) %>%
  mutate(Roost.Density = Roost.Abundance/Roost.Area)

plot_format_count <- data.GHFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(Count.Plot.Density = as.numeric(as.character(Plot.Abundance))/(20*20)) %>% #bats per meter, as 20 squared = 400
  mutate(est.method = "04_Count-based subplot density") %>%
  filter(Count.Plot.Density > 0) #choose occupied plots only

plot_format_kernel <- data.GHFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "03_Kernel-based subplot density") %>%
  filter(Plot.Kernel.Density > 0) #choose occupied plots only

tree_format_2D <- data.GHFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "TWOD_Tree.Density.all", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "02_TWOD_tree-density") %>%
  filter(TWOD_Tree.Density.all > 0) #choose trees plots only

tree_format <- data.GHFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "01_Within-tree packing") %>%
  filter(Tree.Density.subset > 0) #choose occupied plots only


## Create comparison plot:  
comp_1 <- ggplot() + 
  ## Tree scale density (vertical), from tree count/height range:
  #geom_point(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method), size=1) + 
  geom_smooth(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  scale_fill_manual(values=c("#440154FF", "#39568CFF", "#73D055FF", "coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density", "Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  scale_color_manual(values=c("#440154FF","#39568CFF","#73D055FF","coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density","Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  
  theme_bw() +
  background_grid("none")+
  labs(y="Bat density", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=24, angle = 70, hjust = 1),
        axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        plot.title = element_text(size=26), 
        strip.text = element_text(size = 26)) + 
  coord_cartesian(ylim=c(0,20)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(est.method~Site, labeller = labeller(Site = site.labs.ordered, est.method = est.labs), scales="free_y")  # .~ facet into collumns and # ~. facet into rows

jpeg("Output/Figures/Figure 3_GHFF_constrained axis.png", width=1400, height=1800)
ggarrange(comp_1+ggtitle("")+theme(legend.position="none"))
dev.off()


#######################################
############ Little red FF ############
#######################################

## Format the different data frames for plotting:
roost_format <- data.LRFF[,c("site.code", "session", "site.accession","Roost.Abundance", "Roost.Index.Abundance", "Roost.Area", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "05_Roost density") %>%
  filter(!is.na(Roost.Abundance)) %>%
  filter(Roost.Abundance>0) %>%
  mutate(Roost.Density = Roost.Abundance/Roost.Area)

plot_format_count <- data.LRFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(Count.Plot.Density = as.numeric(as.character(Plot.Abundance))/(20*20)) %>% #bats per meter, as 20 squared = 400
  mutate(est.method = "04_Count-based subplot density") %>%
  filter(Count.Plot.Density > 0) #choose occupied plots only

plot_format_kernel <- data.LRFF[,c("site.code", "session", "site.accession", "rep", "Plot.Abundance", "Plot.Kernel.Density", "Ntrees", "Plot.Prop.Trees.Occupied", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "03_Kernel-based subplot density") %>%
  filter(Plot.Kernel.Density > 0) #choose occupied plots only

tree_format_2D <- data.LRFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "TWOD_Tree.Density.all", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  create.Site() %>%
  mutate(est.method = "02_TWOD_tree-density") %>%
  filter(TWOD_Tree.Density.all > 0) #choose trees plots only

tree_format <- data.LRFF[,c("site.code", "session", "site.accession", "rep", "Tree.Abundance.all", "Tree.Abundance.subset", "Tree.Density.subset", "Ntrees", "species")] %>% 
  distinct() %>% ## Remove data that is duplicated per level 
  mutate(est.method = "01_Within-tree packing") %>%
  filter(Tree.Density.subset > 0) #choose occupied plots only 

## Create comparison plot:  
comp_1 <- ggplot() + 
  ## Tree scale density (vertical), from tree count/height range:
  #geom_point(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format, aes(x=as.numeric(as.character(session)), y=Tree.Density.subset, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=tree_format_2D, aes(x=as.numeric(as.character(session)), y=TWOD_Tree.Density.all, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_kernel, aes(x=as.numeric(as.character(session)), y=Plot.Kernel.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method), size=1) + #within-plot averages
  geom_smooth(data=plot_format_count, aes(x=as.numeric(as.character(session)), y=Count.Plot.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  #geom_point(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method), size=1) + 
  geom_smooth(data=roost_format, aes(x=as.numeric(as.character(session)), y=Roost.Density, colour=est.method, fill=est.method), method = "loess", se = TRUE, linetype = "solid", size=1.5) +
  
  scale_fill_manual(values=c("#440154FF", "#39568CFF", "#73D055FF", "coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density", "Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  scale_color_manual(values=c("#440154FF","#39568CFF","#73D055FF","coral","#FDE725FF"), labels = c("Tree-level 3-D density", "Tree-level 2-D density","Subplot-level kernel density","Subplot-level density","Roost-level density"))+
  
  theme_bw() +
  background_grid("none")+
  labs(y="Bat density", x="Survey month", colour="Estimate method", fill="Estimate method")+
  ggtitle("Comparison of methods for density estimation") +
  theme(axis.text.x = element_text(size=24, angle = 70, hjust = 1),
        axis.title.x = element_text(size=26),
        axis.text.y = element_text(size=24),
        axis.title.y = element_text(size=26),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24),
        plot.title = element_text(size=26), 
        strip.text = element_text(size = 26)) + 
  coord_cartesian(ylim=c(0,30)) +
  #scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10, 11, 12, 13), labels=c("Aug 2018","Sept 2018","Oct 2018","Nov 2018","Dec 2018","Jan 2019","Feb 2019","Mar 2019","Apr 2019","May 2019","June 2019", "July 2019", "Aug 2019")) +
  scale_x_continuous(breaks=c(1,4,7,10,13), labels=c("Aug 2018","Nov 2018","Feb 2019","May 2019", "Aug 2019")) +
  facet_grid(est.method~Site, labeller = labeller(Site = site.labs.ordered, est.method = est.labs), scales="free_y")  # .~ facet into collumns and # ~. facet into rows

jpeg("Output/Figures/Figure 3_LRFF_constrained axis.png", width=1400, height=1800)
ggarrange(comp_1+ggtitle("")+theme(legend.position="none"))
dev.off()


##################################################
############# Height vs canopy class #############
##################################################
site.labs <- c("Avondale", "Burleigh", "Canungra", "Clunes", "Lismore", "Redcliffe", "Sunnybank", "Toowoomba")
names(site.labs) <- c("DAVO", "DBUR", "DCAN", "DCLU", "DLIS", "DRED", "DSUN", "DTOW")

crown.labs <- c("Overstory", "Canopy", "Mid-story")
names(crown.labs) <- c("1", "2", "3")

x <- "count"
y <- "height.diff"
ylab <- "Range of vertical tree occupancy (meters)"
xlab <- "Number of bats in tree"
title <- "Number vs height range of bats per tree"

## Facet by site and crown, colour by site
clab <- "site.code"
colour <- "site.code"

p1_site_fac <- ggplot(heights.subset, aes(x=as.numeric(as.character(heights.subset[[x]])))) + 
  theme_bw() +
  #background_grid(major="x", colour.major = "grey95")+
  labs(y=ylab, x=xlab, colour=clab, fill=flab)+
  ggtitle(title) +
  theme(axis.text.x = element_text(size=20, angle = 70, hjust = 1),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        plot.title = element_text(size=22), 
        strip.text = element_text(size = 22)) +
  geom_point(aes(y=heights.subset[[y]]), size=1, colour="black") +
  geom_smooth(aes(y=heights.subset[[y]]), method = "loess", se = TRUE, linetype = "solid", size=0.5, colour="black") +
  #scale_color_manual(values=c("#440154FF", "#404788FF", "#39568CFF", "#238A8DFF","#1F968BFF","#55C667FF","#73D055FF","#FDE725FF"), labels = c("Toowoomba", "Sunnybank","Avondale", "Redcliffe", "Burleigh", "Clunes", "Canungra", "Lismore")) +
  facet_grid(crowngroup~site.code, labeller = labeller(site.code = site.labs, crowngroup=crown.labs), scales="free_x")

### Arrange and save plots:
jpeg("#Output/Figures/Number vs height range of bats per tree_by crown.png", width=1200, height=1000)
#ggarrange(comp_1_tree+ggtitle("")+theme(legend.position="none")+theme(axis.text.x = element_blank(), axis.title.x = element_blank())+coord_cartesian(ylim=c(0, 50)),comp_1_fac+ggtitle("")+theme(legend.position="none")+coord_cartesian(ylim=c(0, 50)), nrow = 2, heights = c(1,1.5), labels = c("A", "B"), font.label = list(size = 24)) #20_Comparison of methods for density estimation_1500 x 800
ggarrange(p1_site_fac+ggtitle("")) #20_Comparison of methods for density estimation_1500 x 800
dev.off()

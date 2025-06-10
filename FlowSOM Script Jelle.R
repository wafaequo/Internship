#Libraries inladen

library(FlowSOM)
library(flowCore)
library(CytoML)
library(flowAI)
library(flowWorkspace)
library(ggcyto)
library(tidyverse)
library(openCyto)
library(knitr)
library(kableExtra)
library(dplyr) 
library(ggplot2)
library(flowWorkspaceData)
library(BiocManager)
library(pheatmap)
library(scattermore)
library(ggpointdensity)
library(readxl)
library(plotly)
library(flowViz)
library(scales)


#increase memory
knitr::opts_chunk$set(cache = TRUE, warning = FALSE, 
                      message = FALSE, cache.lazy = FALSE)

# Read single dataset:
fcs1 <- read.FCS("/Users/jellegerardts/Desktop/Stage/Zuyderland_MNPs/Full panel + MNP BM.fcs",
                 transformation=FALSE, truncate_max_range = FALSE)
colnames(fcs1)

#Compenseren met matrix
  spillover(fcs1)
  fcs1_compensated <- compensate(fcs1, spillover(fcs1)$`$SPILLOVER`) #spillover(fcs1) extracts the spillover matrix from the metadata of the file. 

#Clean
#fcs1_compensated_clean <- flow_auto_qc(fcs1)

#Transform
trans <- estimateLogicle(fcs1_compensated, colnames(fcs1
                                                          [,4:12]))
fcs1_compensated_clean_trans <- transform(fcs1_compensated, trans)
autoplot(fcs1_compensated_clean_trans, x = "PE-A", y = "SSC-A", bins = 200)

# Front scatter - Side Scatter Gating
fcs1_compensated_clean_trans_set <- flowSet(fcs1_compensated_clean_trans)
gs <- GatingSet(fcs1_compensated_clean_trans_set)

#fs_data <- gs_pop_get_data(gs)
#my_gate <- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, 
                                                #channels = c("FSC-A", "SSC-A")))
#gs_pop_add(gs, my_gate, parent = "root", name = "my_gate")
#recompute(gs)

    #Plot the gate
    #p <- autoplot(gs, x = "FSC-A", y = "SSC-A", "my_gate", bins = 200)
    #print(p)

fs_data <- gs_pop_get_data(gs)
Boundary_Gate<- fsApply(fs_data, function(fr) openCyto:::.boundary(fr, 
                                                                   channel = c("FSC-A", "SSC-A"), 
                                                                   min = c(0, 0), max=c(2e5,2e5)))


gs_pop_add(gs, Boundary_Gate, parent = "root", name = "Non-Boundary")
recompute(gs)
autoplot(gs, x = "FSC-A", y = "SSC-A", "Non-Boundary", bins = 200)

#Singlet gating
fs_data <- gs_pop_get_data(gs)
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, 
                                              channels = c("FSC-A", "FSC-H")))
gs_pop_add(gs, singlet_gate, parent = "Non-Boundary", name = "singlets")
recompute(gs)
autoplot(gs, x = "FSC-A", y = "FSC-H", "singlets", bins = 200)

#Gatingset terug naar flowframe veranderen
    subset_data <- gs_pop_get_data(gs, "singlets")
    
  
   # Zet dit om naar een flowFrame
    flow_frame <- as(subset_data, "flowFrame")
    autoplot(flow_frame, x = "V500-C-A", y = "SSC-A",bins = 200)
 
#FlowSOM fuction
names(flow_frame)
fSOM <- FlowSOM(flow_frame, compensate = FALSE, transform = FALSE, 
                toTransform = c(4:16), 
                scale = TRUE,
                colsToUse = c(4, 5, 7, 8, 10, 11, 12), xdim = 6, ydim = 6, nClus = 8,
                seed = 42, rlen = 140)
?FlowSOM
?BuildMST
#Building the selforganizing map
fSOM <- BuildSOM(fSOM, c(4, 5, 7, 8, 10, 11, 12), xdim = 6, ydim = 
                   6, rlen = 140)
str(fSOM$map, max.level = 2)

#Building the minimal spanning tree
fSOM <- BuildMST(fSOM, tSNE= FALSE)
str(fSOM$MST)

#Plot the graph
PlotStars(fSOM, backgroundValues = fSOM$metaclustering, equalNodeSize = TRUE)
?PlotStars
PlotNumbers(fSOM,maxNodeSize = 0.4,equalNodeSize = TRUE)

Plot2DScatters(fSOM, channelpairs = list(c("V500-C-A", "SSC-A")), 
               metaclusters = c(1,2,3,4,5), 
               plotFile = "ps+pmma, NR & Ritdye (rlen = 140, xy =100).png", 
               density = FALSE, 
               centers = FALSE)

#FlowSOM Report
FlowSOMmary(fSOM, plotFile = "FlowSOM pmma + ps (No NR).pdf")


test <- as.data.frame(flow_frame@exprs)
test$`V500-C-A`
fig <- plot_ly(test, x = ~`PE-A`, y = ~`SSC-A`, z = ~`V500-C-A`, type = "scatter3d",
               size = 0.01)

fig

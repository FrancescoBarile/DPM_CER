## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  message = FALSE,
  comment = "#>",
  fig.align = "center",
 fig.width = 7,
fig.height = 5
)

## ----setup-values, include=FALSE----------------------------------------------
n_vec <- c(40, 80, 120, 200)

## ----setup-packages-----------------------------------------------------------
library(tidyverse)
# devtools::install_github("teunbrand/ggh4x")
library(ggh4x) 
library(reporter)

## ----clustering-measures------------------------------------------------------
data(cls_measures, package = "DPmCER")

pp<-ggplot(cls_measures,
           xlab=",",
           aes(x=Model, y=value, fill=Clustering_Measure)) +   theme_bw() +
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank(), legend.position="none",
        legend.title=element_blank() ) +
  geom_violin( )  + 
  facet_grid2(Clustering_Measure~Scenario, 
              labeller = labeller(Clustering_Measure = c("Entropy" = "Entropy",
                                                         "Purity" = "Purity",
                                                         "AdjRandIndex" = "ARI"),
                                  Scenario = c("1"= "low",
                                               "2"= "medium-low",
                                               "3"= "medium",
                                               "4"= "high", 
                                               "5"= "mixed"
                                  ) 
              )
  ) +
  theme( strip.text.x = element_text(size = 15),
         strip.text.y.right  = element_text(size = 15),
         axis.text=element_text(size=15) ) 
pp

## ----distance-measures--------------------------------------------------------
data(distance_measures, package = "DPmCER")

pp <- ggplot(
  distance_measures,
  aes(x = Model, y = value, fill = Distance)
) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank()
  ) +
  geom_violin() +
  facet_grid2(
    Distance ~ Scenario,
    scales = "free_y",
    labeller = labeller(
      Distance = c(
        "KL" = "KL divergence",
        "L1" = paste("L" %p% supsc("1"), "distance")
      ),
      Scenario = c(
        "Scenario 1" = paste0("n = ", n_vec[1]),
        "Scenario 2" = paste0("n = ", n_vec[2]),
        "Scenario 3" = paste0("n = ", n_vec[3]),
        "Scenario 4" = paste0("n = ", n_vec[4])
      )
    )
  ) +
  theme(
    strip.text.x = element_text(size = 15),
    strip.text.y.right = element_text(size = 15),
    axis.text = element_text(size = 15)
  )

pp

## ----core-periphery-----------------------------------------------------------
data(core_periphery, package = "DPmCER")

pp <- ggplot(core_periphery,  xlab=",",
             aes(x=Clustering_Measure, y=value, fill=Clustering_Measure)) +
  geom_violin(alpha = 0.5, bw = 0.015) + 
  scale_fill_manual(name="Clustering_Measure", values=c( rgb(0.9290, 0.6940, 0.1250),
                                                         4,
                                                         rgb(0, 0, 0, 0.5)),
  ) +
  theme_bw() +
  theme(
     axis.title.x=element_blank(),
        axis.text=element_text(size=15),
        axis.title.y=element_blank(), legend.position="none",
        legend.title=element_blank() ) + 
  theme(aspect.ratio = 1) 

pp  


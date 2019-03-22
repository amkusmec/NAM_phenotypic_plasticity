library(tidyverse)
library(stringr)
library(broom)

# Colors for plotting results
colors <- c("royalblue", "orange", "red")

# Vector of phenotype names
phenos <- c("CobDiameter", "CobWeight", "DaystoSilk", "DaysToTassel", "EarDiameter",
            "EarHeight", "EarLength", "EarRankNumber", "EarRowNumber", "EarWeight",
            "GDDAnthesisSilkingInterval", "GDDDaystoSilk", "GDDDaystoTassel",
            "HeightAboveEar", "KernelDepth", "LeafLength", "LeafWidth", 
            "MiddleLeafAngle", "PlantHeight", "TasselLength", 
            "TasselPrimaryBranches", "TotalKernelVolume", "UpperLeafAngle")

# Vector pretty phenotype names for plotting
lbls <- c("Cob Diameter", "Cob Weight", "Days to Silk", "Days to Tassel", "Ear Diameter",
          "Ear Height", "Ear Length", "Ear Rank Number", "Ear Row Number", "Ear Weight",
          "GDD Anthesis Silking Interval", "GDD Days to Silk", "GDD Days to Tassel",
          "Height above Ear", "Kernel Depth", "Leaf Length", "Leaf Width", 
          "Middle Leaf Angle", "Plant Height", "Tassel Length", 
          "Tassel Primary Branches", "Total Kernel Volume", "Upper Leaf Angle")

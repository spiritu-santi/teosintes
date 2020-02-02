  ##################################################################################
  ## Modified script from Fitzpatrick & Keller (2015, Ecology Letters) to model & ## 
  ## map genetic turnover using gradient forests.                                 ##  
  ##                                                                              ##  
  ## This code will continuously be curated and automatized and new               ##
  ## versions and data will be available at https://ramirezbarahona.com/.         ##
  ## Any additional inquery please contact Dr. Aguirre-Liguori at:                ##
  ## jonas_aguirre@hotmail.com                                                    ##
  ##################################################################################  
rm(list=ls())
setwd("/Users/Jonas/Desktop/manuscritos/cambio_clima/PNAS/PNAS_CORRECCIONES/genetic_offset/")
# load libraries that will be used. 
library(raster)
library(gradientForest)


  ##################################################################################
  ## CREATE ACCESORY FUNCTIONS THAT WILL BE NEEDED FOR GRADIENT FOREST            ##
  ##################################################################################
  
  ##################################################################################
  ## Accesory functions to run gradient forest analyses                           ##  
  ## Some functions were obtained from Fitzpatrick and Keller                     ##
  ## (2015,Ecology Letters) available at:                                         ##
  ##                                                                              ##
  ## Code is provided as is, a curated version (annotations) will be available    ##
  ## at www.xxx.com and any additional inquery please contact Dr. Aguirre-Liguori ##
  ## at: jonas_aguirre@hotmail.com                                                ##
  ##################################################################################


  # 1) Function 1: convert_env_trns() constructs the env_trns data.frame that used to predict the 
  # gradient forest model along the landscape across the bioclimatic variables used 
  # this function needs a path to a directory containing the ascii (.asc) with each of the bioclimatic variables
  # brief description: the function reads each layer, extracts the values and appends them to a data frame 
  # who's first columns contains the number of each cell in the layer

  convert_env_trns <- function(path= "./"){ #set path
  file_list <- list.files(path = path,pattern = ".asc")
  ras <- paste(path,file_list[1],sep = "")
  raster_final <- raster(ras)
  temp <- values(raster_final)
  tab_final <- data.frame(cell=1:length(temp))
  nb_layers <- length(file_list)
  for(file in file_list){
    nb_layers <- nb_layers-1
    ras <- raster(paste(path,file,sep = ""))
    ras <- values(ras)
    tab_final <- cbind(tab_final,ras)
    names(tab_final)[ncol(tab_final)] <- file
    print(paste("missing layers:",nb_layers))
    
  }
  
  env_trns <- tab_final
  names_layers <- sub(pattern = ".asc",replacement = "",x = names(env_trns))
  names(env_trns) <-names_layers 
  env_trns <- env_trns[which(!is.na(env_trns$bio_1)),]
  return(env_trns)
  }

  
  
  
  # 2) Function 2: pcaToRaster() is a function ritten by Fitzpatrick and Keller (2015, Ecology Letters)
  # this functions translates the predicted GF model into a RGB raster showing the allelic turnover 
  # across the landscape. For more details see Fitzpatrick and Keller (2015)
  pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
  }


  # 3) Function 3: euclidian_distance()  calculates the euclidian distance between 
  # the predicted turnover in the present and the future data.frames 
  # (it calculates the distance between each bioclimatic layer in the predicted landscapes)
  euclidian_distance <- function(proj_fut=future,pred_pres=present){
  num <- ncol(proj_fut)
  tot <- rep(0,nrow(proj_fut))
  for(i in 1:num ){
    sum <- (proj_fut[,i]-pred_pres[,i])^2
    tot <- tot+sum
  }
  tot <- sqrt(tot)
  return(tot)
  }

  # 4) Function 4: create_raster_dif() creates the raster of euclidian distances
  create_raster_dif <- function(ras=mask,proyeccion=euclidian_distance_model,cells=present_layers$cell){
  mask_temp <- ras
  mask_temp[cells] <- proyeccion
  return(mask_temp)
  }

  ##################################################################################
  ##                         END OF ACCESORY FUNCTIONS                            ##
  ##################################################################################
  
  ##################################################################################
  ## BEGIN SCRIPTS FOR GRADIENT FOREST (GF) MODELING                              ##
  ##################################################################################
  

# part 1) upload datasets and set the objects for GF analyses ------------------------
# load mexicana and parviglumis datasets to run gradient forests. As explained in the main text
# paSNPs, candSNPs were clasified based on Discriminant analyses (Adegenet) performed on candidate
# SNPS described in Aguirre-Liguori et al. (Mol Ecol, 26) ; refSNPs were selected by sampling 500
# SNPs showing FST within the neutral distribution of FST. 
# paSNPs, candSNPs and refSNPs for each subspecies is found in gradient_forest_wild_maize.R file
# loading this file will upload the gfData list that contains a species per slot 
# each slot is a dataframe containing the coordinates, the environmental variables, the allelic frequencies 
# of each SNP and the environmental categories of each population (See main text)

  
  

load("gradient_forest_input.R")

# here is the example for mexicana, for parviglumis modify mexicana
# by parviglumis in all occurrences of the script
  
gfData <- gfData$mexicana
envGF <- gfData[,3:21] # get climate variables

load("/Users/Jonas/Desktop/manuscritos/cambio_clima/PNAS/PNAS_CORRECCIONES/buenos/lista_base_madre.R")
load("/Users/Jonas/Desktop/manuscritos/cambio_clima/PNAS/PNAS_CORRECCIONES/buenos/lista_SNPs_nuevos.R")

nuevos_ref <- setdiff(lista_SNPs$lista_nuevos_ref$mexicana,c(lista_SNPs$SNP_bayenv$mexicana,lista_SNPs$SNP_bayescan$mexicana))
#nuevos_ref <- rownames(lista_base_madre$freq_p)
#nuevos_ref <- sample(rownames(lista_base_madre$freq_p),size = 20000,replace = F)
#nuevos_ref <- lista_SNPs$SNP_bayescan$mexicana




  
frecuencias <- lista_base_madre$freq_p
nuevos_ref <- frecuencias[nuevos_ref,rownames(gfData)]
nuevos_ref <- t(nuevos_ref)

#generate a dataset for each type of SNP
ref_SNPs <- gfData[,grep("refSNP",colnames(gfData))] 
pa_SNPs <- gfData[,grep("paSNP",colnames(gfData))] 
cand_SNPs <- gfData[,grep("candSNP",colnames(gfData))] 

#pa_SNPs <- gfData[,c(grep("paSNP",colnames(gfData)),grep("candSNP",colnames(gfData)))] 

# part 2)  Run gradient forest analyses -------------------------------------------------

maxLevel <- log2(0.368*nrow(envGF)/2) #account for correlations, see ?gradientForest 






print("");print("");print("")

# Fit gf models for each SNP category
gf_paSNPs <- gradientForest(cbind(envGF, pa_SNPs), predictor.vars=colnames(envGF),
                        response.vars=colnames(pa_SNPs), ntree=500, 
                        maxLevel=maxLevel, trace=T, corr.threshold=0.50)
print(gf_paSNPs)
print("");print("");print("")

gf_canSNPs <- gradientForest(cbind(envGF, cand_SNPs), predictor.vars=colnames(envGF),
                               response.vars=colnames(cand_SNPs), ntree=500, 
                               maxLevel=maxLevel, trace=T, corr.threshold=0.50)
print(gf_canSNPs)
print("");print("");print("")

gf_refSNPs <- gradientForest(cbind(envGF, ref_SNPs), predictor.vars=colnames(envGF),
                             response.vars=colnames(ref_SNPs), ntree=500, 
                             maxLevel=maxLevel, trace=T, corr.threshold=0.50)

print(gf_refSNPs)

print("");print("");print("")


gf_nuevos_ref <- gradientForest(cbind(envGF, nuevos_ref), predictor.vars=colnames(envGF),
                                response.vars=colnames(nuevos_ref), ntree=500, 
                                maxLevel=maxLevel, trace=T, corr.threshold=0.50)
print(gf_nuevos_ref)

print("");print("");print("")





# part 3)  plot the GF results  -------------------------------------------------
type = "O"
plot(gf_paSNPs, plot.type=type)
title(sub = "paSNPs")
plot(gf_canSNPs, plot.type=type)
title(sub = "candSNPs")
plot(gf_refSNPs, plot.type=type)
title(sub = "refSNPs")
plot(gf_nuevos_ref, plot.type=type)
title(sub = "nuevos_refSNPs")
#plot R2 distribution for each snp category
boxplot(gf_paSNPs$result,
        gf_canSNPs$result,
        gf_refSNPs$result,
        gf_nuevos_ref$result,
        names = c("paSNPs","candSNPs","refSNPs","nuevos_ref"),
        las=2,
        main="mexicana",
        ylab="R2")

# save the GF functions of each SNP category. This functions describe the contribution of 
# each variable to the function, and allow predicting the turnover across the landscape 
gradient_forests_mexicana <- list(paSNPs=gf_paSNPs,candSNPs=gf_canSNPs,refSNPs=gf_refSNPs,nuevos_ref=gf_nuevos_ref)
save(gradient_forests_mexicana,file = "mexicana/output_GF_mexicana.R")
################################################################################


# part 4) predict allelic turnover -------------------

# first, create a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# function convert_env_trns does this, you need to add a path where the ascii bioclim data is
present_layers <- convert_env_trns(path = "/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/presente/") # set path to the present bioclimatic variables

# second, create a raster layer where the RGB data will be written
raster_final <- raster("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/presente/bio_1.asc") # path to raster layer (i.e. bio_1.asc)
mask <- raster_final

# use the gradient function to project the allelic turnover across the 
# landscape based on the present layers data.frame

bioclimatic <- paste("bio_",1:19,sep = "")
pred_paSNPs <- predict(gf_paSNPs, present_layers[bioclimatic])
pred_candSNPs <- predict(gf_canSNPs, present_layers[bioclimatic])
pred_refSNPs <- predict(gf_refSNPs, present_layers[bioclimatic]) # remove cell column before transforming
pred_nuevos_refSNPs <- predict(gf_nuevos_ref, present_layers[bioclimatic]) # remove cell column before transforming

# plot and save (raster file .tif) allelic turnover for each SNP category
#pa_RGBmap <- pcaToRaster(pred_paSNPs, mask, present_layers$cell)
#plotRGB(pa_RGBmap,main="paSNPs")
#writeRaster(pa_RGBmap, paste(path,"mexicana/mexicana_pa_SNPs_map.tif",sep = ""), format="GTiff", overwrite=TRUE)

#cand_RGBmap <- pcaToRaster(pred_candSNPs, mask, present_layers$cell)
#plotRGB(cand_RGBmap,main="candSNPs")
#writeRaster(cand_RGBmap,paste(path,"mexicana/mexicana_cand_SNPs_map.tif",sep = ""), format="GTiff", overwrite=TRUE)

#ref_RGBmap <- pcaToRaster(pred_refSNPs, mask, present_layers$cell)
#plotRGB(ref_RGBmap,main="refSNPs")
#writeRaster(ref_RGBmap, paste(path,"mexicana_ref_SNPs_map.tif",sep = ""), format="GTiff", overwrite=TRUE)


# part 5) estimate genetic offset  ---------------------------

# for this part, in addition to the present projection, you need to project
# the gradient forest function to a dataframe of future layers (use convert_env_trns)
# containing extracted raster data of a given climatic scenario
# in this case we estimated genetic offset for 8 future climatic variables (see main text)
layer_50_cc_45 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_50/cc45bi50/") #set path to future layer containing 50_cc_45 layers
layer_50_cc_85 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_50/cc85bi50/") #set path to future layer containing 50_cc_85 layers
layer_50_mc_45 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_50/mc45bi50/") #set path to future layer containing 50_mc_45 layers
layer_50_mc_85 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_50/mc85bi50/") #set path to future layer containing 50_mc_85 layers

layer_70_cc_45 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_70/cc45bi70/") #set path to future layer containing 70_cc_45 layers
layer_70_cc_85 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_70/cc85bi70/") #set path to future layer containing 70_cc_85 layers
layer_70_mc_45 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_70/mc45bi70/") #set path to future layer containing 70_mc_45 layers
layer_70_mc_85 <- convert_env_trns("/Users/Jonas/Desktop/manuscritos/cambio_clima/capas_climaticas/futuro/capas_70/mc85bi70/") #set path to future layer containing 70_mc_85 layers

# next predict FUTURE env. variables

proj_paSNPs_50_cc_45 <- predict(gf_paSNPs, layer_50_cc_45[,bioclimatic])
proj_paSNPs_50_cc_85 <- predict(gf_paSNPs, layer_50_cc_85[,bioclimatic])
proj_paSNPs_50_mc_45 <- predict(gf_paSNPs, layer_50_mc_45[,bioclimatic])
proj_paSNPs_50_mc_85 <- predict(gf_paSNPs, layer_50_mc_85[,bioclimatic])

proj_candSNPs_50_cc_45 <- predict(gf_canSNPs, layer_50_cc_45[,bioclimatic])
proj_candSNPs_50_cc_85 <- predict(gf_canSNPs, layer_50_cc_85[,bioclimatic])
proj_candSNPs_50_mc_45 <- predict(gf_canSNPs, layer_50_mc_45[,bioclimatic])
proj_candSNPs_50_mc_85 <- predict(gf_canSNPs, layer_50_mc_85[,bioclimatic])

proj_refSNPs_50_cc_45 <- predict(gf_refSNPs, layer_50_cc_45[,bioclimatic])
proj_refSNPs_50_cc_85 <- predict(gf_refSNPs, layer_50_cc_85[,bioclimatic])
proj_refSNPs_50_mc_45 <- predict(gf_refSNPs, layer_50_mc_45[,bioclimatic])
proj_refSNPs_50_mc_85 <- predict(gf_refSNPs, layer_50_mc_85[,bioclimatic])

proj_nuevo_refSNPs_50_cc_45 <- predict(gf_nuevos_ref, layer_50_cc_45[,bioclimatic])
proj_nuevo_refSNPs_50_cc_85 <- predict(gf_nuevos_ref, layer_50_cc_85[,bioclimatic])
proj_nuevo_refSNPs_50_mc_45 <- predict(gf_nuevos_ref, layer_50_mc_45[,bioclimatic])
proj_nuevo_refSNPs_50_mc_85 <- predict(gf_nuevos_ref, layer_50_mc_85[,bioclimatic])


proj_paSNPs_70_cc_45 <- predict(gf_paSNPs, layer_70_cc_45[,bioclimatic])
proj_paSNPs_70_cc_85 <- predict(gf_paSNPs, layer_70_cc_85[,bioclimatic])
proj_paSNPs_70_mc_45 <- predict(gf_paSNPs, layer_70_mc_45[,bioclimatic])
proj_paSNPs_70_mc_85 <- predict(gf_paSNPs, layer_70_mc_85[,bioclimatic])

proj_candSNPs_70_cc_45 <- predict(gf_canSNPs, layer_70_cc_45[,bioclimatic])
proj_candSNPs_70_cc_85 <- predict(gf_canSNPs, layer_70_cc_85[,bioclimatic])
proj_candSNPs_70_mc_45 <- predict(gf_canSNPs, layer_70_mc_45[,bioclimatic])
proj_candSNPs_70_mc_85 <- predict(gf_canSNPs, layer_70_mc_85[,bioclimatic])

proj_refSNPs_70_cc_45 <- predict(gf_refSNPs, layer_70_cc_45[,bioclimatic])
proj_refSNPs_70_cc_85 <- predict(gf_refSNPs, layer_70_cc_85[,bioclimatic])
proj_refSNPs_70_mc_45 <- predict(gf_refSNPs, layer_70_mc_45[,bioclimatic])
proj_refSNPs_70_mc_85 <- predict(gf_refSNPs, layer_70_mc_85[,bioclimatic])

proj_nuevo_refSNPs_70_cc_45 <- predict(gf_nuevos_ref, layer_70_cc_45[,bioclimatic])
proj_nuevo_refSNPs_70_cc_85 <- predict(gf_nuevos_ref, layer_70_cc_85[,bioclimatic])
proj_nuevo_refSNPs_70_mc_45 <- predict(gf_nuevos_ref, layer_70_mc_45[,bioclimatic])
proj_nuevo_refSNPs_70_mc_85 <- predict(gf_nuevos_ref, layer_70_mc_85[,bioclimatic])


# calculate euclidean distance between current and future genetic spaces  

genOffset_paSNPs_50_cc_45 <- euclidian_distance(proj_fut =proj_paSNPs_50_cc_45,pred_pres = pred_paSNPs)
genOffset_paSNPs_50_cc_85 <- euclidian_distance(proj_fut =proj_paSNPs_50_cc_85,pred_pres = pred_paSNPs)
genOffset_paSNPs_50_mc_45 <- euclidian_distance(proj_fut =proj_paSNPs_50_mc_45,pred_pres = pred_paSNPs)
genOffset_paSNPs_50_mc_85 <- euclidian_distance(proj_fut =proj_paSNPs_50_mc_85,pred_pres = pred_paSNPs)

genOffset_candSNPs_50_cc_45 <- euclidian_distance(proj_fut =proj_candSNPs_50_cc_45,pred_pres = pred_candSNPs)
genOffset_candSNPs_50_cc_85 <- euclidian_distance(proj_fut =proj_candSNPs_50_cc_85,pred_pres = pred_candSNPs)
genOffset_candSNPs_50_mc_45 <- euclidian_distance(proj_fut =proj_candSNPs_50_mc_45,pred_pres = pred_candSNPs)
genOffset_candSNPs_50_mc_85 <- euclidian_distance(proj_fut =proj_candSNPs_50_mc_85,pred_pres = pred_candSNPs)

genOffset_refSNPs_50_cc_45 <- euclidian_distance(proj_fut =proj_refSNPs_50_cc_45,pred_pres = pred_refSNPs)
genOffset_refSNPs_50_cc_85 <- euclidian_distance(proj_fut =proj_refSNPs_50_cc_85,pred_pres = pred_refSNPs)
genOffset_refSNPs_50_mc_45 <- euclidian_distance(proj_fut =proj_refSNPs_50_mc_45,pred_pres = pred_refSNPs)
genOffset_refSNPs_50_mc_85 <- euclidian_distance(proj_fut =proj_refSNPs_50_mc_85,pred_pres = pred_refSNPs)


genOffset_nuevo_refSNPs_50_cc_45 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_50_cc_45,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_50_cc_85 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_50_cc_85,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_50_mc_45 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_50_mc_45,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_50_mc_85 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_50_mc_85,pred_pres = pred_nuevos_refSNPs)


genOffset_paSNPs_70_cc_45 <- euclidian_distance(proj_fut =proj_paSNPs_70_cc_45,pred_pres = pred_paSNPs)
genOffset_paSNPs_70_cc_85 <- euclidian_distance(proj_fut =proj_paSNPs_70_cc_85,pred_pres = pred_paSNPs)
genOffset_paSNPs_70_mc_45 <- euclidian_distance(proj_fut =proj_paSNPs_70_mc_45,pred_pres = pred_paSNPs)
genOffset_paSNPs_70_mc_85 <- euclidian_distance(proj_fut =proj_paSNPs_70_mc_85,pred_pres = pred_paSNPs)

genOffset_candSNPs_70_cc_45 <- euclidian_distance(proj_fut =proj_candSNPs_70_cc_45,pred_pres = pred_candSNPs)
genOffset_candSNPs_70_cc_85 <- euclidian_distance(proj_fut =proj_candSNPs_70_cc_85,pred_pres = pred_candSNPs)
genOffset_candSNPs_70_mc_45 <- euclidian_distance(proj_fut =proj_candSNPs_70_mc_45,pred_pres = pred_candSNPs)
genOffset_candSNPs_70_mc_85 <- euclidian_distance(proj_fut =proj_candSNPs_70_mc_85,pred_pres = pred_candSNPs)

genOffset_refSNPs_70_cc_45 <- euclidian_distance(proj_fut =proj_refSNPs_70_cc_45,pred_pres = pred_refSNPs)
genOffset_refSNPs_70_cc_85 <- euclidian_distance(proj_fut =proj_refSNPs_70_cc_85,pred_pres = pred_refSNPs)
genOffset_refSNPs_70_mc_45 <- euclidian_distance(proj_fut =proj_refSNPs_70_mc_45,pred_pres = pred_refSNPs)
genOffset_refSNPs_70_mc_85 <- euclidian_distance(proj_fut =proj_refSNPs_70_mc_85,pred_pres = pred_refSNPs)

genOffset_nuevo_refSNPs_70_cc_45 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_70_cc_45,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_70_cc_85 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_70_cc_85,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_70_mc_45 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_70_mc_45,pred_pres = pred_nuevos_refSNPs)
genOffset_nuevo_refSNPs_70_mc_85 <- euclidian_distance(proj_fut =proj_nuevo_refSNPs_70_mc_85,pred_pres = pred_nuevos_refSNPs)


## create raster of euclidian distances (genetic offset)
ras_diff_paSNPs_50_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_50_cc_45,cells = present_layers$cell)
ras_diff_paSNPs_50_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_50_cc_85,cells = present_layers$cell)
ras_diff_paSNPs_50_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_50_mc_45,cells = present_layers$cell)
ras_diff_paSNPs_50_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_50_mc_85,cells = present_layers$cell)

ras_diff_canSNPs_50_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_50_cc_45,cells = present_layers$cell)
ras_diff_canSNPs_50_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_50_cc_85,cells = present_layers$cell)
ras_diff_canSNPs_50_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_50_mc_45,cells = present_layers$cell)
ras_diff_canSNPs_50_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_50_mc_85,cells = present_layers$cell)

ras_diff_refSNPs_50_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_50_cc_45,cells = present_layers$cell)
ras_diff_refSNPs_50_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_50_cc_85,cells = present_layers$cell)
ras_diff_refSNPs_50_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_50_mc_45,cells = present_layers$cell)
ras_diff_refSNPs_50_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_50_mc_85,cells = present_layers$cell)


ras_diff_nuevo_refSNPs_50_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_50_cc_45,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_50_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_50_cc_85,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_50_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_50_mc_45,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_50_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_50_mc_85,cells = present_layers$cell)


ras_diff_paSNPs_70_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_70_cc_45,cells = present_layers$cell)
ras_diff_paSNPs_70_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_70_cc_85,cells = present_layers$cell)
ras_diff_paSNPs_70_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_70_mc_45,cells = present_layers$cell)
ras_diff_paSNPs_70_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_paSNPs_70_mc_85,cells = present_layers$cell)

ras_diff_canSNPs_70_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_70_cc_45,cells = present_layers$cell)
ras_diff_canSNPs_70_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_70_cc_85,cells = present_layers$cell)
ras_diff_canSNPs_70_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_70_mc_45,cells = present_layers$cell)
ras_diff_canSNPs_70_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_candSNPs_70_mc_85,cells = present_layers$cell)

ras_diff_refSNPs_70_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_70_cc_45,cells = present_layers$cell)
ras_diff_refSNPs_70_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_70_cc_85,cells = present_layers$cell)
ras_diff_refSNPs_70_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_70_mc_45,cells = present_layers$cell)
ras_diff_refSNPs_70_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_refSNPs_70_mc_85,cells = present_layers$cell)

ras_diff_nuevo_refSNPs_70_cc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_70_cc_45,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_70_cc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_70_cc_85,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_70_mc_45 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_70_mc_45,cells = present_layers$cell)
ras_diff_nuevo_refSNPs_70_mc_85 <-create_raster_dif(ras = mask,proyeccion = genOffset_nuevo_refSNPs_70_mc_85,cells = present_layers$cell)


st_paSNPs <- stack(ras_diff_paSNPs_50_cc_45,
                   ras_diff_paSNPs_50_cc_85,
                   ras_diff_paSNPs_50_mc_45,
                   ras_diff_paSNPs_50_mc_85,
                   ras_diff_paSNPs_70_cc_45,
                   ras_diff_paSNPs_70_cc_85,
                   ras_diff_paSNPs_70_mc_45,
                   ras_diff_paSNPs_70_mc_85)


st_canSNPs <- stack(ras_diff_canSNPs_50_cc_45,
                   ras_diff_canSNPs_50_cc_85,
                   ras_diff_canSNPs_50_mc_45,
                   ras_diff_canSNPs_50_mc_85,
                   ras_diff_canSNPs_70_cc_45,
                   ras_diff_canSNPs_70_cc_85,
                   ras_diff_canSNPs_70_mc_45,
                   ras_diff_canSNPs_70_mc_85)

st_refSNPs <- stack(ras_diff_refSNPs_50_cc_45,
                   ras_diff_refSNPs_50_cc_85,
                   ras_diff_refSNPs_50_mc_45,
                   ras_diff_refSNPs_50_mc_85,
                   ras_diff_refSNPs_70_cc_45,
                   ras_diff_refSNPs_70_cc_85,
                   ras_diff_refSNPs_70_mc_45,
                   ras_diff_refSNPs_70_mc_85)

st_nuevo_refSNPs <- stack(ras_diff_nuevo_refSNPs_50_cc_45,
                    ras_diff_nuevo_refSNPs_50_cc_85,
                    ras_diff_nuevo_refSNPs_50_mc_45,
                    ras_diff_nuevo_refSNPs_50_mc_85,
                    ras_diff_nuevo_refSNPs_70_cc_45,
                    ras_diff_nuevo_refSNPs_70_cc_85,
                    ras_diff_nuevo_refSNPs_70_mc_45,
                    ras_diff_nuevo_refSNPs_70_mc_85)

stack_all_mexicana <- list(st_paSNPs=st_paSNPs,
                           st_canSNPs=st_canSNPs,
                           st_refSNPs=st_refSNPs,
                           st_nuevo_refSNPs=st_nuevo_refSNPs)


save(stack_all_mexicana,file = paste(path,"mexicana/stack_all_mexicana.R",sep = ""))


# part 6) extract population information ------------------------------------------------

# get coordinates
coords <- gfData[,c("X","Y")]

diff_paSNPs_50_cc_45 <- extract(ras_diff_paSNPs_50_cc_45,coords)
diff_paSNPs_50_cc_85 <- extract(ras_diff_paSNPs_50_cc_85,coords)
diff_paSNPs_50_mc_45 <- extract(ras_diff_paSNPs_50_mc_45,coords)
diff_paSNPs_50_mc_85 <- extract(ras_diff_paSNPs_50_mc_85,coords)

diff_canSNPs_50_cc_45 <- extract(ras_diff_canSNPs_50_cc_45,coords)
diff_canSNPs_50_cc_85 <- extract(ras_diff_canSNPs_50_cc_85,coords)
diff_canSNPs_50_mc_45 <- extract(ras_diff_canSNPs_50_mc_45,coords)
diff_canSNPs_50_mc_85 <- extract(ras_diff_canSNPs_50_mc_85,coords)

diff_refSNPs_50_cc_45 <- extract(ras_diff_refSNPs_50_cc_45,coords)
diff_refSNPs_50_cc_85 <- extract(ras_diff_refSNPs_50_cc_85,coords)
diff_refSNPs_50_mc_45 <- extract(ras_diff_refSNPs_50_mc_45,coords)
diff_refSNPs_50_mc_85 <- extract(ras_diff_refSNPs_50_mc_85,coords)

diff_nuevo_refSNPs_50_cc_45 <- extract(ras_diff_nuevo_refSNPs_50_cc_45,coords)
diff_nuevo_refSNPs_50_cc_85 <- extract(ras_diff_nuevo_refSNPs_50_cc_85,coords)
diff_nuevo_refSNPs_50_mc_45 <- extract(ras_diff_nuevo_refSNPs_50_mc_45,coords)
diff_nuevo_refSNPs_50_mc_85 <- extract(ras_diff_nuevo_refSNPs_50_mc_85,coords)

diff_paSNPs_70_cc_45 <- extract(ras_diff_paSNPs_70_cc_45,coords)
diff_paSNPs_70_cc_85 <- extract(ras_diff_paSNPs_70_cc_85,coords)
diff_paSNPs_70_mc_45 <- extract(ras_diff_paSNPs_70_mc_45,coords)
diff_paSNPs_70_mc_85 <- extract(ras_diff_paSNPs_70_mc_85,coords)

diff_canSNPs_70_cc_45 <- extract(ras_diff_canSNPs_70_cc_45,coords)
diff_canSNPs_70_cc_85 <- extract(ras_diff_canSNPs_70_cc_85,coords)
diff_canSNPs_70_mc_45 <- extract(ras_diff_canSNPs_70_mc_45,coords)
diff_canSNPs_70_mc_85 <- extract(ras_diff_canSNPs_70_mc_85,coords)

diff_refSNPs_70_cc_45 <- extract(ras_diff_refSNPs_70_cc_45,coords)
diff_refSNPs_70_cc_85 <- extract(ras_diff_refSNPs_70_cc_85,coords)
diff_refSNPs_70_mc_45 <- extract(ras_diff_refSNPs_70_mc_45,coords)
diff_refSNPs_70_mc_85 <- extract(ras_diff_refSNPs_70_mc_85,coords)

diff_nuevo_refSNPs_70_cc_45 <- extract(ras_diff_nuevo_refSNPs_70_cc_45,coords)
diff_nuevo_refSNPs_70_cc_85 <- extract(ras_diff_nuevo_refSNPs_70_cc_85,coords)
diff_nuevo_refSNPs_70_mc_45 <- extract(ras_diff_nuevo_refSNPs_70_mc_45,coords)
diff_nuevo_refSNPs_70_mc_85 <- extract(ras_diff_nuevo_refSNPs_70_mc_85,coords)


genetic_offset_results_mexicana <- data.frame(coords,
                               diff_paSNPs_50_cc_45=diff_paSNPs_50_cc_45,
                               diff_paSNPs_50_cc_85=diff_paSNPs_50_cc_85,
                               diff_paSNPs_50_mc_45=diff_paSNPs_50_mc_45,
                               diff_paSNPs_50_mc_85=diff_paSNPs_50_mc_85,
                               diff_canSNPs_50_cc_45=diff_canSNPs_50_cc_45,
                               diff_canSNPs_50_cc_85=diff_canSNPs_50_cc_85,
                               diff_canSNPs_50_mc_45=diff_canSNPs_50_mc_45,
                               diff_canSNPs_50_mc_85=diff_canSNPs_50_mc_85,
                               diff_refSNPs_50_cc_45=diff_refSNPs_50_cc_45,
                               diff_refSNPs_50_cc_85=diff_refSNPs_50_cc_85,
                               diff_refSNPs_50_mc_45=diff_refSNPs_50_mc_45,
                               diff_refSNPs_50_mc_85=diff_refSNPs_50_mc_85,
                               diff_nuevo_refSNPs_50_cc_45=diff_nuevo_refSNPs_50_cc_45,
                               diff_nuevo_refSNPs_50_cc_85=diff_nuevo_refSNPs_50_cc_85,
                               diff_nuevo_refSNPs_50_mc_45=diff_nuevo_refSNPs_50_mc_45,
                               diff_nuevo_refSNPs_50_mc_85=diff_nuevo_refSNPs_50_mc_85,
                               
                               diff_paSNPs_70_cc_45=diff_paSNPs_70_cc_45,
                               diff_paSNPs_70_cc_85=diff_paSNPs_70_cc_85,
                               diff_paSNPs_70_mc_45=diff_paSNPs_70_mc_45,
                               diff_paSNPs_70_mc_85=diff_paSNPs_70_mc_85,
                               diff_canSNPs_70_cc_45=diff_canSNPs_70_cc_45,
                               diff_canSNPs_70_cc_85=diff_canSNPs_70_cc_85,
                               diff_canSNPs_70_mc_45=diff_canSNPs_70_mc_45,
                               diff_canSNPs_70_mc_85=diff_canSNPs_70_mc_85,
                               diff_refSNPs_70_cc_45=diff_refSNPs_70_cc_45,
                               diff_refSNPs_70_cc_85=diff_refSNPs_70_cc_85,
                               diff_refSNPs_70_mc_45=diff_refSNPs_70_mc_45,
                               diff_refSNPs_70_mc_85=diff_refSNPs_70_mc_85,
                               diff_nuevo_refSNPs_70_cc_45=diff_nuevo_refSNPs_70_cc_45,
                               diff_nuevo_refSNPs_70_cc_85=diff_nuevo_refSNPs_70_cc_85,
                               diff_nuevo_refSNPs_70_mc_45=diff_nuevo_refSNPs_70_mc_45,
                               diff_nuevo_refSNPs_70_mc_85=diff_nuevo_refSNPs_70_mc_85)


save(genetic_offset_results_mexicana,file = paste(path,"genetic_offset_results_mexicana.R",sep = ""))


# part 7) plot map of genetic offset ---------------------------------------

# set the area that will be plotted.
ext<- list(x=range(coords$X),y=range(coords$Y))
ext[[1]][1] <- ext[[1]][1]-1
ext[[1]][2] <- ext[[1]][2]+1
ext[[2]][1] <- ext[[2]][1]-1
ext[[2]][2] <- ext[[2]][2]+1

ext <- extent(ext)

st_paSNPs <- crop(st_paSNPs,ext)
st_canSNPs <- crop(st_canSNPs,ext)
st_refSNPs <- crop(x = st_refSNPs,ext)
st_nuevo_refSNPs <- crop(x = st_nuevo_refSNPs,ext)

#get coordinates of as many populations for the given subspecies, i.e. mexicana. 
# a large list of localities can be downloaded from CONABIO. 
library(gdata)
maize <- read.table("teosintes_subspecies.csv",sep=",",header=T) # access the CONABIO webpage to download a list of teosintes localities
maize <-  maize[grep("mexicana",maize$Taxa),]
all_paSNPs <- extract(st_paSNPs,maize[,c("Longitud","Latitud")])
all_canSNPs <- extract(st_canSNPs,maize[,c("Longitud","Latitud")])
all_refSNPs <- extract(st_refSNPs,maize[,c("Longitud","Latitud")])
all_nuevo_refSNPs <- extract(st_nuevo_refSNPs,maize[,c("Longitud","Latitud")])


#set layers
layers <- c("50_cc_45","50_cc_85","50_mc_45","50_mc_85","70_cc_45","70_cc_85","70_mc_45","70_mc_85")
# plot the genetic offset of each future projection; the maps are relativez based on the maximum genomic offset value
# found in the subspecies locations. We relativize based on all genomic offset combined to have a relative value across all models
maximum_gen_off <- max(c(all_paSNPs,all_canSNPs,all_refSNPs,all_nuevo_refSNPs),na.rm = T)

for(i in 1:length(layers)){
  st <- st_paSNPs[[i]]
  st[st>=maximum_gen_off]<-NA
  st <- st/maximum_gen_off
  cols <- heat.colors(4)
  cols <- c("white",cols[4:1])
  plot(st,col=colorRampPalette(c(cols))(100),main=paste("paSNPs:",layers[i]))
  puntos <- round(extract(st,coords),digits = 1)*100
  points(coords,bg=colorRampPalette(c(cols))(100)[puntos],pch=21,cex=1.5,lwd=2)
}

for(i in 1:length(layers)){
  st <- st_canSNPs[[i]]
  st[st>=maximum_gen_off]<-NA
  st <- st/maximum_gen_off
  cols <- heat.colors(4)
  cols <- c("white",cols[4:1])
  plot(st,col=colorRampPalette(c(cols))(100),main=paste("canSNPs:",layers[i]))
  puntos <- round(extract(st,coords),digits = 1)*100
  points(coords,bg=colorRampPalette(c(cols))(100)[puntos],pch=21,cex=1.5,lwd=2)
}

for(i in 1:length(layers)){
  st <- st_refSNPs[[i]]
  st[st>=maximum_gen_off]<-NA
  st <- st/maximum_gen_off
  cols <- heat.colors(4)
  cols <- c("white",cols[4:1])
  plot(st,col=colorRampPalette(c(cols))(100),main=paste("refSNPs:",layers[i]))
  puntos <- round(extract(st,coords),digits = 1)*100
  points(coords,bg=colorRampPalette(c(cols))(100)[puntos],pch=21,cex=1.5,lwd=2)
}

for(i in 1:length(layers)){
  st <- st_nuevo_refSNPs[[i]]
  st[st>=maximum_gen_off]<-NA
  st <- st/maximum_gen_off
  cols <- heat.colors(4)
  cols <- c("white",cols[4:1])
  plot(st,col=colorRampPalette(c(cols))(100),main=paste("refSNPs:",layers[i]))
  puntos <- round(extract(st,coords),digits = 1)*100
  points(coords,bg=colorRampPalette(c(cols))(100)[puntos],pch=21,cex=1.5,lwd=2)
}

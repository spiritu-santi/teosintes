  ################################################################################
  ## Code to estimate the distribution of adaptive alleles and estimate the     ## 
  ## potential migration of populations                                         ## 
  ##                                                                            ##    
  ## This code will continuously be curated and automatized and new             ##
  ## versions and data will be available at https://ramirezbarahona.com/.       ##
  ## Any additional inquery please contact Dr. Aguirre-Liguori at:              ##    
  ## jonas_aguirre@hotmail.com                                                  ##
  ################################################################################

# part 1) load libraries
library(raster)
library(gdistance)
library(gdata)


# part 2) run the function that will automize the run for each future projections.
# for details of the function, screen the code within the function. 
function_projection_migration_adaptive_alleles <- function(year=50, # set the year
                                                           rcp=45,  # set the RCP model
                                                           circ.mod="cc", # set circulation model
                                                           subspecies="parviglumis", #set subspecies
                                                           path="./"){ #set path to climatic variables
  subspecies=subspecies
  layer <- paste(year,circ.mod,rcp,sep = "_") # this object will be used to access columns from gradient forest results / name plots
  #modify the next paste function to assure that you have the correct path to the folder that has the
  # future layers
  layer2 <- paste(path,year,"/",circ.mod,rcp,"bi",year,"/",sep = "")
  layer3 <- paste(circ.mod,rcp,"bi",year,sep = "") ## assure that this paste function sets the name of each layer according to world clim
  
  # load genetic_offset_results_subspecies.R file (results from gradient_forest) and input gradient_forest_wild_maize.R
  load("genetic_offset_results_mexicana.R") # this file is created with the "run_genetic_offset" function and has the genetic offset for each climatic variable
  load("gradient_forest_input.R") #input containing the coordinates, variables, allelic frequencies for SNPS of each cateogry and environmental categories of populations
  
  gfData <- gfData[[subspecies]]
  env_cateogories <- gfData[c("X","Y","enviromental_categories")]
  names_SNPs <- sub("paSNPS.","",names(gfData)[grep("paSNP",names(gfData),fixed = T)])
  
  # for the next part of the script it is important to run with maxent the distribution of the allelic alleles See text for description
  # The asc files used for this manuscript can be asked to Dr. Aguirre-Liguori at jonas_aguirre@hotmail.com or found at https://ramirezbarahona.com/
  
  #set path to the maxent distribution of adaptive alleles projected to a particular layer
  path <- paste("./",subspecies,"/",layer,"/",sep = "") #set path
  
  # for loop will generate the raster values for each SNP, will binerize them based on the 
  # log value of the minimum value of mid_cold populations. Finally it will stack them in st()object
  st_final <- stack() # stack all present layers
  st_final_f <- stack() # stack all future layers
  for(j in 1:length(names_SNPs)){
    par(mfrow = c(2,2))
    print(paste("Runing SNP:",names_SNPs[j]))
    print(paste("missing",length(names_SNPs)-j,"SNPs"))
    loci <- names_SNPs[j]
    locus_file <- paste(path,loci,"_avg.asc",sep = "")
    locus_file_f <- paste(path,loci,"_",layer3,"_avg.asc",sep = "")
    ras <- raster(locus_file)
    plot(ras,main=paste("present",loci))
    ras_f <- raster(locus_file_f)
    plot(ras_f,main=paste("future",loci))
    minimum <- extract(ras,env_cateogories[c("X","Y")])
    minimum <- min(minimum[which(env_cateogories$enviromental_categories!="cold")])
    ras[ras<minimum]<-0
    ras[ras>minimum]<-1
    plot(ras,main=paste("present",loci))
    ras_f[ras_f<minimum]<-0
    ras_f[ras_f>minimum]<-1
    plot(ras_f,main=paste("Future",loci))
    st_final <- stack(st_final,ras)
    st_final_f <- stack(st_final_f,ras_f)
  }
  
  ### get and plot all area predicted for at least 5 SNPs
  #sum all stacks
  sum_stack <- st_final[[1]]
  for(i in 2:length(names(st_final))){
    sum_stack <- sum_stack+st_final[[i]]  
  }
  sum_stack[sum_stack<5]<-0
  sum_stack[sum_stack>0]<-1
  
  sum_stack_f <- st_final_f[[1]]
  for(i in 2:length(names(st_final_f))){
    sum_stack_f <- sum_stack_f+st_final_f[[i]]  
  }
  
  #binarize based on the  presence of at least 5 predicted layers. 
  sum_stack_f[sum_stack_f<5]<-0
  sum_stack_f[sum_stack_f>0]<-3
  
  par(mfrow = c(1,2))
  plot(sum_stack,main="present")
  plot(sum_stack_f,main="future")
  
  # create raster showing the current, future and overlaped distirbution
  all_stack <-sum_stack+sum_stack_f
  
  # create directory were all results will be written or plotted
  dir.create(paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,sep = ""))
  pdf(paste("./",subspecies,"/",layer,"/overlap_models_before_cutting.pdf",sep = ""))
  par(mfrow=c(1,1))
  plot(all_stack,main=layer,col=colorRampPalette(c("gray88","deepskyblue","orange","black"))(4))
  dev.off()
  
  
  #Next part: Cut layers based on environmental ranges 
  # run function to obtain the cells that contain information of the raster layer
  convert_env_trns <- function(path= "/Volumes/Seagate 4TB/capas_climaticas/futuro/capas_70/cc45bi70/ascii/"){
    # read a raster layer to define size (nrow) of the data.frame
    file_list <- list.files(path = path,pattern = ".asc")
    ras <- paste(path,file_list[1],sep = "")
    raster_final <- raster(ras)
    temp <- values(raster_final)
    tab_final <- data.frame(cell=1:length(temp))
    nb_layers <- length(file_list)
    # loop extracting each raster layer, getting the 
    # values for each grid and appending to the dataframe
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
  
  
  
  # convert to a dataframe all the variables for the future (path to future layers)
  # and to present (path to present layers)- all layers can be downloaded from Worldclim
  # this step is mainly used to get all the cell numbers of the raster layers
  layers_future <- convert_env_trns(layer2)
  layers_present <-convert_env_trns("/media/OS/Users/Jonas/Desktop/cambio_clima/capas_climaticas/presente/") 
  
  ras_pres <- sum_stack 
  ras_fut <- sum_stack_f
  
  points_pres <- rasterToPoints(ras_pres)
  points_pres <- data.frame(cells=layers_present$cell,X=points_pres[,1],Y=points_pres[,2],Value=points_pres[,3])
  points_future <- rasterToPoints(ras_fut)
  points_future <- data.frame(cells=layers_future$cell,X=points_future[,1],Y=points_future[,2],Value=points_future[,3])
  points_pres <- points_pres[which(points_pres$Value==1),]
  points_future <- points_future[which(points_future$Value==3),]
 
  
  path_pres <- "/media/OS/Users/Jonas/Desktop/cambio_clima/capas_climaticas/presente/"
  path_fut <- layer2
  file.list <- list.files(path_pres,".asc")
  
  # get the environmental variables por each pixel of the present and future models
  for(i in 1:length(file.list)){
    print(paste("missing",length(file.list)-i))
    nam <- paste(path_pres,file.list[i],sep = "")
    ras <- raster(nam)
    ras <- extract(ras,points_pres[,c("X","Y")])
    points_pres$temp <- ras
    names(points_pres)[ncol(points_pres)] <- sub(".asc","",file.list[i])
    
    nam <- paste(path_fut,file.list[i],sep = "")
    ras <- raster(nam)
    ras <- extract(ras,points_future[,c("X","Y")])
    points_future$temp <- ras
    names(points_future)[ncol(points_future)] <- sub(".asc","",file.list[i])
    
  }
  
  points_future$time <- "future"
  points_pres$time <- "present"
  
  # perform a PCA on all the variables of the present and future layers,
  # this way we can compare the models based on the same PC ranges
  pca <- rbind(points_future,points_pres)
  PC <- prcomp(pca[,grep("bio",names(pca),fixed = T)],center = T,scale. = T)
  print(summary(PC))
  PC <- PC$x 
  PC <- PC[,1:6]
  pca <- data.frame(pca,PC)
  points_future <- data.frame(points_future,pca[which(pca$time=="future"),c("PC1","PC2","PC3","PC4","PC5","PC6")])
  points_pres <- data.frame(points_pres,pca[which(pca$time=="present"),c("PC1","PC2","PC3","PC4","PC5","PC6")])
  
  # create a raster layer for each PC value for the present and future layers
  mask <- raster("/media/OS/Users/Jonas/Desktop/cambio_clima/capas_climaticas/presente/bio_1.asc")
  mask[mask>0]<- -20
  mask1 <- mask2 <- mask3 <- mask4 <- mask5 <- mask6  <- mask
  
  mask1[ points_pres$cells] <- points_pres$PC1
  mask2[ points_pres$cells] <- points_pres$PC2
  mask3[ points_pres$cells] <- points_pres$PC3
  mask4[ points_pres$cells] <- points_pres$PC4
  mask5[ points_pres$cells] <- points_pres$PC5
  mask6[ points_pres$cells] <- points_pres$PC6
  
  stack_present <- stack(mask1,mask2,mask3,mask4,mask5,mask6)
  
  mask <- raster(paste(layer2,"bio_1.asc",sep = ""))
  mask[mask>0]<- -20
  mask1 <- mask2 <- mask3 <- mask4 <- mask5 <- mask6  <- mask
  
  mask1[ points_future$cells] <- points_future$PC1
  mask2[ points_future$cells] <- points_future$PC2
  mask3[ points_future$cells] <- points_future$PC3
  mask4[ points_future$cells] <- points_future$PC4
  mask5[ points_future$cells] <- points_future$PC5
  mask6[ points_future$cells] <- points_future$PC6
  
  stack_future <- stack(mask1,mask2,mask3,mask4,mask5,mask6)
  
  # get the PC ranges for warm and mid-warm adapted populations. These ranges will be used
  # to cut the raster layers based on the ranges
  pcs_pobs <- extract(stack_present,env_cateogories[,1:2])
  pcs_pobs <- data.frame(env_cateogories,pcs_pobs)
  if(length(which(pcs_pobs$bio_1.1==-20))>0){
    pcs_pobs <- pcs_pobs[-which(pcs_pobs$bio_1.1==-20),]
  }
  names(pcs_pobs)[4:ncol(pcs_pobs)] <- paste("PC",1:(ncol(pcs_pobs)-3),sep = "")
  warm <- pcs_pobs[which(pcs_pobs$enviromental_categories=="warm" | pcs_pobs$enviromental_categories=="mid_warm"),]
  cold <- pcs_pobs[-which(pcs_pobs$enviromental_categories=="warm" | pcs_pobs$enviromental_categories=="mid_warm"),]
  
  pc1 <- range(warm$PC1)
  pc2 <- range(warm$PC2)
  pc3 <- range(warm$PC3)
  pc4 <- range(warm$PC4)
  pc5 <- range(warm$PC5)
  pc6 <- range(warm$PC6)
  pc_values <- data.frame(PC1=pc1,PC2=pc2,PC3=pc3,PC4=pc4,PC5=pc5,PC6=pc6)
  temp <- grep("PC",names(pcs_pobs),fixed = T)
  for(i in temp[-c(1)]){
    plot(pcs_pobs[,temp[1]],pcs_pobs[,i],xlab=names(pcs_pobs)[temp[1]],ylab=names(pcs_pobs)[i])
    points(warm[,temp[1]],warm[,i],col="red",pch=19)
    points(cold[,temp[1]],cold[,i],col="blue",pch=19)
    abline(v=pc_values[,names(pcs_pobs)[temp[1]]])
    abline(h=pc_values[,names(pcs_pobs)[i]]) 
  }
  
  temp <- rasterToPoints(stack_present)
  temp <- data.frame(cell=layers_present$cell,temp)
  names(temp)[4:ncol(temp)]<-c(paste("PC",1:(ncol(temp)-3),sep = ""))
  temp <-temp[-which(temp[,4]==-20),]
  print(dim(temp))
  temp2 <- temp[which(temp$PC1>=pc1[1] & temp$PC1<=pc1[2] & 
                   temp$PC2>=pc2[1] & temp$PC2<=pc2[2] &
                   temp$PC3>=pc3[1] & temp$PC3<=pc3[2] & 
                   temp$PC4>=pc4[1] & temp$PC4<=pc4[2] ),]
  
  temp2$val <- 1
  
  mask[mask>0]<- -20
  mask[ temp2$cell] <- temp2$val
  mask_pres <- mask
  plot(mask_pres,main="Present")
  
  temp <- rasterToPoints(stack_future)
  temp <- data.frame(cell=layers_future$cell,temp)
  names(temp)[4:ncol(temp)]<-c(paste("PC",1:(ncol(temp)-3),sep = ""))
  temp <-temp[-which(temp[,4]==-20),]
  temp2 <- temp[which(temp$PC1>=pc1[1] & temp$PC1<=pc1[2] & 
                   temp$PC2>=pc2[1] & temp$PC2<=pc2[2] &
                   temp$PC3>=pc3[1] & temp$PC3<=pc3[2] & 
                   temp$PC4>=pc4[1] & temp$PC4<=pc4[2]),] 
  
  temp2$val <- 1
  
  mask[mask>0]<- -20
  mask[ temp2$cell] <- temp2$val
  mask_fut <- mask
  plot(mask_fut,main="Future")
  mask_pres[mask_pres<=0]<-0
  mask_pres[mask_pres>0]<-2
  
  mask_fut[mask_fut<=0]<-0
  mask_fut[mask_fut>0]<-4
  
  temp <- mask_fut+mask_pres
  
  temp2 <- crop(temp,extent(x = list(x=c(-105.85,-95.18),y=c(15.65,22.61))))
  
  # plot the models cut based on the PC ranges 
  pdf(paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,"/overlap_models_cut.pdf",sep = ""))
  plot(temp2,main="Reduced model",col=colorRampPalette(c("gray88","deepskyblue","orange","black"))(4))
  colors <- extract(temp2,env_cateogories[,c("X","Y")])
  colors[colors==0]<-"gray88"
  colors[colors=="2"] <- "deepskyblue"
  colors[colors=="4"] <- "orange"
  colors[colors=="6"] <- "black"
  points(env_cateogories[,c("X","Y")],pch=21,bg=colors)
  dev.off()
  
  areas <- summary(factor(values(temp2)))
  if(is.na(areas["6"])){
    areas_pres <- areas["2"]
    areas_fut <- areas["4"]
    
  }else{
    areas_pres <- areas["2"]+areas["6"]
    areas_fut <- areas["4"]+areas["6"]
  }
  
  nom <- paste(year,"_",circ.mod,"_",rcp,sep="")
  print(nom)
  print(paste("areas_pres:",areas_pres))
  print(paste("areas_fut:",areas_fut))
  areas <- list(pres=areas_pres,fut=areas_fut)
  
  
  # import a dataset containting teosintes populations (for example CONABIO) 
  # and estimate how many popualtions are predicted in the future
  maices <- read.xls("~/Escritorio/clima/sin_duplicados/buenos_fitzpatrick/BaseMaicesNativos.xls")
  maices <- maices[grep(subspecies,maices$Taxa,fixed = T),]
  env_cateogories$projection <- extract(temp,env_cateogories[,c("X","Y")])
  all_projection <- extract(temp,maices[,c("Longitud" ,"Latitud")])
  all_projection <- summary(factor(all_projection))
  pop_projection <- summary(factor(env_cateogories$projection))
  
  if(is.na(all_projection["6"])){
    predicted=as.vector(all_projection["4"])
  }else{
    predicted=as.vector(all_projection["4"]+all_projection["6"]) 
  }
  
  not_predicted <- as.vector(all_projection["0"]+all_projection["2"])
  all_projection <- list(predicted=predicted,not_predicted=not_predicted)
  
  
  
  
  if(is.na(pop_projection["6"])){
    predichos=as.vector(pop_projection["4"])
  }else{
    predichos=as.vector(pop_projection["4"]+pop_projection["6"]) 
  }
  
  no_predichos <- as.vector(pop_projection["0"]+pop_projection["2"])
  pop_projection <- list(predichos=predichos,no_predichos=no_predichos)
  projections <- list(populations=pop_projection,all=all_projection)
  
  superposition <- temp2
  
  # plot the relationship between allelic frquencies the environmental category and projection layer.
  pdf(paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,"/boxplots_categorias.pdf",sep = ""))
  for(i in grep("paSNP",names(gfData))){
    par(mfrow = c(1,2),mar=c(6,5,2,2))
    freq <- gfData[,i]
    plot(freq~factor(env_cateogories$projection),las=2,xlab="",ylab="Frecuency")
    title(main=names(gfData)[i],cex.main=0.6)
    plot(freq~factor(env_cateogories$enviromental_categories),las=2,xlab="",ylab="Frecuency")
    title(main=names(gfData)[i],cex.main=0.6)
  }
  dev.off()
  graphics.off()
  
  all_data <- data.frame(gfData[,-ncol(gfData)],projection=env_cateogories$projection)
  save(all_data,file = paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,"/summary.R",sep = ""))
  
  # get the number of populations (projections)-sampled and from conabio - and the number pixels being predicted at the present and future layerss
  results_projection <- list(projections=projections,areas=areas)
  
  
  ######## final part, the next code estimates migration potential por each population
  # the first part is to set the migration layer, that is the area predicted
  # for both the present and the future. We assume a population 
  # will be able to move actively from present to future
  migration_layer <- superposition
  migration_layer[migration_layer>0]<-1

  migrations <- vector(mode = "list",length = nrow(gfData))
  
  pdf(paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,"/migraciones.pdf",sep = ""))
  
  # the loop analyzes population by population and crops it to search for close regions where the individuals might migrate if thre is a migration layer
  for(i in 1:nrow(gfData)){
    pop <- as.vector(t(gfData[i,c("X","Y")]))
    name <- rownames(gfData[i,])
    ext <- extent(x = list(x=c(pop[1]-0.5,pop[1]+0.5),y=c(pop[2]-0.5,pop[2]+0.5)))
    temp <- crop(superposition,ext)
    plot(temp,main=name,col=colorRampPalette(c("gray88","deepskyblue","orange","black"))(4))
    points(pop[1],pop[2],col="red",pch=4,lwd=4)

    # get all the coordinates predicted in the future: values 4 and 6
    coords <- temp
    coords <- rasterToPoints(coords)
    coords <- coords[which(coords[,3]>=4),]
    # the next two conditionals treat a scenario where the are no migration matrixes and assumes the populations 
    # will no be able to migrate
    if(length(dim(coords))==0){
      mig_capacity <- c(rep(Inf,7))
      names(mig_capacity) <- c("min","10%","20%","30%","40%","50%","max")
      migrations[[i]] <- mig_capacity
      names(migrations)[i] <- name
      next
    }
    if(dim(coords)[1]==0){
      mig_capacity <- c(rep(Inf,7))
      names(mig_capacity) <- c("min","10%","20%","30%","40%","50%","max")
      migrations[[i]] <- mig_capacity
      names(migrations)[i] <- name
      next
    }
    
    # if there is a migration matrix, then the code continues and uses landscape genetics to test
    # potential regions of migration. 
    # first it creates a data frame that contains all the points of potential migration (those only predcted in the future)
    # and the last row containing the coordinate ot the population
    coords <- data.frame(coords)
    coords <- data.frame(x=coords$x,y=coords$y)
    #points(puntos,col="red")
    coords <- rbind(coords,c(x=pop[1],y=pop[2]))
    rownames(coords)<-c(paste("point",1:(nrow(coords)-1),sep = "_"),name)
    
    # from the superposition cropped (temp) we binarize it to 
    # have the potential migration layer
    migration_layer <- temp
    migration_layer[migration_layer>1]<-1
    migration_layer[migration_layer==0]<-NA

    # next we create a transition layer all across the potential migration layer (present and future layers)
    tr_t <- transition(x = migration_layer,transitionFunction = function(x) 1/mean(x),directions = 16)
    #we calculate all the distances between potential points and the population
    distances <- costDistance(x = tr_t,fromCoords = as.matrix(coords))
    distances <- as.matrix(distances)
    distances <- distances[nrow(distances),]
    distances <- distances[-length(distances)]
    distances <- distances
    # we estimate the distribution of potential migrations 
    probs_migration <- quantile(x = distances,probs = c(0.1,0.2,0.3,0.4,0.5))
    colors <- heat.colors(6)
    if(min(distances)==Inf){
      finals <- coords[which(distances==min(distances)),]
      points(finals,cex=0.5,col="black",pch=19)
      distances <- c(min(distances),probs_migration,max(distances))
      names(distances)[c(1,length(distances))] <- c("Min","Max")
      migrations[[i]] <- distances
      names(migrations)[i] <- name
      next
    }
    # we plot based on colors the possible areas of migration, based on the distance estimated with the Gdistance package
    finals <- coords[which(distances==min(distances)),]
    points(finals,cex=1,pch=21,bg=colors[1])
    finals <- coords[which(distances==probs_migration[1]),]
    points(finals,cex=1,pch=21,bg=colors[2])
    finals <- coords[which(distances==probs_migration[2]),]
    points(finals,cex=1,pch=21,bg=colors[3])
    finals <- coords[which(distances==probs_migration[3]),]
    points(finals,cex=1,pch=21,bg=colors[4])
    finals <- coords[which(distances==probs_migration[4]),]
    points(finals,cex=1,pch=21,bg=colors[5])
    finals <- coords[which(distances==probs_migration[5]),]
    points(finals,cex=1,pch=21,bg=colors[6])
    
    
    distances <- c(min(distances),probs_migration,max(distances))
    names(distances)[c(1,length(distances))] <- c("Min","Max")
    migrations[[i]] <- distances
    names(migrations)[i] <- name
  }
  
  # we plot the possible migration distance for each population
  mig <- as.data.frame(migrations)
  names(mig)[which(mig[1,]==Inf)]
  if(length(which(mig[1,]!=Inf))!=0){
    barplot(as.matrix(mig[1,which(mig[1,]!=Inf)]),las=2,main="Min",ylim=c(0,40))
  }else{plot(1,type="n",main="Min");text(1,"No migrations")}
  if(length(which(mig[2,]!=Inf))!=0){
    barplot(as.matrix(mig[2,which(mig[2,]!=Inf)]),las=2,main="10%",ylim=c(0,40))
  }else{plot(1,type="n",main="10%");text(1,"No migrations")}
  if(length(which(mig[3,]!=Inf))!=0){
    barplot(as.matrix(mig[3,which(mig[3,]!=Inf)]),las=2,main="20%",ylim=c(0,40))
  }else{plot(1,type="n",main="20%");text(1,"No migrations")}
  if(length(which(mig[4,]!=Inf))!=0){
    barplot(as.matrix(mig[4,which(mig[4,]!=Inf)]),las=2,main="30%",ylim=c(0,40))
  }else{plot(1,type="n",main="30%");text(1,"No migrations")}
  if(length(which(mig[5,]!=Inf))!=0){
    barplot(as.matrix(mig[5,which(mig[5,]!=Inf)]),las=2,main="40%",ylim=c(0,40))
  }else{plot(1,type="n",main="40%");text(1,"No migrations")}
  if(length(which(mig[6,]!=Inf))!=0){
    barplot(as.matrix(mig[6,which(mig[6,]!=Inf)]),las=2,main="50%",ylim=c(0,40))
  }else{plot(1,type="n",main="50%");text(1,"No migrations")}
  
  dev.off()
  # we save the distribution of potential migration for each population for a given environemtal layer
  save(migrations,file = paste("~/Escritorio/clima/functions_migration/",subspecies,"/",layer,"/migrations.R",sep = ""))
  return(results_projection)
}



### after running the function, now we create, for each future layer, a list containing the variables for each model



cc_50_45 <- list(year=50,rcp=45,circ.mod="cc",subspecies="mexicana")
cc_50_85 <- list(year=50,rcp=85,circ.mod="cc",subspecies="mexicana")
mc_50_45 <- list(year=50,rcp=45,circ.mod="mc",subspecies="mexicana")
mc_50_85 <- list(year=50,rcp=85,circ.mod="mc",subspecies="mexicana")
cc_70_45 <- list(year=70,rcp=45,circ.mod="cc",subspecies="mexicana")
cc_70_85 <- list(year=70,rcp=85,circ.mod="cc",subspecies="mexicana")
mc_70_45 <- list(year=70,rcp=45,circ.mod="mc",subspecies="mexicana")
mc_70_85 <- list(year=70,rcp=85,circ.mod="mc",subspecies="mexicana")

# for each layer, we run the function
cc_50_45 <- function_projection_migration_adaptive_alleles(year = cc_50_45$year,rcp = cc_50_45$rcp,circ.mod = cc_50_45$circ.mod,subspecies=cc_50_45$subspecies,path=path)
cc_50_85 <- function_projection_migration_adaptive_alleles(year = cc_50_85$year,rcp = cc_50_85$rcp,circ.mod = cc_50_85$circ.mod,subspecies=cc_50_85$subspecies,path=path)
mc_50_45 <- function_projection_migration_adaptive_alleles(year = mc_50_45$year,rcp = mc_50_45$rcp,circ.mod = mc_50_45$circ.mod,subspecies=mc_50_45$subspecies,path=path)
mc_50_85 <- function_projection_migration_adaptive_alleles(year = mc_50_85$year,rcp = mc_50_85$rcp,circ.mod = mc_50_85$circ.mod,subspecies=mc_50_85$subspecies,path=path)
cc_70_45 <- function_projection_migration_adaptive_alleles(year = cc_70_45$year,rcp = cc_70_45$rcp,circ.mod = cc_70_45$circ.mod,subspecies=cc_70_45$subspecies,path=path)
cc_70_85 <- function_projection_migration_adaptive_alleles(year = cc_70_85$year,rcp = cc_70_85$rcp,circ.mod = cc_70_85$circ.mod,subspecies=cc_70_85$subspecies,path=path)
mc_70_45 <- function_projection_migration_adaptive_alleles(year = mc_70_45$year,rcp = mc_70_45$rcp,circ.mod = mc_70_45$circ.mod,subspecies=mc_70_45$subspecies,path=path)
mc_70_85 <- function_projection_migration_adaptive_alleles(year = mc_70_85$year,rcp = mc_70_85$rcp,circ.mod = mc_70_85$circ.mod,subspecies=mc_70_85$subspecies,path=path)

#we save the results
mexicana <- list(cc_50_45=cc_50_45,cc_50_85=cc_50_85,cc_70_45=cc_70_45,cc_70_85=cc_70_85,mc_50_45=mc_50_45,mc_50_85=mc_50_85,mc_70_45=mc_70_45,mc_70_85=mc_70_85)
save(mexicana,file = "~/Escritorio/clima/sin_duplicados/buenos_fitzpatrick/mexicana/areas_mexicana.R")



####################################
### WB MODEL OUTPUTS PROCESSING ####
####################################

WBModel.output.process <- function ( runID             = Jaar_model,
                                               FileWildBoarMatO =    "./Input/MatrixLandscape.csv", #
                                               outputfolder     = "./Output/Example_run/",  #
                                               maxit            = aantal_iterations,
                                               MaxDays          = aantal_jaar*365,
                                               mytif_name       =  "./Input/Habitat_suitability.tif",
                                               mode             = "noncumulative",
                                               outputindex      = ""
                                               ) {
  #-----------------------------------
  
  
  library(raster)
  library(RColorBrewer)
  library(dplyr)
  
  NAMETG    <- paste0(outputfolder, runID,"-DayOutPopMat.txt")
  NAMETA    <- paste0(outputfolder, runID,"-DayOutAniMat.txt")
  NAMETP    <- paste0(outputfolder, runID,"-DayOutPixels.txt")
  NAMETP7    <- paste0(outputfolder, runID,"-DayOutPixels7.txt")
  NAMEYP     <- paste0(outputfolder, runID,"-YearPixels.txt")
  DayOutPopMat2    <- read.table(NAMETG,sep=" ")
  DayOutAniMat2    <- read.table(NAMETA,sep=" ")
  YearPix  <- read.table(NAMEYP, sep=" ")
  
  
  
  ## POPULATION EVOLUTION IN TIME - ALL
  # # Code can be uncommented to obtain the graph of the population evolution over time
  # mypdfpop <- paste0(outputfolder, "PopGraph-", runID, "-", outputindex, ".pdf")
  # pdf(mypdfpop)
  # myiter <- 1  # iteration
  # plot(1:MaxDays, DayOutAniMat2[,myiter],
  #      xlab='Time (days)', ylab='Wild boar population (# animals)',
  #      main= "Wild boar population over time, by iteration",
  #      lwd=2, type='l',
  #      ylim = c(0, max(DayOutAniMat2)))
  # abline(v = (0:(MaxDays/365))*365, col='blue', lwd=2)
  # for (i in 2:maxit){
  #   myiter <- i
  #   lines(1:MaxDays, DayOutAniMat2[, myiter],
  #         lwd = 2,
  #         type='l')
  # }
  # 
  # DayOutAniMat2.no0 <- DayOutAniMat2
  # DayOutAniMat2.no0[DayOutAniMat2 == 0] <- NA
  # DayOutAniMat2.no0$Mean <- apply(DayOutAniMat2.no0, 1, mean, na.rm=T)
  # 
  # lines(1:dim(DayOutAniMat2.no0)[1], DayOutAniMat2.no0$Mean, pch=20, type='l', col="green",lwd=3 )
  # 
  # legend('topleft', inset= 0.02,
  #        legend=c("Each iteration", "Mean value", "1st January"),
  #        col=c("black", "green", "blue"), lty=1, lwd= 2, cex=0.6, bg="white")
  # 
  # dev.off()
  
  
  ####################### 
  
  ## MAP 
  
  WBMatO <- read.table(FileWildBoarMatO, sep=";", header = T, dec=',')
  #names(WBMatO) <- c('ID', 'XCoord', 'YCoord', 'value', 'E', 'F', 'G', 'H', 'I','J', 'K', 'L') 
  mymap <- raster(mytif_name)
  mymap <- mymap*1
  
  
  ######## HEATMAP ############
  original.mapValues <- mymap@data@values 
  dayskip   <- apply(YearPix, 2, FUN = function(x) (which(x <1)[1:(1+MaxDays/365)])) 
  
  #hmp# mypdffileHM <- paste0(outputfolder, "Heatmap_", mode, "-", runID,"-YP-CATEG-10_40_70-bis-", outputindex, ".pdf")  
  #hmp# pdf(mypdffileHM)
  #hmp# plot(mymap, col=c('grey', 'yellow', 'green', 'dark green'), main= "Day 0")  
  # #these lines #hmp# can be uncommented to obtain the pdf with heatmaps
  
  mytimeafterinfHM <- 1:(1+(MaxDays/365))
  
  #### SUMMARY FILE with occupied pixels
  AllPixels      <- 1:length(mymap@data@values)
  AllPixelsModel <- cellFromXY(mymap, WBMatO[ , c('XCoord', 'YCoord')])
  MySummary      <- matrix(nrow = length(AllPixels), ncol= 2+length(mytimeafterinfHM))
  MySummary[,1]  <- AllPixels
  MySummary[AllPixelsModel, 2] <- WBMatO[,1]
  
    for (mytimeafterinf in  mytimeafterinfHM){
    pixelsD1 <- list()
    for (i in 1:maxit){
      if (mode == "noncumulative") {
        dayskipbase <- ifelse (mytimeafterinf == 1, 1, dayskip[mytimeafterinf-1, i])
        pixelsD1[[i]] <- YearPix[[i]][dayskipbase:dayskip[mytimeafterinf,i]]
      }
      
      if (mode == "cumulative") {
        pixelsD1[[i]] <- YearPix[[i]][1:dayskip[mytimeafterinf,i]]
      }
    }
    
    pixelsD1 <- lapply(pixelsD1, FUN = function(x) unique(x))
    pixelsD1 <- unlist(pixelsD1)
    pixelsD1 <- pixelsD1[-which(pixelsD1 == 0)]
    mapcells <- cellFromXY(mymap, WBMatO[pixelsD1, c('XCoord', 'YCoord')])
    infPixels.alliter <- mapcells  
    infPixels.alliter <- data.frame(table(infPixels.alliter))
    infPixels.alliter$infPixels.alliter <- as.numeric(as.character(infPixels.alliter$infPixels.alliter))
    infPixels.alliter$Freq <- round((100*infPixels.alliter$Freq)/ maxit)  # 
    
    if(mytimeafterinf == 1) { startpixels <- infPixels.alliter[,1]}
    MySummary[infPixels.alliter$infPixels.alliter, (2+which(mytimeafterinfHM == mytimeafterinf))] <- infPixels.alliter$Freq
    infPixels.alliter <- infPixels.alliter[-which(infPixels.alliter[,1] %in% startpixels),]
    infPixels.alliter <- infPixels.alliter %>%
      mutate(FreqClass = case_when(Freq < 10 ~ 1,
                                   Freq < 40 ~ 2,
                                   Freq < 70 ~ 3,
                                   Freq <= 100 ~ 4) )
    
        colfunc     <- colorRampPalette(c("pink", "red", "royalblue"))
    mynewcolors <- colfunc(4)  
    
    # Code to enables map colors to be consistent
    tmp1 <- sort(unique(infPixels.alliter$FreqClass))
    infPixels.alliter$FreqClass2 <- infPixels.alliter$FreqClass
    for (i in 1:length(tmp1)){
      infPixels.alliter$FreqClass2 <- replace(infPixels.alliter$FreqClass2, infPixels.alliter$FreqClass2 == tmp1[i], i)
    }
    
    mymap@data@values <- original.mapValues
    
    for (i in 1: dim(infPixels.alliter)[1]){
      mymap@data@values[infPixels.alliter[i,1]] <- (5+infPixels.alliter$FreqClass2[i])
      mymap@data@values[startpixels] <- 5
    }
    tmp1 <- sort(unique(infPixels.alliter$FreqClass))
    mynewcolors2 <- mynewcolors[tmp1] 
    
    #if (1 %in% (mymap@data@values)){
    #hmp# plot(mymap, col=c('grey', 'yellow', 'green','dark green', 'black', mynewcolors2),
    #hmp#      main = paste0("+ ", round(mytimeafterinf-1), " year" ), legend=FALSE)
    #}
    # if (!1 %in% (mymap@data@values)){
    #   plot(mymap, col=c('orange', 'grey', mynewcolors2),
    #        main = paste0("+ ", round(mytimeafterinf-1), " year" ), legend=FALSE)
    # }
    
   }
  
  #plot(c(3, 10, 50, 100), rep(1,4), col=mynewcolors, pch=19, cex=3, ylab='', xlab="Proportion of iterations with wild boars",
  #hmp# plot(c(5, 25, 50, 85), rep(1,4), col=mynewcolors, pch=19, cex=3, ylab='', xlab="Proportion of iterations with wild boars",
  #hmp#     
  #hmp#     xaxt="n", yaxt="n", xlim=c(0,100))
  #xtick<-c(5, 25, 75)
    #hmp# xtick<-c(0,10, 40, 70,100)
    #hmp#  text(x=xtick,   par("usr")[3], 
    #hmp#       labels = xtick, xpd = TRUE, pos = 1) #, srt = 45, pos = 1, xpd = TRUE)
  
    #hmp# dev.off()
  
  MySummary <- data.frame(MySummary)
  names(MySummary) <- c("PixelsTif", "PixelsModel", sapply((mytimeafterinfHM-1), FUN = function(x) paste0("Year ", x)))
  write.csv2(MySummary, paste0(outputfolder, "PropPixelsOccupied_",mode, "-", runID, "-", outputindex,".csv"), row.names = F)
}




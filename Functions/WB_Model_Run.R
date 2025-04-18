######################################################################
## Wild boar spread version 2025, 
## Programmed by Xavier Simons, based on WB ASF model of Tariq Halasa

######################################################################

## Function to generate n PERT beta random variates with min a, mode l, and max
## a.
rpert <- function(n,a,l,b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
  }
  a+(b-a)*rbeta(n,v,w)
}

## Model function
## b.
WBModel <- function(  MaxIterations,             # NUmber of iterations
                      MaxDays,                   # Number of days to simulate (2 years)
                      Detailed = TRUE,           # TRUE for saving detailed output
                      runID,                     # model run ID
                      
                      # POPULATION INITIALIZATION
                      FemaleInitMin,             # Initialization of the number of females per group     
                      FemaleInitMax,
                      AgeProbab,                 # Initialization of age of adults (from 2 to 10)
                      FemPropGroup,              # % of adult females in a group  
                      MalPropGroup,              # % of adult males in a group
                      SubAPropGroup,             # % of subadults in a group
                      PigletPropGroup,           # % of piglets in a group
                      
                      # YEARLY SURVIVAL PROBABILITIES
                      SurvivalProbAdF,           # for adult Females
                      SurvivalProbAdM,           # for adult males       
                      SurvivalProbSAdF,          # for subadult Females
                      SurvivalProbSAdM,          # for subadult males
                      SurvivalProbPigF,          # for piglets females 
                      SurvivalProbPigM,          # for piglets males
                      
                      # REPRODUCTION
                      propAdultRepro,            # proporiton of adults females that will reproduce 
                      propSubadRepro,            # proporiton of subadults females that will reproduce
                      propPigletRepro,           # proporiton of piglets females that will reproduce
                      RepProbList,               # probability of birth by week of the year 
                      NumOfSprProbList,          # probability of the number of offspring per litter
                      
                      # Max caring capacity for the number of adult females/pixel by landscape category (1, 2 and 3). Min/mode/max used in pert function.
                      MinCap1,
                      ModeCap1,
                      MaxCap1,
                      MinCap2,
                      ModeCap2,
                      MaxCap2,
                      MinCap3,
                      ModeCap3,
                      MaxCap3,
                      
                      # DISPERSION
                      ProbSplitMSA = 1/(6*7),              # Probability of male splitting per group per day of the week (7days) of the total spliting period (6 weeks)
                      SettlingProb,              # Probability to settle for a male in a habitat cell if it was only suitable habitate (0.5) or suitable habitate with females (1.0)
                      PmoveCatFem,               # Probability of females settling in category 1 or 2 pixels
                      HabitatProb,               # The probability to select a new pixel while male walking (we model that they prefer to walk in a habitat cell) 
                      DirectionList =  list(c(5,6,7),c(7,8,12), c(10,11,12), c(5,9,10)), ## North, ## East, ## South, and ## West
                      FencePorosity,             # Fence porosity
                      
                      # File pathways
                      FileWildBoarMat,           # WBMat file, with info on pixels ID, location an neighboring pixels
                      FileStartPixels,           # file of occupied pixels at the initialization of the model
                      FileBarriersPixels,        # file with ID of barriers pixels
                      MyOutputFolder            # path to folder to save model outputs 
                      
){
  
  WBMat <- read.table(FileWildBoarMat, sep=";", header = T, dec=',')
  StartPixels <- read.csv2(FileStartPixels, sep=';')
  StartPixels <- StartPixels[,2]
  InitialSizeGroup <- length(StartPixels)
  FencedPix <- read.csv2(FileBarriersPixels)
  FencedPix <- FencedPix$PixelID
  
  # reorganise category for test, keep only 3 category : 1 = habitat, 2 = accessible, 3 = inaccessible
  WBMat$NewValue <- NA
  WBMat$NewValue[WBMat$value == 1]        <- "inaccessible"  ## inaccessible
  WBMat$NewValue[WBMat$value == 2]        <- "accessible"  ## accessible
  WBMat$NewValue[WBMat$value == 3]        <- "suitable"   
  WBMat$NewValue[WBMat$value == 4]        <- "habitat" ## habitat
  WBMat$NewValue <- factor(WBMat$NewValue, levels=c("habitat", "suitable","accessible","inaccessible"))
  WBMat$value <- as.numeric(WBMat$NewValue)
  WBMat <- WBMat[, -which(names(WBMat)=="NewValue")]
  
  
  WBMat <- cbind(WBMat,0)
  library(data.table)
  
  ## inititate the daily summary matrices and variables
  DayOutPopMat  <- matrix(0,ncol=MaxIterations,nrow=MaxDays)
  DayOutAniMat  <- matrix(0,ncol=MaxIterations,nrow=MaxDays)
  YearPixels  <- matrix(0,ncol=MaxIterations,nrow=(MaxDays*InitialSizeGroup*20/365))
  cumDeath             <- rep(0,MaxIterations)
  
  for(Iter in 1:MaxIterations){
    set.seed(Iter)
    tmp.YearPixels   <- c()
    ## Initialize the population matrix
    PopMatWB <- matrix(0,ncol=12,nrow=50000)
    ## Add the habitat max caring capacity (animal/pixel) by landscape category (1, 2 and 3)
    WBMat[WBMat[,4]==1,13] <- ceiling(rpert(sum(WBMat[,4]==1),MinCap1,ModeCap1,MaxCap1)) 
    WBMat[WBMat[,4]==2,13] <- ceiling(rpert(sum(WBMat[,4]==2),MinCap2,ModeCap2,MaxCap2)) 
    WBMat[WBMat[,4]==3,13] <- ceiling(rpert(sum(WBMat[,4]==3),MinCap3,ModeCap3,MaxCap3)) 
    
    ## depict the home pixels for the populations, 
    TMPHP <- (WBMat[,1] %in% StartPixels) 
    homePixelsAll   <- sample(WBMat[TMPHP,1],InitialSizeGroup)
    
    ### Initialize the wild boar groups. Structure is based on Merta et al. (2015)
    ## initial size matrix
    InitSizeMat <- matrix(0,ncol=12)
    for(i in 1:InitialSizeGroup){
      ## the initial distribution of the groups was based on Merta et al. (2015). They were calculated relative to the average
      ## female percentage. for instance males were average % males divided by average % females. 
      females       <- round(runif(1,FemaleInitMin, FemaleInitMax)) 
      males         <- ceiling(females*MalPropGroup/FemPropGroup)
      subAdults     <- round(females*SubAPropGroup/FemPropGroup)
      Piglets       <- round(females*PigletPropGroup/FemPropGroup)
      GroupID       <- rep(i,sum(females,males,subAdults,Piglets))
      Sex           <- c(rep(1,females),rep(0,males),rbinom(subAdults,1,prob=0.5),rbinom(Piglets,1,prob=0.5))
      AgeCat        <- c(rep(3,females),rep(3,males),rep(2,subAdults),rep(1,Piglets))
      Dam           <- c(rep(0,sum(females,males)),sample(1:females,subAdults,rep=T),sample(1:females,Piglets,rep=T))
      # females deliver between Jan and June, so the age of the piglets will be between 183 and 365 days.
      # the same thing for sub-adults but with 365 days extra.
      tmpAgeSubA    <- sample((183:365)+365,females,rep=T)
      tmpAgePig     <- sample((183:365),females,rep=T)
      Age           <- c(sample(2:10,size=sum(females,males),T,prob=AgeProbab)*365,tmpAgeSubA[Dam[(females+males+1):(females+males+subAdults)]],
                         tmpAgePig[Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)]])
      # Adjust the Dam number to fit the actual ID of the DAMs
      Dam[(females+males+1):(females+males+subAdults)] <- Dam[(females+males+1):(females+males+subAdults)] + (dim(InitSizeMat)[1]) - 1 # -1 because we have extra row at start
      Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)] <- Dam[(females+males+subAdults+1):(females+males+subAdults+Piglets)] + (dim(InitSizeMat)[1]) - 1
      Breed         <- rep(0,sum(females,males,subAdults,Piglets))
      HomePixel     <- rep(homePixelsAll[i],sum(females,males,subAdults,Piglets))
      CurrPixel     <- HomePixel
      infectStatus  <- rep(0,sum(females,males,subAdults,Piglets))  
      SplitStatus   <- rep(0,sum(females,males,subAdults,Piglets))
      SplitMale     <- rep(0,sum(females,males,subAdults,Piglets))
      IDs           <- (max(InitSizeMat[,1])+1):(max(InitSizeMat[,1])+sum(females,males,subAdults,Piglets))
      
      InitMatWBPop  <- cbind(IDs,GroupID,Sex,AgeCat,Age,Breed,HomePixel,CurrPixel,infectStatus,SplitStatus,Dam,SplitMale) 
      InitSizeMat   <- rbind(InitSizeMat,InitMatWBPop)
      
    }
    InitSizeMat <- InitSizeMat[-1,] ## this is just to remove the first line that includes zeros. the first line of zeros is necessary for initialization to keep track of IDs.
    PopMatWB[1:dim(InitSizeMat)[1],] <- InitSizeMat
    
    # A matrix to store male subadult groups to split.
    GroupsToSplit        <- matrix(numeric(0),ncol=4)  
    
    ## Store the groups that have split this year
    SplittedGroups  <- numeric(0)
    cumDeathPar     <- 0
    Year            <- 1
    gTime           <- 0
    Criteria        <- TRUE
    TMPOutYesterday <- FALSE
    OnlyOnce        <- TRUE             
    
    # START DAILY LOOP
    while((gTime < MaxDays)){   
      
      MatSizeExp <- sum(PopMatWB[,1]==0)
      if(MatSizeExp<10000){
        ExtraRows <- matrix(rep(0,10000),ncol=12)
        PopMatWB <- rbind(PopMatWB,ExtraRows)
      } 
      
      ###### IMPORTANT TO REMEMBER ################
      ## You need to remember to update the spliting status per year. 
      # maybe other variables have to be updated yearly, like mortality and reproduction and so on....
      ##################################################################################################
      
      gTime <- gTime + 1
      PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
      
      if(gTime > 365) Year <- ceiling(gTime/365)  ## Count the year 
      PopMatWB[PopMatWB[,5]>0,5] <- PopMatWB[PopMatWB[,5]>0,5] + 1
      PopMatWB[PopMatWB[,5]>=365,4]     <- 2
      PopMatWB[PopMatWB[,5]>=(365*2),4] <- 3
      
      ## initiate variables and matrix that have to be initiated on yearly basis
      if(gTime %in% c(1, (366 * 1:(MaxDays/365)))){
        
        ### Mortality probabilities   
        ProbMortFact <- rnorm(1,1,0.05)
        ProbMortAdF <- 1 - (ProbMortFact*SurvivalProbAdF)^(1/365)
        ProbMortAdM <- 1 - (ProbMortFact*SurvivalProbAdM)^(1/365)
        ProbMortAdM[ProbMortAdM<0] <- 0.0001
        ProbMortSAdF <- 1 - (ProbMortFact*SurvivalProbSAdF)^(1/365)
        ProbMortSAdM <- 1 - (ProbMortFact*SurvivalProbSAdM)^(1/365)
        ProbMortSAdM[ProbMortSAdM<0] <- 0.0001
        ProbMortPigF <- 1 - (ProbMortFact*SurvivalProbPigF)^(1/365)
        ProbMortPigF[ProbMortPigF<0] <- 0.0001
        ProbMortPigM <- 1 - (ProbMortFact*SurvivalProbPigM)^(1/365)
        ProbMortPigM[ProbMortPigM<0] <- 0.0001
        
        
        ## Reproduction   
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[,2]),function(x) 1:sum(PopMatWB[,2]==x)))  
        yo <- data.table(PopMatWB)
        myfemAd <- yo[ , .(myfemales = rep(sum( (V3==1)*(V4==3)), .N)), V2 ][,myfemales]  
        myfemSubad <- yo[ , .(myfemSubad = rep(sum( (V3==1)*(V4==2)), .N)), V2 ][,myfemSubad]
        myfemPiglet <- yo[ , .(myfemPiglet = rep(sum( (V3==1)*(V4==1)), .N)), V2 ][,myfemPiglet]
        NumTotRepro <- ceiling(propAdultRepro*myfemAd + propSubadRepro*myfemSubad + propPigletRepro*myfemPiglet)   
        rm(yo)
        
        # A specific proportion of female adults, subadiults and piglets are allowed to breed 
        AniAllowBreed   <- AnimalsWInG<=NumTotRepro & PopMatWB[,3]==1 
        PopMatWB[AniAllowBreed,6]    <- (sample(RepProbList$Week,sum(AniAllowBreed),rep=T,prob=RepProbList$Prob)*7) + ((Year-1)*365)
        
        # Reset Male Splitting
        PopMatWB[,12]          <- 0 
        IndexRem               <- which(PopMatWB[,2]%in%GroupsToSplit[,1])
        PopMatWB[IndexRem,]    <- 0
        PopMatWB               <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        cumDeathPar            <- cumDeathPar + length(IndexRem)
        GroupsToSplit          <- matrix(numeric(0),ncol=4) 
      }
      
      ### Model Mortality ### 
      AdultsFToDie    <- which(PopMatWB[,4]==3&PopMatWB[,3]==1)
      AdultsFToDie    <- AdultsFToDie[rbinom(length(AdultsFToDie),1,ProbMortAdF)==1]
      AdultsMToDie    <- which(PopMatWB[,4]==3&PopMatWB[,3]==0)
      AdultsMToDie    <- AdultsMToDie[rbinom(length(AdultsMToDie),1,ProbMortAdM)==1]
      
      SubAdultsFToDie <- which(PopMatWB[,4]==2&PopMatWB[,3]==1)
      SubAdultsFToDie <- SubAdultsFToDie[rbinom(length(SubAdultsFToDie),1,ProbMortSAdF)==1]
      SubAdultsMToDie <- which(PopMatWB[,4]==2&PopMatWB[,3]==0)
      SubAdultsMToDie <- SubAdultsMToDie[rbinom(length(SubAdultsMToDie),1,ProbMortSAdM)==1]
      
      PigsFToDie      <- which(PopMatWB[,4]==1&PopMatWB[,3]==1)
      PigsFToDie      <- PigsFToDie[rbinom(length(PigsFToDie),1,ProbMortPigF)==1]
      PigsMToDie      <- which(PopMatWB[,4]==1&PopMatWB[,3]==0)
      PigsMToDie      <- PigsMToDie[rbinom(length(PigsMToDie),1,ProbMortPigM)==1]
      
      TooOldAni       <- which(PopMatWB[,5]==(11*365))
      ToDieNormal <- unique(c(AdultsFToDie,AdultsMToDie,SubAdultsFToDie,SubAdultsMToDie,PigsFToDie,PigsMToDie,
                              TooOldAni))
      
      if(length(ToDieNormal)>0){
        PopMatWB[ToDieNormal,] <- 0
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        
        cumDeathPar <- cumDeathPar+length(ToDieNormal)
      }
      
      ### Model Reproduction
      ## on daily basis, from Jan to end June, check if there are animals to deliver and make them deliver.
      if(gTime%in%PopMatWB[,6]){
        DelIndex        <- which(PopMatWB[,6]%in%gTime)
        NumOfSpring     <- sample(NumOfSprProbList$Number,length(DelIndex),rep=T,prob=NumOfSprProbList$Prob)
        DelIndex        <- DelIndex[NumOfSpring>0]
        NumOfSpring     <- NumOfSpring[NumOfSpring>0]
        if(length(NumOfSpring)>0){
          IDNew              <- (max(PopMatWB[,1])+1):(max(PopMatWB[,1])+sum(NumOfSpring))
          IndexLocation      <- (sum(PopMatWB[,1] > 0)+1):(sum(PopMatWB[,1] > 0)+length(IDNew))
          PopMatWB[IndexLocation,1]  <- IDNew
          PopMatWB[IndexLocation,2]  <- rep(PopMatWB[DelIndex,2],NumOfSpring)
          PopMatWB[IndexLocation,3]  <- rbinom(sum(NumOfSpring),1,0.5)
          PopMatWB[IndexLocation,4]  <- 1
          PopMatWB[IndexLocation,5]  <- 1
          PopMatWB[IndexLocation,6]  <- 0
          PopMatWB[IndexLocation,7]  <- rep(PopMatWB[DelIndex,7],NumOfSpring)
          PopMatWB[IndexLocation,8]  <- rep(PopMatWB[DelIndex,7],NumOfSpring)
          PopMatWB[IndexLocation,9]  <- 0
          PopMatWB[IndexLocation,10] <- 0
          PopMatWB[IndexLocation,11] <- rep(PopMatWB[DelIndex,1],NumOfSpring)
          PopMatWB[IndexLocation,12] <- 0
          
          PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
          
        }
      }
      
      ### Model females groups splitting
      ## This code checks whether splitting of subadult females may happen. it occurs only once in week 28 Kramer-Schadt et al. (2009)
      if(gTime%in%((28*7+((Year-1)*365))+0:6)){
        
        mytmporder   <- PopMatWB[,4]
        mytmporder[mytmporder == 1] <- 11 ## piglets
        mytmporder[mytmporder == 3] <- 13 ## adults  ## mytmp order enables to put order 1°) adults, 2°) piglets, 3°) subaults. ==> Breedcap is a max capa for females (adults, subadults and piglets included!)
        PopMatWB        <- PopMatWB[order(PopMatWB[,2],PopMatWB[,3],mytmporder,decreasing=T),]  
        AnimalsWInG     <- unlist(sapply(unique(PopMatWB[,2]),function(x) 1:sum(PopMatWB[,2]==x)))
        BreedCapCells   <- c(WBMat[PopMatWB[,7],13],rep(0,sum(PopMatWB[,7]==0)))   
        ## Identify the animals within the groups that will split. females subadults > breading capacity, did not split before 
        GroupSplitNum   <- cbind(sort(unique(PopMatWB[,2])),(tapply((AnimalsWInG>BreedCapCells & PopMatWB[,4]==2 & PopMatWB[,3]==1 & PopMatWB[,10]==0),  PopMatWB[,2], sum)))
        GroupSplitNum   <- GroupSplitNum[GroupSplitNum[,2]>=2,,drop=FALSE]  ### !!! need to be at least 2 females subadults to split
        ## Idetify groups where there are only males. to allow females to join groups where there are only males.
        PixelMalesOnly  <- cbind(sort(unique(PopMatWB[,2])),(tapply((PopMatWB[,3]==0),PopMatWB[,2], all)))
        PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly[,2]==1,,drop=FALSE]
        
        GroupsSplitShD  <- sapply(PopMatWB[match(GroupSplitNum[,1],PopMatWB[,2]),7],function(x) {
          Distance <- sqrt((WBMat[x,2]-WBMat[,2])^2 + (WBMat[x,3]-WBMat[,3])^2)/1000
          tmp1     <- which(Distance<=9)
          any(WBMat[tmp1,4]==1 & !(tmp1%in%PopMatWB[!PopMatWB[,2]%in%PixelMalesOnly[,1],7]))})
        GroupSplitNumSD <- cbind(GroupSplitNum,GroupsSplitShD)
        GroupSplitNumSD <- GroupSplitNumSD[GroupSplitNumSD[,3]==1,,drop=FALSE]
        ## Splitting status
        PopMatWB[PopMatWB[,2]%in%GroupSplitNumSD[,1],10] <- 1
        ## Make short term split to happen
        if(dim(GroupSplitNumSD)[1]>0){
          OriginPixel     <- PopMatWB[match(GroupSplitNumSD[,1],PopMatWB[,2]),7]
          
          TargetPixel <- numeric(0)
          PixelMalesOnly  <- unique(PopMatWB[PopMatWB[,2]%in%PixelMalesOnly[,1],7])
          PixelMalesOnly  <- PixelMalesOnly[PixelMalesOnly>0]
          for(xx in OriginPixel){
            Distance <- sqrt((WBMat[xx,2]-WBMat[,2])^2 + (WBMat[xx,3]-WBMat[,3])^2)/1000
            tmp1      <- which(Distance<=9)
            # The pixel that is <9km should be 1) of category at least suitable AND not occupied OR 2) with only males AND at least suitable
            tmp2      <- tmp1[((WBMat[tmp1,4] %in% c(1,2)) & !(tmp1%in%PopMatWB[,7])) | (tmp1 %in% PixelMalesOnly & (WBMat[tmp1,4] %in% c(1,2)))]
            tmpcat <- WBMat[tmp2,4]
            
            if(length(tmp2)>1)  tmp3 <- sample(tmp2,1, prob = PmoveCatFem[tmpcat])
            if(length(tmp2)==1) tmp3 <- tmp2
            if(length(tmp2)==0) tmp3 <- 0
            TargetPixel <- c(TargetPixel,tmp3)
          }
          GroupSplitNumSD <- GroupSplitNumSD[TargetPixel>0,,drop=FALSE]
          OriginPixel     <- OriginPixel[TargetPixel>0]
          TargetPixel     <- TargetPixel[TargetPixel>0]
          GroupSplitNumSD <- cbind(GroupSplitNumSD,TargetPixel) 
          
          ## This part of the code to allow movement
          if(dim(GroupSplitNumSD)[1]>0){
            # A list to keep track for the pixels where the pigs have been. Make sure that this information is exported on daily basis, because the 
            # list will be re-initiate every day splitting may happen.
            PixelsMoved <- as.list(matrix(0,ncol=length(TargetPixel)))
            CurrentPos <- OriginPixel
            for(i in 1:length(TargetPixel)){
              Trail <- 0
              prevEdge <- CurrentPos[i]
              while(CurrentPos[i]!=TargetPixel[i] & Trail < 10){
                Trail      <- Trail + 1
                Edges      <- unlist(WBMat[CurrentPos[i],5:12])
                Edges      <- Edges[Edges>0&Edges!=prevEdge]
                EdgesCat <- WBMat[Edges, 4]
                Fenced   <-  Edges %in% FencedPix 
                Fenced   <- rbinom(length(Fenced), Fenced, (1-FencePorosity))  
                Edges    <- Edges[EdgesCat %in% 1:3 & Fenced == 0]    
                
                if(length(Edges)>0){
                  DistEdges  <- sqrt((WBMat[TargetPixel[i],2]-WBMat[Edges,2])^2 + (WBMat[TargetPixel[i],3]-WBMat[Edges,3])^2)/1000
                  prevEdge   <- CurrentPos[i]
                  NewPosition<- Edges[DistEdges==min(DistEdges)]
                  if(length(NewPosition)>1) NewPosition<-sample(NewPosition,1)
                  CurrentPos[i] <- NewPosition
                  # Here we keep track of the pixels where the pigs moved to.
                  PixelsMoved[[i]]<- c(PixelsMoved[[i]],NewPosition)
                }
              }
            }
            names(PixelsMoved) <- GroupSplitNumSD[,1]
            ## Here we assume that groups that did not find the way, did not actually split as the edges was not connecting.
            IndexNotSplit   <- which(CurrentPos!=TargetPixel)
            if (length(IndexNotSplit) > 0) {GroupSplitNumSD <- GroupSplitNumSD[-IndexNotSplit,,drop=FALSE]}
            ## Make females from different groups moving to a new pixel to form a new group.
            if(dim(GroupSplitNumSD)[1]>0){
              newGroupIDs         <- (max(PopMatWB[,2])+1):(max(PopMatWB[,2])+dim(GroupSplitNumSD)[1])
              if(sum(duplicated(GroupSplitNumSD[,4]))>0){
                for(l in unique(GroupSplitNumSD[,4])){
                  TEMP  <- which(GroupSplitNumSD[,4]==l)
                  if(length(TEMP)>1){
                    newGroupIDs[TEMP] <- newGroupIDs[TEMP[1]]
                  }
                }
              }
              
              for(b in 1:dim(GroupSplitNumSD)[1]){
                ## if there is a male group already in the pixel and is not splitting, then make them one group with the new females.
                ## Notice the code above allows them to come into a new pixel only if it is empty or there are only males. so dont worry tat IndexMalComb1
                ## is matching with all PopMatWB[,7]
                IndexMalComb1  <- GroupSplitNumSD[b,4]%in%PopMatWB[,7]
                if(IndexMalComb1){
                  IndexMalComb2  <- unique(PopMatWB[PopMatWB[,7]==GroupSplitNumSD[b,4],2])%in%GroupsToSplit[,1]
                  if (length(IndexMalComb2) > 1) {IndexMalComb2 <- FALSE}
                  if(!IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2  & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                    MalesInPixel     <- PopMatWB[,3]==0  & PopMatWB[,7] == GroupSplitNumSD[b,4] & !PopMatWB[,2]%in% GroupsToSplit[,1]
                    newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)+sum(MalesInPixel)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF)+sum(MalesInPixel))
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF|MalesInPixel,8]  <- NewHomePixel
                  }
                  if(IndexMalComb2){
                    AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2  & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                    newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)) 
                    NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF))
                    PopMatWB[AnimalsCanSplitF,2]  <- newGroupIDsAn
                    PopMatWB[AnimalsCanSplitF,7]  <- NewHomePixel
                    PopMatWB[AnimalsCanSplitF,8]  <- NewHomePixel
                  }   
                }
                if(!IndexMalComb1){
                  AnimalsCanSplitF <- AnimalsWInG>BreedCapCells & PopMatWB[,3]==1 & PopMatWB[,4]==2 & PopMatWB[,2]%in%GroupSplitNumSD[b,1]
                  newGroupIDsAn    <- rep(newGroupIDs[b],sum(AnimalsCanSplitF)) 
                  NewHomePixel     <- rep(GroupSplitNumSD[b,4],sum(AnimalsCanSplitF))
                  PopMatWB[AnimalsCanSplitF,2]  <- newGroupIDsAn
                  PopMatWB[AnimalsCanSplitF,7]  <- NewHomePixel
                  PopMatWB[AnimalsCanSplitF,8]  <- NewHomePixel
                }   
              }
            }
          } 
        }
      }#Closes short distance splitting
      
      ## after the end of female splitting, we reset female splitting.
      if(gTime%in%((28*7+((Year-1)*365))+7)){
        PopMatWB[,10]    <- 0
        SplittedGroups   <- numeric(0)
        PixelsMoved      <- as.list(matrix(0,ncol=1))
      }
     
      ##########################################
      ######################
      #Model Male splitting
      ##############################
      # First we determine the period of splitting of males, as defined in Lange et al., (2012). Males may find a pixel to live in otherwise they will keep
      # wondering around until they either die or find a place to live in  (fin juin - fin juillet)
      if(gTime%in%((25*7+((Year-1)*365)):(30*7+((Year-1)*365)))){
        ## sort the matrix
        PopMatWB    <-  PopMatWB[order(PopMatWB[,2],PopMatWB[,3],PopMatWB[,5],decreasing=T),]
        ## define the males as adults and subadults
        # Identify the groups where subadult males may split. 
        ## Calculate nuber of adult males
        NumbAdultMInG <-  cbind(sort(unique(PopMatWB[,2])),tapply(PopMatWB[,3]==0&PopMatWB[,4]==3,PopMatWB[,2],sum))
        ## No need for subadult mails to split if there is need for them in the group
        ## select groups in which there are adults males and extract their number, with same legnth as PopMatWB 
        tmpIndexSG    <- NumbAdultMInG[match(PopMatWB[,2],NumbAdultMInG[,1]),2]
        # Groups that may split must not be a new splilt groups and have more than 2 adult males
        GroupSubAM    <- unique(PopMatWB[!(PopMatWB[,2]%in%SplittedGroups) & PopMatWB[,3]==0 & PopMatWB[,4]==2 & PopMatWB[,12]==0 & tmpIndexSG>2,2])
        if(length(GroupSubAM)>0) 
          GroupSubAM <- GroupSubAM[rbinom(length(GroupSubAM),1,ProbSplitMSA)==1]
        if(length(GroupSubAM)>0){
          # Identify the subadault males that may split. We must ensure that the groups have adult males before the subadults are split.
          # no need for the subadults to find new areas if they can nurish and breed in their home
          tmpSAM <- which(PopMatWB[,2]%in%GroupSubAM & PopMatWB[,3]==0 & PopMatWB[,4]==2 )# & TEMPMat[match(PopMatWB[,2],TEMPMat[,1]),2]>=2)
          # check which ones that would split
          if(length(tmpSAM)>0) {
            ## make subadults to split based on probability from 25 to 75% (Truve et al., 2004).
            tmpSAM           <-tmpSAM[rbinom(length(tmpSAM),1,prob=runif(length(tmpSAM),0.25,0.75))==1]
            AinmFromGToSplit <- cbind(sort(unique(PopMatWB[tmpSAM,2])),tapply(tmpSAM,PopMatWB[tmpSAM,2],length))
            AinmFromGToSplit <- AinmFromGToSplit[AinmFromGToSplit[,2]>1,,drop=FALSE]
            tmpSAM           <-tmpSAM[PopMatWB[tmpSAM,2] %in% AinmFromGToSplit[,1]]   ## here tmpsam is the index of the animals that will split (subadults males)
            
            if(length(tmpSAM)>0) {# do not worry, if the dimension is > 0 then there will be at least 2 animals from one group to split.
              PopMatWB[tmpSAM,12] <- 1
              ## Store the group number of the groups that have splitted
              SplittedGroups <- c(SplittedGroups,unique(PopMatWB[tmpSAM,2]))
              # number of spliting males per group
              AinmFromGToSplit   <- cbind(AinmFromGToSplit,   sample(1:4, dim(AinmFromGToSplit)[1], replace=T),    
                                          sapply(AinmFromGToSplit[,1],  function(x) 
                                            PopMatWB[PopMatWB[,2] == x,c(8)][1]))    
              
              # the third column, informs about the direction
              # which should be dipicted from the DirectionList. each group would walk in a specific direction
              # if they cannot find a connecting edge in that direction, a random cell is then selected.
              # make them to prefer to move through habitat cells rather than fields.
              indexsNGroups  <- (max(PopMatWB[,2])+1):(max(PopMatWB[,2])+dim(AinmFromGToSplit)[1])
              indexNGPigs    <- rep(indexsNGroups,AinmFromGToSplit[,2])
              indexNewGNum   <- sort(which(PopMatWB[,2]%in%AinmFromGToSplit[,1] & PopMatWB[,12] == 1),decreasing = T)
              PopMatWB[indexNewGNum,2] <- indexNGPigs
              AinmFromGToSplit[,1] <- indexsNGroups
              GroupsToSplit <-rbind(GroupsToSplit,AinmFromGToSplit)      
            }
          }
        }
      }
      
      if(dim(GroupsToSplit)[1]>0){
        ## First we check that all groups still exist in the population
        Checktmp         <- which(!GroupsToSplit[,1]%in%PopMatWB[,2])
        if(length(Checktmp)>0) GroupsToSplit[-Checktmp,,drop=FALSE]
        ## determine the number of pixels they may move per day and the direction
        NumbPixMovesTod  <- round(rpert(dim(GroupsToSplit)[1],0,2,4)) # Lange et al (2015) on average 4 km/day
        MovedPixMale <- as.list(matrix(0,ncol=length(NumbPixMovesTod>0)))
        ToRemovejj <- numeric(0)
        if(any(NumbPixMovesTod>0)){
          for(i in 1:max(NumbPixMovesTod)){
            if(any(PopMatWB[,12]==1) & any(NumbPixMovesTod>=i)){
              IndexGroup   <- which(NumbPixMovesTod>=i)
              if(length(IndexGroup)>0){
                for(jj in IndexGroup){
                  tmp1 <- unlist(WBMat[GroupsToSplit[jj,4],DirectionList[[GroupsToSplit[jj,3]]]])
                  tmp1 <- tmp1[tmp1>0]
                  tmp2 <- WBMat[tmp1,4]
                  Fenced     <- tmp1 %in% FencedPix  
                  Fenced     <- rbinom(length(Fenced), Fenced, (1-FencePorosity))  
                  tmp3 <- which(tmp2>0 & tmp2<4 & Fenced ==0)  
                  tmp4 <- tmp1[tmp3]
                  if(length(tmp4)>0){
                    if(length(tmp4)==1) tmp5 <- tmp4
                    if(length(tmp4)> 1) tmp5 <- sample(tmp4,1,prob=HabitatProb[tmp2[tmp3]])
                  }
                  if(length(tmp4)==0){
                    Tmp1 <- unlist(WBMat[GroupsToSplit[jj,4],5:12])
                    Tmp1 <- Tmp1[Tmp1>0]
                    Tmp2 <- WBMat[Tmp1,4]
                    Fenced     <- Tmp1 %in% FencedPix  
                    Fenced     <- rbinom(length(Fenced), Fenced, (1-FencePorosity)) 
                    Tmp3 <- which(Tmp2>0 & Tmp2<4 & Fenced ==0)
                    Tmp4 <- Tmp1[Tmp3]
                    if(length(Tmp4)>1)  tmp5 <- sample(Tmp4,1,prob=HabitatProb[Tmp2[Tmp3]])
                    if(length(Tmp4)==1) tmp5 <- Tmp4
                  }
                  #### here we add the list with the moved pixels
                  MovedPixMale[[jj]] <- c(MovedPixMale[[jj]],tmp5)
                  # is this only a suitable habitat pixel or a suitable habitat and has a female(s) with no or 1 male
                  
                  TMPIndPix <- (WBMat[tmp5,4] %in% c(1,2) & sum(PopMatWB[PopMatWB[,7]==tmp5,3]==1)>0  &
                                  sum(PopMatWB[PopMatWB[,7]==tmp5,3]==0)<2) 
                  # the first part of the c() represent a good habitat cell and not occupied
                  tmp6      <- which(c(WBMat[tmp5,4] %in% c(1,2) & !(tmp5%in%PopMatWB[,7]), TMPIndPix))
                  
                  # we decide whether they will settle in this pixel or not based on a random process
                  # Notice here that the animals that die during splitting do not affect the splitting. We checked above that all groups in splitting martix have animals
                  # in the Population matrix
                  if(length(tmp6)>0){
                    tmp7 <- rbinom(1,1,prob=SettlingProb[max(tmp6)])==1
                    if(tmp7&!TMPIndPix){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1)
                      PopMatWB[indexspPigs,7] <- tmp5
                      PopMatWB[indexspPigs,8] <- tmp5
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    
                    ## If the group will settle with a pre-existing female group, then their group number will be the same as the females
                    if(tmp7&TMPIndPix){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1)
                      IndexSelPigs            <- which(PopMatWB[,7]%in%tmp5 & !PopMatWB[,2]%in%GroupsToSplit[jj,1] & !PopMatWB[,2]%in%GroupsToSplit[GroupsToSplit[,4]%in%PopMatWB[,7],1])
                      newGroupIDM             <- rep(PopMatWB[IndexSelPigs,2][1],length(indexspPigs))  
                      PopMatWB[indexspPigs,2] <- newGroupIDM
                      PopMatWB[indexspPigs,7] <- tmp5
                      PopMatWB[indexspPigs,8] <- tmp5
                      PopMatWB[indexspPigs,12]<- 0
                      ToRemovejj              <- c(ToRemovejj,jj)
                    }
                    if(!tmp7){
                      indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 )
                      PopMatWB[indexspPigs,8] <- tmp5
                      GroupsToSplit[jj,4]     <- tmp5
                    }
                  }
                  
                  if(length(tmp6)==0) {
                    indexspPigs             <- which(PopMatWB[,2]==GroupsToSplit[jj,1] & PopMatWB[,12]==1 )
                    PopMatWB[indexspPigs,7] <- tmp5
                    PopMatWB[indexspPigs,8] <- tmp5 ## we set the home pixel as the current pixel for splitting groups as they do not have home yet
                    GroupsToSplit[jj,4]     <- tmp5
                  }
                  
                }
              }
            }
          }
          names(MovedPixMale) <- GroupsToSplit[,1]
          if(length(ToRemovejj)>0) GroupsToSplit  <- GroupsToSplit[-ToRemovejj,,drop=FALSE]
        }
      }
      
      #######################################################################
      # ASF MODULE was here
      #####################################################################
      
      ## Here we re-initiate the daily matrix for males and females movements to insure that there is no bleed from previous days
      MovedPixMale <- as.list(matrix(0,ncol=1))
      PixelsMoved <- as.list(matrix(0,ncol=1))
      
      ## Make the daily summary
      DayOutPopMat[gTime,Iter] <- length(unique(PopMatWB[,2]))-1
      DayOutAniMat[gTime,Iter] <- sum(PopMatWB[,2]>0)
      
      ## Make yearly summary
      if (gTime %in% c(1, (365 * 1:(MaxDays/365)))){
        tmp.YearPixels <- c(tmp.YearPixels, unique(PopMatWB[,8]))
      }
      
      
      print(c(Iter,gTime))  
    }#while(gTime
    
    ## Make the summaries per iteration
    YearPixels[1:length(tmp.YearPixels),Iter]  <- tmp.YearPixels
    
  }#Close for(Iter in...)
  
  NAMETG    <- paste0(MyOutputFolder, runID,"-DayOutPopMat.txt")
  NAMETA    <- paste0(MyOutputFolder, runID,"-DayOutAniMat.txt")
  NAMETPY    <- paste0(MyOutputFolder, runID,"-YearPixels.txt")
  
  write.table(DayOutPopMat,NAMETG,sep=" ",col.names = F,row.names=F)
  write.table(DayOutAniMat,NAMETA,sep=" ",col.names = F,row.names=F)
  write.table(YearPixels,NAMETPY,sep=" ",col.names = F,row.names=F)
 
}


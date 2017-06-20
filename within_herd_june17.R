#Creation of own rpert function
rpert <- function(n,a,l,b) {
  mu <- (a+4*l+b)/6
  if (mu==l) v <- w <- 3 else {
    v <- (mu-a)*(2*l-a-b)/(l-mu)/(b-a)
    w <- v*(b-mu)/(mu-a)
  }
  a+(b-a)*rbeta(n,v,w)
}

#Building a function to generate the herd
GenerateAHerd <- function(NAnimals=7266){ #Number of animals set to 7266


Bes <- data.frame(Animal=1:NAnimals,Infected=FALSE,Status=1) #Create a data frame with this no. of animals


Lit_Table <- c(358 , 505 , 652 , 799 , 946 ,1093 ,1240 ,1387) #Vector for 8 litters (sow age in days upon farrowing)

#Introduction of different animal agegroups/types ("Type"). The categories are: 1 = Sows, 2 = WF, 3 = G, 4 = Boars
Type=rep_len(rep(c(1:4), c(23,322,1,0)),length(Bes$Animal)) #Type==4 (boars) is currently not included
Bes$Type=rep_len(rep(c(1:4), c(23,322,1,0)),length(Bes$Animal)) 

#Agedistribution in days (Dist) within the different agegroups
Dist <- c('rep_len(rep(seq(from=239,to=1388, by=7), each=23),sum(Bes$Type==1))','rep_len(rep(seq (from=1, to=170, by=7), each=322),sum(Bes$Type==2))','239','150')
#Distribution of batch numbers within the different agegroups. Set to one fot all gilts and boars, since batch production is not relevant for these
Dist2 <- c('rep_len(rep(seq(from=1,to=21, by=1), each=23),sum(Bes$Type==1))','rep_len(rep(seq (from=1, to=1000, by=1), each=12),sum(Bes$Type==2))', 'rep_len(rep(seq (1, 1000, by=1), each=5),sum(Bes$Type==3))','1')

#Assign the ages and batch numbers defined above to the animals
Bes$Age <- 0 
for(i in unique(Bes$Type)){
  Index <- Bes$Type==i
  Bes$Age[Index] <- eval(parse(text=rep(Dist[i],sum(Index))))
  Bes$teamo[Index] <- eval(parse(text=rep(Dist2[i],sum(Index))))
}

#Apply mixed sow age within the batches
for (teamo in 1:21) {
  teamdiff <- rep_len(7*(teamo-1), 23)
  Sowage <- sample(seq(240, 1269, by=147), 23, rep=T, prob=c(0.22, 0.20, 0.17, 0.14, 0.11, 0.085, 0.05, 0.025)) + teamdiff
  Bes$Age[Bes$teamo==teamo&Bes$Type==1] <- Sowage
}

#Define and assign slaughterage to the slaughterpigs
sage <- c(rep(165, times=40),rep(c(158, 172), times=8), rep(c(151, 179), times=2)) #Andel restgrise sat til ca. 20%
Bes$Slaughterage <- ifelse(Bes$Type==2, sample(sage, replace=TRUE, nrow(Bes)),0)
                                                 
#Define housing unit (stable)  based on age on a given day
Tmate <-c(240:244) #First stay in the mating unit - time per cycle = 5 days, starting from 8 months of age. Sows/gilts are moved to the gestation unit 4 weeks after mating
Scycmat<-rep(0:7, each=length(Tmate))#Sow cycle number, repeats depend on length of stay in the relevant unit
Sclength<-147 #Length of a full sow cycle (No. of dsyd between the sow is returning to a given stable unit)
Mate<-Tmate+Scycmat*Sclength #Vector including the days in a sows' life spend in the mating unit (based on the sows (adjusted) age in days)

Tgest <- c(245:353) #First stay in the gestation unit - it is assumed that the duration of the gestation period is 114 days, and that the sow is moved to the farrowing stable, approximately 5 days before farrowing. Time per cycle in thsi unit = 109 days
Scycges<-rep(0:7, each=length(Tgest)) #Sow cycle number, repeats depend on length of stay in the relevant unit
Gest<-Tgest+Scycges*Sclength #Vector including the days in a sows' life spend in the gestation unit (based on the sows (adjusted) age in days)

TFarr <-c(354:386) #First stay in the farrowing unit - time per cycle = 33 days
Scycfar<-rep(0:7, each=length(TFarr)) #Sow cycle number, repeats depend on length of stay in the relevant unit
Farr<-TFarr+Scycfar*Sclength #Vector including the days in a sows' life spend in the farrowing unit (based on the sows (adjusted) age in days)

#Extension of the stay in the mating unit to include the first four weeks of the gestation period
Lit_no <- c(1:8) #for a given sow (no to be confused with the global litter no: "litter") #This variable is also used in rel. to culling
Farrowdays <- (Lit_no-1)*Sclength+358 #Day of farrowing
Matdays <- Farrowdays - 114 #Days of insemination
Matxtra <- rep(Matdays, each=28) + 1:28
Mate <- c(Mate, Matxtra)
Gest <- Gest[!(Gest%in%Mate)] 

#Division into stables (age-based)
Bes$stable <- ifelse(Bes$Type==2&Bes$Age>=80, 5, ifelse(Bes$Type==2&Bes$Age<80&Bes$Age>28, 4, 
                                                        ifelse(Bes$Type==2&Bes$Age<29, 3, ifelse(Bes$Type==1&Bes$Age%in%Gest, 2, ifelse(Bes$Type==1&Bes$Age%in%Farr, 3, 1)))))
#Categories for Stable 1 = Mate, 2 = Gest, 3 = Farr, 4 = Wean, 5 = Fini.

#In the beginning litternumber + pen no. etc are set to zero
Bes$litter <- 0
Bes$dam[Bes$Type %in% c(1,3)]    <- Bes$Animal[Bes$Type %in% c(1,3)]
Bes$Sti    <- 0
Bes$LCount <- 0
Bes$LCount[Bes$Type==1] <- sapply(Bes$Age[Bes$Type==1],function(x) sum(x>=Lit_Table))

#To balance the herd at start - 3 extra finisher batches are added
for(i in c(1:3)) {
  tmp10  <- as.numeric(max(Bes$Animal)+1:322)
  tmp11 <- as.data.frame(cbind(tmp10,FALSE,1, 2,148+((i-1)*7),max(Bes$teamo+1),0,5,0,9999,0,0)) 
  names(tmp11) <- names(Bes)
  Bes <- rbind(Bes,tmp11)
}

Slaage <- sample(sage, replace=TRUE, 644)
Bes$Slaughterage[Bes$teamo %in% c(max(Bes$teamo), (max(Bes$teamo)-1))] <- Slaage

if (any(Bes$Age>Bes$Slaughterage)) {
  rem1 <- Bes$Animal[Bes$Type==2&Bes$Age>Bes$Slaughterage]
  if (length(rem1)>0) {
    Bes <- Bes[-which(Bes$Animal%in%rem1),] 
  }
} 

Bes$dam[Bes$Type %in% c(1,3)]    <- Bes$Animal[Bes$Type %in% c(1,3)]
Bes$dam[Bes$Type==2 & Bes$Age>28]    <- 9999 #No need for assigning a dam to these, 9999 is used as "dummy" number

#Assignment of dam id to the piglets (which sow, were they born by)
for (i in 0:3) {
  sid <- Bes$Animal[Bes$Type==1 & Bes$Age %in% (Farrowdays+1+(7*i))]
  Bes$dam[Bes$Type==2&Bes$Age==((7*i)+1)] <- rep(sid, each=14) #14 is calibrated to match the 322 animals
 }

sid <- Bes$Animal[Bes$Type==1 & Bes$Age %in% c(Farrowdays+29, 240)] #Needs to be different from above, as it bspecifically neeeds to match the 240 olds
Bes$dam[Bes$Type==2&Bes$Age==29] <- rep(sid, each=14)

#Addition of extra pigs to team 1 to balance batch size (only needed at simulation initiation, because it is the first batch - later batches will receive sows that have undergone re-insemination)
for(i in c(1:3)) {
  tmp8  <- as.numeric(max(Bes$Animal+1))
  tmp9 <- as.data.frame(cbind(tmp8,FALSE,1, 1,240+147*(i-1),1,0,1,0,tmp8,0,(i-1))) 
  names(tmp9) <- names(Bes)
  Bes <- rbind(Bes,tmp9)
}

return(Bes)
 }#End of GenerateAHerd



##################
### UpdateHerd ###
##################
#updateHerd <- function(runID="MRSA",MaxIterations=1,Days=26,DayOption=1,seed=1){
MaxIterations <- 25 #No. of interations to run
Days          <- 3650 #No. of days to run
DayOption     <- 1
seed          <- 10 #Specify seed (if want to run with this)
runID         <- "Test_output" #Name to identify output file
DiseaseStart  <- 1460 #Day after simualation start, where the infection is first introduced

set.seed(seed)

#Prepare output file
Outputherd <- matrix(numeric(0),ncol=65)
NAME <- paste(runID,"MRSA.txt",sep="-") 
write.table(Outputherd,NAME,sep=" ")

Bes <- GenerateAHerd()

for(iter in 1:MaxIterations){
    
  
  print(iteration) #Print iteration no. if want to follow run

  #set various baseline values to be used later
  gTime=0
  TBes<-Bes
  iteration  <- iter
  basefrss   <- NULL
  amsti      <- NULL
  d          <- as.data.frame(NULL)
  fasec      <- 0
  wsec       <- 7
  finsec     <- 0
  TBes$nsr   <- 0
  TBes$ns    <- 0
  TBes$nsow  <- TBes$dam
  TBes$nsd   <- 0
  TBes$ri    <- 0
  TBes$rid   <- 0
   
  TBes$teama <-TBes$teamo #Adjusted batch no is equal to the original batch no. at simulation start
  TBes$PACount <-TBes$Age #PACount = adjusted age counter, wheres Age = chronological age counter
  
  TBes$snoutcont  <- 0
  
  TBes$ProbInfInd <- 0 #Probability of infection
  TBes$dcd        <- 0 #Duration of carriage
  TBes$Load       <- 0 #
  
  TBes$sec        <- 0 #Housing section set to zero here, but specified later
  TBes$sti        <- 0
  
  h <- seq(0,3650, by=14) #Will need extension if the model should ever need to run for more than 10 years

  secfrach      <-0.20 #Fraction of between pen spread, the between section spread will be equal to in workintensive units
  secfrac       <- 0.15 #As above for other units
  stafrac       <- 0.02 #Fraction of between pen spread, the between stables spread will be eqaul to

#Calibration factors (index in the name indicates stable followed by type, i.e. 11=sows in the mating stable (animal type 1 in stable unit 1))

  kal11 <- 1
  kal13 <- 1
  kal21 <- 1
  kal31 <- 0.5
  kal32 <- 0.5
  kal42 <- 1.5
  kal52 <- 1 

#Transmission parameters with use of risk antimicrobials - inputs for distributions
  BetaWDW      <-0.17931 #Based on Broens et al., 2012a: Value for post-weaning, ab:+ and pIP=1 and Broens, 2011 (Duration of infection=17.4 days)
  BetaWDWMin   <-0.06609
  BetaWDWMax   <-0.48621
  BetaWDF      <-0.17931 #Currently the same as for weaners
  BetaWDFMin   <-0.06609
  BetaWDFMax   <-0.48621
  BetaWDO      <-0.17931 #Within pen for others (gilts and pregnant sows sharing pen in the gestation stable). Currently the same as for finishers
  BetaWDOMin   <-0.06609
  BetaWDOMax   <-0.48621
  BetaBPDW     <-0.03448 #The value from Broens for post-weaning, ab:+ and pIP=0 
  BetaBPDWMin  <-0.01954
  BetaBPDWMax  <-0.06092
  BetaBPDF     <-0.03448 #Same as for weaners
  BetaBPDFMin  <-0.01954
  BetaBPDFMax  <-0.06092
  BetaBPDO     <-0.03448 #This value needs to be somewhere between BPDF and BPDW (currently the same value). 
  BetaBPDOMin  <-0.01954
  BetaBPDOMax  <-0.06092
  BetaBPDX     <-0.08966 #Pre-weaning values 
  BetaBPDXMin  <-0.03736
  BetaBPDXMax  <-0.21667
  BetaSOD      <-0.46437 #Currently the same as between piglets
  BetaSOMin    <-0.12471 #Currently the same as between piglets
  BetaSOMax    <-1.73103 #Currently the same as between piglets
  BetaPID      <-0.46437 #Broens value for pre-weaning pigs, ab:+, pIP=1
  BetaPIMin    <-0.12471  
  BetaPIMax    <-1.73103 

#Transmission paramters with use of risk antimicrobials
# BetaWDW      <-0.0711 #%AB, Based on Broens et al., 2012a: Value for post-weaning, ab:% and pIP=1 and Broens, 2011 (Duration of infection=17.4 days)
# BetaWDWMin   <-0.0345 #%AB   
# BetaWDWMax   <-0.1425 #%AB
# BetaWDF      <-0.0711 #%AB #Currently the same as for weaners
# BetaWDFMin   <-0.0345 #%AB
# BetaWDFMax   <-0.1425 #%AB
# BetaWDO      <-0.0711 #%AB #Within pen for others (gilts and pregnant sows sharing pen in the gestation stable). Currently the same as for finishers
# BetaWDOMin   <-0.0345 #%AB
# BetaWDOMax   <-0.1425 #%AB
# BetaBPDW     <-0.01379 #%AB #The value from Broens for post-weaning, ab:+ and pIP=0
# BetaBPDWMin  <-0.01034 #%AB
# BetaBPDWMax  <-0.01782 #%AB
# BetaBPDF     <-0.01379 #%AB #Same as for weaners
# BetaBPDFMin  <-0.01034 #%AB
# BetaBPDFMax  <-0.01782 #%AB
# BetaBPDO     <-0.01379 #%AB #This value needs to be somewhere between BPDF and BPDW (current the same value).  
# BetaBPDOMin  <-0.01034 #%AB
# BetaBPDOMax  <-0.01782 #%AB
# BetaBPDX     <-0.03506 #%AB #Pre-weaning values
# BetaBPDXMin  <-0.01954 #%AB
# BetaBPDXMax  <-0.06322 #%AB
# BetaSOD      <-0.18161 #Currently the same as between piglets
# BetaSOMin    <-0.06552 #Currently the same as between piglets
# BetaSOMax    <-0.50690 #Currently the same as between piglets
# BetaPID      <-0.18161 #Broens value for pre-weaning pigs, ab:+, pIP=1 0.18161  0.06552	0.50690
# BetaPIMin    <-0.06552 
# BetaPIMax    <-0.50690 

#Probabilities for tranmission during day 1 in the life of a piglet
  ProbPND      <-0.76 #Perinatal + all other transmission at day zero. Probability of the offspring being pos given that the dam is pos, not transmission rate. Based on Verhegghe et al., 2013. (Broens: 0.81 or 0.78)
  ProbPNMin    <-0.56 #Based on the error bars in Verhegghe et al., 2013 (which I assume illustrate the 95% CI's even if it is not stated)
  ProbPNMax    <-0.91 #Based on the error bars in Verhegghe et al., 2013 (which I assume illustrate the 95% CI's even if it is not stated)
  
  ProbNPND     <-0.35 #Perinatal + all other transmission at day zero. Probability of the offspring being pos given that the dam is neg., not transmission rate. Based on Verhegghe et al., 2013. (Broens: 0.81 or 0.78)
  ProbNPNMin   <-0.26 #Based on the error bars in Verhegghe et al., 2013 (which I assume illustrate the 95% CI's even if it is not stated)
  ProbNPNMax   <-0.46 #Based on the error bars in Verhegghe et al., 2013 (which I assume illustrate the 95% CI's even if it is not stated)

#Distributions for transmission parameters
  BetaWPW <- rpert(1,BetaWDWMin, BetaWDW,BetaWDWMax) #Within pen for weaners, mean value used as most likely value, min and max: guestimate based on R0 95% CIs
  BetaWPF <- rpert(1,BetaWDFMin,BetaWDF,BetaWDFMax)  #Within pen for finishers, mean value used as most likely value, min and max: guestimate based on R0 95% CIs
  BetaWPO <- rpert(1,BetaWDOMin,BetaWDO,BetaWDOMax) #Within pen for other post-weaning pigs housed in groups (gilts and pregnant sows)
  BetaBPW <- rpert(1,BetaBPDWMin, BetaBPDW, BetaBPDWMax) #Between pens for weaners
  BetaBPF <- rpert(1,BetaBPDFMin,BetaBPDFMin,BetaBPDFMax) #Between pens for finishers
  BetaBPO <- rpert(1,BetaBPDOMin,BetaBPDO, BetaBPDOMax) #Between pens for other post-weaning pigs housed in groups (gilts and pregnant sows)
  BetaBPX <- rpert(1,BetaBPDXMin,BetaBPDX, BetaBPDXMax) #Between pens for pens housing mixed ageagroups (sow+ piglets) in the farrow stable
  BetaPI  <- rpert(1,BetaPIMin, BetaPID,BetaPIMax) #Between piglets housed within the same pen
  BetaSO  <- rpert(1,BetaSOMin, BetaSOD,BetaSOMax) #Between a sow and the piglets suckling at her (no difference between own offspring or piglets she might be functioning as nursery sow for) #All spread at day zero including perinatal transmission included in PND
  BetaBSE <- secfrac*BetaBPX 
  BetaBSEh <- secfrach*BetaBPX #Between sections in work-intensive units (farrowing and mating)
  BetaBSE0 <- 0 #Used for the gestation stable, where there is no sectioning
  BetaBSTA <-stafrac*BetaBPX 
  PrevCuttOff  <- rpert(1,0.5,0.7,1.0)
  ProbPershigh <- rpert(1,0.5,0.75,1)
  ProbPersLow  <- rpert(1,0.01,0.10,0.40)

#Other paramters
  prob.pers.mean <- 0.24 #Mean probability of having the potential to become a persistent shedder #Partly based on Espinosa-Gongora et al., 2015
  prob.pers.sd <- 0.01 #Standard deviation for the above - a guess
  prob.pers <- rnorm(1,prob.pers.mean, prob.pers.sd)
  prob.pers[prob.pers<0] <- 0 #Truncation of the distribution (negative values set to zero)
  DurCarMed <- 7.5 #Median duration of MRSA carriage= 7.5 # / Mean = 10.3 days / SD= 7.7  based on Broens et al., 2012 (no distinction between persistent and intermittent carriers, defined as positive if testet positive once (17.4 if testing positive at two consequtive samplings are required))
  DurCarMin <- 1 #Min duration of MRSA carriage. Broens et al 2012 (10.3 scenario)
  DurCarMax <- 26 #Max duration of MRSA carriage. Broens et al 2012 (10.3 scenario)
  NewInf <- numeric(0)
  InfPens <- numeric(0)
  TBes$Infected <- FALSE
  Pers.carriage <- sample(2:4,length(TBes$Animal), rep=TRUE, prob=c(1-prob.pers,prob.pers*0.99, prob.pers*0.01))
  TBes$Status <-  Pers.carriage #Pigs are assigned the potential to either become intermittent shedders upon exposure (Status=2),
                                #or to become persistent shedders without (Satus=3) or with (Satus=4) probability of recovery
  TBes$Sec <- 0
  TBes$BetaLUF <- 0
  TBes$BetaWPT  <- 0
  TBes$BetaBPT  <- 0
  TBes$BetaBST  <- 0
  TBes$BetaBStT <- 0
  TBes$NumWPInf   <- 0
  TBes$NumBPInf   <- 0
  TBes$NumBSInf   <- 0
  TBes$NumWPAll   <- 0
  TBes$NumBPAll   <- 0
  TBes$NumBSAll   <- 0
  TBes$ShedStatus <- 0 #Current shedder status for the pig
  NewInfected <- NULL

  #Repetions from the generate herd function, that also need to be run every day
  #Define housing unit (stable)  based on age on a given day
  Tmate <-c(240:244) #First stay in the mating unit - time per cycle = 5 days, starting from 8 months of age. Sows/gilts are moved to the gestation unit 4 weeks after mating
  Scycmat<-rep(0:7, each=length(Tmate))#Sow cycle number, repeats depend on length of stay in the relevant unit
  Sclength<-147 #Length of a full sow cycle (No. of dsyd between the sow is returning to a given stable unit)
  Mate<-Tmate+Scycmat*Sclength #Vector including the days in a sows' life spend in the mating unit (based on the sows (adjusted) age in days)

  Tgest <- c(245:353) #First stay in the gestation unit - it is assumed that the duration of the gestation period is 114 days, and that the sow is moved to the farrowing stable, approximately 5 days before farrowing. Time per cycle in thsi unit = 109 days
  Scycges<-rep(0:7, each=length(Tgest)) #Sow cycle number, repeats depend on length of stay in the relevant unit
  Gest<-Tgest+Scycges*Sclength #Vector including the days in a sows' life spend in the gestation unit (based on the sows (adjusted) age in days)

  TFarr <-c(354:386) #First stay in the farrowing unit - time per cycle = 33 days
  Scycfar<-rep(0:7, each=length(TFarr)) #Sow cycle number, repeats depend on length of stay in the relevant unit
  Farr<-TFarr+Scycfar*Sclength #Vector including the days in a sows' life spend in the farrowing unit (based on the sows (adjusted) age in days)

  #Extension of the stay in the mating unit to include the first four weeks of the gestation period
  Lit_no <- c(1:8) #for a given sow (no to be confused with the global litter no: "litter") #This variable is also used in rel. to culling
  Farrowdays <- (Lit_no-1)*Sclength+358 #Day of farrowing
  Matdays <- Farrowdays - 114 #Days of insemination
  Matxtra <- rep(Matdays, each=28) + 1:28
  Mate <- c(Mate, Matxtra)
  Gest <- Gest[!(Gest%in%Mate)]

  #Define and assign slaughterage to the slaughterpigs
  sage <- c(rep(168, times=40),rep(c(161, 175), times=20), rep(c(154, 182), times=10))


### Start of daily loop

  #DaysToSelectInf <- 1:DayOption

 while(gTime<=Days & length(TBes$Animal[TBes$Infected==0])>1){
      gTime <- gTime + 1
      TBes$Age <- TBes$Age + 1 #Cronological age of the pig
      TBes$PACount <- TBes$PACount + 1 #Activity-adjusted age counter
            
      if(sum(is.na(TBes))>0) warning('NAs have been introduced')
     
     
### Farrowing
Lit_no <- c(1:8) #for a given sow (no to be confused with the global litter no: "litter") #This variable is also used in rel. to culling

Farrowdays <- (Lit_no-1)*Sclength+358 #We assumed fixed farrowing days (but the will be some variation due to adjustment of the agecounter)

#Generation of ne litters (For newborns batchnumber is determined by litter)
  IndexSow <- TBes$PACount%in%Farrowdays & TBes$Type==1
  if(sum(IndexSow)>0){
  TBes$LCount[match(TBes$Animal[IndexSow],TBes$Animal)] <- TBes$LCount[match(TBes$Animal[IndexSow],TBes$Animal)] + 1
  for(i in TBes$Animal[IndexSow]){
    tmp  <- as.numeric(max(TBes$Animal)+(1:round(rnorm(n=1, mean=15.9,sd=1.5)))) #Baseret on no. of weaned piglets per litter = 13.8, lifeborn per litter = 15,9
    tmp2 <- as.data.frame(cbind(tmp,FALSE,1,2,0,(max(TBes$teamo[TBes$Type==2])+1), sample(sage, replace=TRUE, length(tmp)),3,(max(TBes$litter)+1),i,0,0,0,0,i,0,0,0,(max(TBes$teama[TBes$Type==2])+1),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
    names(tmp2) <- names(TBes) #Piglets are temporaily all placed in pen 0"
    TBes <- rbind(TBes,tmp2) 
   }
  }

###Cross-fostering / creation of nursery sows to ensure that each sow has a max. of 14 piglets to take care of
# It is assumed that two-step nursery sows are used
tmp7 <- NULL
tmp7 <- subset(TBes, TBes$Type==2&TBes$Age==1) #Get datalines for all day-old piglets on the farm
tmp7 <- tmp7[order(tmp7$litter),] #Sort them
tmp7$ind <- ave(tmp7$Animal, tmp7$litter,  FUN = seq_along) 

ForDistAge1 <- subset(tmp7, tmp7$ind>14&tmp7$litter>0&tmp7$Age==1) #Day-old piglets for distribution among nursery sows
fsti <- 1:45

if(dim(ForDistAge1)[1]>0){
   
  ### Sows that have 8 days old piglets and are selected as nursery sows are getting these replaced by day-old piglets.
   SowsToGetPigletsAge1 <- TBes$Animal[TBes$PACount%in%(Farrowdays+8)&TBes$Type==1&!TBes$ns%in%c(1,2)] 
   NoPigletsFroDistAge1  <- tmp7$Animal[tmp7$ind>14]
   SowsGotAge1         <- SowsToGetPigletsAge1[1:(ceiling(length(NoPigletsFroDistAge1)/14))]
   AnimalSeq            <- rep(1:ceiling(length(NoPigletsFroDistAge1)/14),each=14) #It is assumed taht each sow can maximum handle 14 piglets
   AnimalSeq            <- AnimalSeq[1:length(NoPigletsFroDistAge1)]
   Tmp                  <- sapply(unique(AnimalSeq),function(x)sum(AnimalSeq==x))
   SowPig               <- cbind(NoPigletsFroDistAge1 ,rep(SowsGotAge1,times=Tmp))

   Newsec <- unique(TBes$sec[TBes$Animal %in% NoPigletsFroDistAge1]) #Identify section in the farrowing tsbale the nursery sows will be housed in
   EmptyPens <- fsti[!fsti%in%unique(TBes$sti[TBes$sec==Newsec & TBes$stable==3])] #Identify empty pens

## Removal in case of "pig overflow" - will introduce a warning
   if(length(EmptyPens)<length(SowsGotAge1 )){
     SowsGotAge12 <- SowsGotAge1  
     SowsGotAge1   <- SowsGotAge1[1:length(EmptyPens)]
     SowsGotAge12 <- SowsGotAge12[(length(SowsGotAge1 )+1):length(SowsGotAge12)]
     ToRemove       <- SowPig[SowPig[,2]%in%SowsGotAge12,1]
     SowPig         <- SowPig[SowPig[,2]%in%SowsGotAge1,]
     ToRemAnimals   <- match(ToRemove,TBes$Animal)
     TBes           <- TBes[-ToRemAnimals,]
     #warning("piglets were removed")
    }

    
   TBes$sec[match(SowsGotAge1 ,TBes$Animal)] <- Newsec #Put the nursery sows intro the relevant section
   ToTakePens <- EmptyPens[1:length(SowsGotAge1 )] #Select empty pens and put the new nursery sows into them
   TBes$sti[match(SowsGotAge1 ,TBes$Animal)] <- ToTakePens 
   TBes$ns[match(SowsGotAge1 ,TBes$Animal)]  <- 1 #"Nursery sow indicator": If=1, the pig is currently used as nursery sow
   TBes$nsd[match(SowsGotAge1 ,TBes$Animal)] <- gTime + 28 #"Nursery sows duration": Untill which gTime step will the sow remain a nursery sow
   TBes$PACount[match(SowsGotAge1 ,TBes$Animal)] <- TBes$PACount[match(SowsGotAge1 ,TBes$Animal)] - 7 #Down-adjust the activity-based agecounter for the nursery sows by one week (since they now got one week younger piglets compared to before)
     
   Tmp2            <- match(SowPig[,1],TBes$Animal)
   TBes$nsow[Tmp2] <- SowPig[,2] #This variables indicates who have been nursery sows for a given pig (can be rqual to TBes$dam = the biological mother animal)
   
   for(i in unique(SowPig[,2])){
      Index <- SowPig[,2]==i
      Pigs <- match(SowPig[Index,1],TBes$Animal)
      TBes$sti[Pigs] <-  TBes$sti[match(i,TBes$Animal)]
   }

   #Same procedure for the 8 day old piglets, that now will have to be taken over by other sows
   ForDistAge8 <- NULL
   ForDistAge8 <- TBes[TBes$dam%in%SowsGotAge1  & TBes$dam==TBes$nsow & TBes$Age%in%c(7,8),c(1,15)] #Changed from 7 to 8
  
   if(dim(ForDistAge8)[1]>0){  
   SowsToGetPigletsAge8  <- TBes$Animal[TBes$PACount%in%(Farrowdays+29) & TBes$Type==1 & !TBes$ns%in%c(1,2)]

   TMP3 <- rep(1:length(unique(ForDistAge8$nsow)),times=sapply(unique(ForDistAge8$nsow),function(x)sum(ForDistAge8$nsow==x)))
   ForDistAge8 <- cbind(ForDistAge8,SowsToGetPigletsAge8[TMP3])
   
   TBes$nsow[match(ForDistAge8[,1],TBes$Animal)] <- ForDistAge8[,3]
   
   TBes$sec[match(ForDistAge8[,3],TBes$Animal)] <- unique(TBes$sec[match(ForDistAge8[,1],TBes$Animal)])
   TBes$sti[match(ForDistAge8[,3],TBes$Animal)] <- TBes$sti[match(ForDistAge8[,1],TBes$Animal)]
   TBes$ns[match(ForDistAge8[,3],TBes$Animal)]  <- 1
   TBes$nsd[match(ForDistAge8[,3],TBes$Animal)] <- gTime + 21
   TBes$PACount[match(ForDistAge8[,3],TBes$Animal)] <- TBes$PACount[match(ForDistAge8[,3],TBes$Animal)] - 21
  } 
   TBes$teama <- ifelse(TBes$ns==1,22, TBes$teamo) #Nursery sows are temporaily given batch no. 22
 }

if (any(TBes$nsd==gTime&TBes$ns==1)) {
  nm <- Farrowdays+29 
  TBes$teamo[TBes$nsd==gTime&TBes$ns==1] <- TBes$teama[TBes$PACount %in% nm & TBes$Type==1 & TBes$ns%in%c(0,9999) & TBes$teama<22] [1]  
  TBes$teama[TBes$nsd==gTime&TBes$ns==1] <- TBes$teamo[TBes$nsd==gTime&TBes$ns==1]
  TBes$ns[TBes$nsd==gTime&TBes$ns==1]   <- 0
}

SowsGotAge1  <- NULL

### Sows needing re-insemination
ProbOm <- 0.12 #Probability of an inseminated sow not farrowing after first insemination ( in principle 0.88 for 1 and 8 litters, 0.90 for the rest)

TBes$ri[TBes$rid==gTime] <- 0 #If -1 this indicates that the sow is reinseminated (used as indicator, when "moving" pigs to a new team)
Matdays <- Farrowdays - 114 #The length of the gestation period might have to be adjusted to 118 days all the way thorughout the model
heatcontroldays <- Matdays + 21 #3 weeks after insemination it is checked whether any of the sows have become in heat again and therefore need re-insemination
pregcontroldays <- Matdays + 28 #4 weeks after insemination it is checked whether the sows are pregnant

probreins <- c(1.00, 0.70, 0.55, 0.50, 0.40, 0.30, 0.20, 0) #Probability of re-insemination #Depends on litter count (prob. for 1-8 litters are from vsp, 1.0 for gilts (0 litters) are an assumption)

Sow_mat <- TBes[TBes$Type==1&TBes$stable==1&TBes$PACount%in%heatcontroldays,1] #Sows in the mating stable ready for insemination #ns=-1 is not excluded, since they are inseminated on the day they are identified

notreins <- NULL

if (length(Sow_mat)>0) {
 if(any(TBes$PACount[TBes$Type==1]%in%heatcontroldays)){ #
   tmp5 <- Sow_mat[runif(length(Sow_mat))>=(1-ProbOm)] #Binomial process coded, as an alternative to rbinom
   tmp6 <- TBes[TBes$Animal %in% tmp5,c(1,12)] #Column 11 contains LCount, which indicates what value to use in the probreins vector
   if (nrow(tmp6)==1) {
     tmp6$oml <- rbinom(1,1,probreins[(tmp6[1,2]+1)])
   } else {
     tmp6$oml <- rbinom(tmp6[,1],1,probreins[(tmp6[,2]+1)]) 
   }
  
   if(length(tmp6[tmp6$oml==1,1])>0){
     TBes$ri[TBes$Animal%in%tmp6[tmp6$oml==1,1]]  <- -1 #Just put in to indicate newcomers to a team
     TBes$rid[TBes$Animal%in%tmp6[tmp6$oml==1,1]] <- 1 + gTime  #Just to keep the animals tagged for one day untill moved to the right section and pen
     TBes$PACount[TBes$Animal%in%tmp6[tmp6$oml==1,1]] <- TBes$PACount[TBes$Animal%in%tmp6[tmp6$oml==1,1]] - 21
     VectorX <- c(19:21,1:18) 
     for(i in 1:length(tmp6[tmp6$oml==1,1])){
       TBes$teamo[TBes$Animal%in%tmp6[tmp6$oml==1,1][i]] <- VectorX[TBes$teamo[TBes$Animal%in%tmp6[tmp6$oml==1,1][i]]]
       TBes$teama[TBes$Animal%in%tmp6[tmp6$oml==1,1][i]] <- TBes$teamo[TBes$Animal%in%tmp6[tmp6$oml==1,1][i]]     }
    }
    TBes$sec[TBes$ri<0] <- 0
    TBes$sti[TBes$ri<0] <- 0
   }
   notreins <- tmp5[-which(tmp5%in%tmp6[tmp6$oml==1,1])]
 }


### Movement between stables
TBes$stable <- ifelse(TBes$teama==22&TBes$Type==1, 3, 
                  ifelse(TBes$ri==-1, 1, 
                      ifelse(TBes$Type==2&TBes$PACount>=80, 5, 
                            ifelse(TBes$Type==2&TBes$PACount<80&TBes$PACount>28, 4, 
                                  ifelse(TBes$Type==2&TBes$PACount<29, 3, #Changed so that day 28 now is the weaning stable
                                           ifelse(TBes$Type==1&TBes$PACount%in%Gest, 2, 
                                                  ifelse((TBes$Type==1&TBes$PACount%in%Farr), 3, 1)))))))

#Categories for Stable 1 = Mati, 2 = Gest, 3 = Farr, 4 = Wean, 5 = Fini.

### "Left-over pigs": Approx. 20% 'slow growing leftover-pigs' are to remain in the weaner stable at the time where movement from weaner to finisher stable should usually have taken place.
#Movement is postponed a week for each time a pig is sampled
#Pigs are ssumed to be moved to the finisher stable when the reach an age of 81 days 
ProbRWF <- 0.20 #Probability of remaining in the weaner stable, when the rest of the team is moved to the finisher stable
if (any(TBes$PACount==80)) {
  WeFi <- TBes$Animal[TBes$PACount==80&TBes$Type==2&TBes$Age<110] #Inserted age maximum, as a pig else in principle could be sampled to remain in the stable again and again
  RWeFi <- sample(WeFi,ProbRWF*length(WeFi), replace=FALSE) #Inserted vhronological age counter maximum, as a pig else in principle could be sampled to remain again and again
  TBes$teamo <- ifelse(TBes$Animal%in%RWeFi, TBes$teamo + 1, TBes$teamo + 0)
  TBes$PACount <- ifelse(TBes$Animal%in%RWeFi, TBes$PACount - 7, TBes$PACount)
  TBes$stable[TBes$Animal%in%RWeFi] <- 4 
  TBes$sec[TBes$Animal%in%RWeFi] <- -1 #Currently, these are put in a separate section of the weaner stable, but then leftover pigs from different teams will be mixed there
  penrest <- rep(1:ceiling((length(RWeFi))/30), each=30)
  TBes$sti[TBes$Animal%in%RWeFi] <- penrest[1:length(RWeFi)] #It is assumed that all pigs leaves after 7 days or is selected again, so no update of empty pens are needed
}

###Pigs are sent for slaughter during 3-5 weeks at one specific weekday. The majority of pigs will go at an estimated age around 168 days (24 weeks)
#Is based on PACount and not Age. Slaugtherage has already been assigned: a) When the herd was generated or b) when the litter was generated
#Leftover pigs are moved to the buffertsbale at the same time

#ugt <- seq(from=3, to=Days, by=7)
ugt <- round(seq(from=3, to=Days, by=3.5)) #Days where pigs are sent for slaughter 

Sla <- NULL
if (gTime%in%ugt&length(TBes$Animal[TBes$PACount >= TBes$Slaughterage & TBes$Type==2])>0) {
  Sla <- TBes[TBes$PACount>=TBes$Slaughterage&TBes$Type==2,1]
  if (length(Sla)>0) {
    TBes <- TBes[-which(TBes$Animal%in%Sla),] 
  }
  
  #Leftover pigs are put together and moved into the buffer stable (which is just a separate section in the finisher stable named section -1 (to avoid that this sec ever gets selected as max)):
  #It needs to be ensured that this is only done for pigs close to slaughterage (chronological age limit put to 158)
  per_sec <- aggregate(data=TBes[TBes$stable==5,], Animal ~ sec, function(x) length(unique(x))) #Gives the no. of animals per section
  colnames(per_sec) <- c("sec", "No_of_animals")
  #halfempt  <- max(ifelse(per_sec$No_of_animals<140, per_sec$sec, 0)) 
  halfempt  <- unique(ifelse(per_sec$No_of_animals<140, per_sec$sec, 0))
  halfempt <- halfempt[halfempt>0]
  #halfempt  <- as.vector(ifelse(per_sec$No_of_animals<32&TBes[TBes$Animal%in%per_sec$Animal,]>TBes$Slaughterage[TBes$Animal%in%per_sec$Animal], per_sec$sec, 0)) #Dvs. kun svin med justeret PACount kan ende i bufferstalden.
  movtobuff  <- TBes$Animal[TBes$stable==5&TBes$sec%in%halfempt&TBes$PACount>158]
  
  if (length(movtobuff)>150) {
    buffage <- TBes[TBes$Animal %in% movtobuff, c(1,20)]
    agegrad <- buffage[order(-buffage$PACount),]
    notoberem <- length(movtobuff)-150
    Sla2 <- agegrad$Animal[1:notoberem]
    TBes <- TBes[-which(TBes$Animal%in%Sla2),] 
  }
  
  
  if (length(movtobuff)>0) {
    TBes$sec[TBes$Animal%in%movtobuff] <- -1 #The buffer section if referred to as sec= -1
    TBes$sti[TBes$Animal%in%movtobuff] <-  0 # In case of movetobuff consisting of pigs from more than one pen, I am not keeping track of which ones are from which pen, but just assume they are sorted after size 
    Bsti <- 1:10
    EmptyBufPens <- Bsti[!Bsti%in%unique(TBes$sti[TBes$stable==5&TBes$sec==-1])] #identify empty pens in the buffer section
    pensneeded <- ceiling(length(movtobuff)/15) #How many pens are needed (with 15 animals per pen)
    if (pensneeded > length(EmptyBufPens)) {
      per_pen <- aggregate(data=TBes[TBes$stable==5&TBes$sec==-1,], Animal ~ sti, function(x) length(unique(x)))
      pens_to_be_emptied <- pensneeded - length(EmptyBufPens)
      less_fill_pens <- per_pen$sti[order(per_pen$Animal)] #less than full pens?
      pens_empt <- less_fill_pens[1:pens_to_be_emptied]
      Sla3 <- TBes$Animal[TBes$stable==5&TBes$sec==-1&TBes$sti%in%pens_empt]
      if (length(Sla3)>0) {
        TBes <- TBes[-which(TBes$Animal%in%Sla3),] 
      }
      EmptyBufPens <- Bsti[!Bsti%in%unique(TBes$sti[TBes$stable==5&TBes$sec==-1])]
    }
    pensassigned <- rep(EmptyBufPens[1:pensneeded], each=15)
    TBes$sti[TBes$Animal%in%movtobuff] <-  pensassigned[1:length(movtobuff)] #The pigs have now been put into the empty pens
  } 
}



### Replacement of removed sows by gilts
#With a goal of 23 farrowings per team 0.22*23 -> approx. 5 first parity sows per team
#This means that (as a start) we will model a weekly purchase of 4 gilts (none of the sows will have reached 8 parities, so in the beginning replacement will be slow)- these gilts are bought approx. 7 weeks before they are expected to become ready for the firt insemination
#It is assumed that the gilts will come in heat and be inseminated for the first time, when they are approx. 240 days old (->age at purchase = approx. 191 dage)

### Purchase of gilts: Gilts are in separate teams untill first insemination. 
#Definition of h moved to the very start of the update herd code (also used for defining placement into pens)
if (gTime%in%ugt & length(TBes$Animal[TBes$Type==3])<50) { #50 gilts is in principle a bit too much, but nescessary due to variation
  tmp3 <- as.numeric(max(TBes$Animal)+(1:(50-length(TBes$Animal[TBes$Type==3]))))
  tmp4 <- as.data.frame(cbind(tmp3,FALSE,1,3,191,(max(TBes$teamo[TBes$Type==3])+1),0,1,0,tmp3,0,0,0,0,tmp3,0,0,0,(max(TBes$teamo[TBes$Type==3])+1),191,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
  names(tmp4) <- names(TBes) 
  TBes <- rbind(TBes,tmp4) 
}

### First insemination of gilts / integration into the sow teams
#Expected age when the gilt come in oestrus for the first time: 238-242 days
#k is no. of sows already available for insemination in the mating section
k <- TBes$Animal[TBes$stable==1&TBes$Type==1&TBes$PACount%in%(Matdays-3)] #Gilts available for insemination #Evaluation day changed to matdays-3 to ensure that this will not be evaluated untill after culling of sows
l <- NULL

if (length(k)>0 & length(k)<23 & any(TBes$PACount%in%(Matdays-3))) { #New gilts are added to a sow team if it consistes of less than 23 sows
  samplingpool <- sort(TBes$Animal[TBes$Type==3&TBes$Age>=238]) #Gilts to select from
  l <- samplingpool[1:(26-length(k))]
  l <- l[!is.na(l)]
  TBes$Type[TBes$Animal%in%l] <- 1 #Update status for selected gilts from gilts (type=3) to sows (type=1)
  TBes$teamo[TBes$Animal%in%l] <- rep(unique(TBes$teamo[TBes$Animal%in%k]), length(l)) #Update their batch no
  TBes$teama[TBes$Animal%in%l] <- TBes$teamo[TBes$Animal%in%l] #Update their adjusted batch no.
  #TBes$PACount[TBes$Animal%in%l] <- rep(min(TBes$PACount[TBes$Animal%in%k]), length(l))
  TBes$PACount[TBes$Animal%in%l] <- rep((Matdays[1]-3), length(l)) #not really realisitc, but solves problems with adjusting to a far to high age, which is unrealistic too
  TBes$ns[TBes$Animal%in%l] <- 1000 #ns=1000 is temporaily used as tag for imported gilts
}

##################################
#Placement into sections and pens

#Housing in the farrowing unit

movedays <- Farrowdays - 4 

if (gTime==1) { #Assignment of places in the farrowing unit from the beginning
  teama <- sort(unique(TBes[TBes$Type==1&TBes$stable==3&TBes$teama<22,19]), decreasing=FALSE) #Type needs to be set to one, as the piglets have diff team numbers#Sorted to ensure that team 22 (nursery sows from diff teams) will always end up together in section 6
  FSec <- 1:5 #Sections available for sows that are not nursery sows
  for (i in 1:length(teama)) {
    TBes$sec[TBes$stable==3&TBes$Type==1&TBes$teama==teama[i]] <- i
    TBes$sti[TBes$stable==3&TBes$Type==1&TBes$teama==teama[i]] <- 1:length(TBes$Animal[TBes$stable==3&TBes$Type==1&TBes$teama==teama[i]])
    }
  TBes$sec[TBes$stable==3&TBes$Type==2] <- TBes$sec[match(TBes$nsow[TBes$stable==3&TBes$Type==2],TBes$Animal)]
  TBes$sti[TBes$stable==3&TBes$Type==2] <- TBes$sti[match(TBes$nsow[TBes$stable==3&TBes$Type==2],TBes$Animal)]
}

#Update of places when new are filled/emptied
if (any(TBes$PACount%in%movedays & TBes$Type==1 | TBes$PACount==0)) { #Husk at team 22 skal kunne opdateres på andre dage end movedays
  TBes$sec[TBes$PACount%in%movedays & TBes$Type==1&TBes$stable==3&TBes$ns %in% c(0,9999)] <- 0 #Nulstilling af sektion og dyr for myligt indflyttede dyr
  TBes$sti[TBes$PACount%in%movedays & TBes$Type==1&TBes$stable==3&TBes$ns %in% c(0,9999)] <- 0
  #FSec <- 1:5 #Sections available for sows that are not nursery sows
  EmptyFSec <- FSec[!FSec%in%unique(TBes$sec[TBes$stable==3&TBes$Type==1])]
  TBes$sec[TBes$PACount%in%movedays & TBes$Type==1&TBes$stable==3&TBes$ns %in% c(0,9999)] <- EmptyFSec[1]
  EmptyFPens <- fsti[!fsti%in%unique(TBes$sti[TBes$stable==3&TBes$Type==1&TBes$sec==EmptyFSec[1]])]
  incomfarpigs <- TBes$Animal[TBes$PACount%in%movedays & TBes$Type==1&TBes$stable==3&TBes$ns %in% c(0,9999)&TBes$sec==EmptyFSec[1]]
  if (length(EmptyFPens)<length(incomfarpigs)) { #If the no. of animals exceeds stable capacity, these are removed for sale or culled - shouldn't be nescessary in the farrowing stable, these sows should be removed before mating
    fasold <- sample(incomfarpigs, (length(incomfarpigs)-length(EmptyFPens)), rep=F)
    if (length(fasold)>0) {
      TBes <- TBes[-which(TBes$Animal%in%fasold),] 
      incomfarpigs <- incomfarpigs[!incomfarpigs %in% fasold] 
    }
  }
  TakenFPens <- 1:length(incomfarpigs)
  TBes$sti[TBes$Animal %in% incomfarpigs]  <- TakenFPens
} 

#.... updates when piglets are born..
if (any(TBes$PACount%in%Farrowdays&TBes$stable==3)) {
   TestPen <- merge(TBes[TBes$stable==3&TBes$PACount==0,c(1:24)],TBes[TBes$stable==3&TBes$PACount%in%(Farrowdays),c(15,25:39)], by=c("nsow"), all.x="TRUE") [, union(names(TBes[TBes$stable==3&TBes$PACount==0,c(1:24)]), names(TBes[TBes$PACount%in%Farrowdays,c(16,25:39)]))] 
   TBes <- rbind(TBes[-which(TBes$Animal%in%TestPen$Animal),], TestPen)
   newborns <- TBes$Animal[TBes$PACount==0&TBes$stable==3&TBes$Type==2] #Riskikerer at flytte daggamle grise, der er flyttet til en ammeso tilbage, hvis ikke der indekseres på dette
   #TBes$sec[match(newborns, TBes$Animal)] <- TBes$sec[match(TBes$nsow[TBes$Animal %in% newborns],TBes$Animal[TBes$ns==0])]
   TBes$sec[match(newborns, TBes$Animal)] <- TBes$sec[match(TBes$nsow[TBes$Animal %in% newborns],TBes$Animal)]
   TBes$sti[match(newborns, TBes$Animal)] <- TBes$sti[match(TBes$nsow[TBes$Animal %in% newborns],TBes$Animal)]
  }

#########################################################
#Housing in the gestation unit

gmovedays <- Matdays + 29 #Age where the gestation stable can be enteres
gmoveoutdays <- Farrowdays - 5
gsti <- c(1:12)

if (gTime==1&any(TBes$stable==2)) { #Assignment of places in the gestation unit from the beginning
  teama <- unique(TBes$teama[TBes$stable==2&TBes$Type==1])
  gsti <- seq_along(teama)  
  TBes$sti[TBes$stable==2&TBes$Type==1] <- gsti[match(TBes$teama[TBes$stable==2&TBes$Type==1], teama)] #if want to run with sectioning
  #TBes$sec[TBes$stable==2&TBes$Type==1] <- TBes$sec[TBes$stable==2&TBes$Type==1] #if want to run with sectioning
  TBes$sec[TBes$stable==2&TBes$Type==1] <- 10 #No sectioning in the gestation unit, so all pigs will be assigned section 10
}

if (any(TBes$PACount%in%gmovedays&TBes$Type==1&TBes$stable==2)) { #Need to reset these to zero, when the animals are moved into a new stable or else the pens inhabited in the old stable, will also be taken as occupied in the new stable.
  TBes$sti[TBes$PACount%in%gmovedays&TBes$Type==1&TBes$stable==2] <- 0
  TBes$sec[TBes$PACount%in%gmovedays&TBes$Type==1&TBes$stable==2] <- 0
}

#Update of places when new are filled/emptied
if (any(TBes$PACount%in%gmovedays&TBes$Type==1&TBes$stable==2) | any(TBes$PACount%in%gmoveoutdays&TBes$Type==1&TBes$stable==2)) {
  gteamin <- unique(TBes$teama[TBes$PACount%in%gmovedays&TBes$Type==1&TBes$stable==2])
  EmptyGPens <- gsti[!gsti%in%unique(TBes$sti[TBes$stable==2&TBes$Type==1])]
  TakenGPens <- EmptyGPens[1:length(gteamin)] #Indicate the pen we want to take, not the pens already taken
  TBes$sti[TBes$teama%in%gteamin&TBes$stable==2&TBes$Type==1] <- TakenGPens #It is assumed that in the gestation stable, there is only one pen per section
  #TBes$sec[TBes$teama%in%gteamin & TBes$stable==2&TBes$Type==1] <-  TBes$sti[TBes$teama%in%gteamin & TBes$stable==2] #If want to run with sectioning in the gestation stable 
  TBes$sec[TBes$teama%in%gteamin & TBes$stable==2&TBes$Type==1] <- 10 #No sectioning
}

######################################################################################################
##Housing in the mating unit
mssti <- c(1:40) #40 pens for individual housing of sows #Sows ready for mating + sows that failed to get pregnant, and are waiting for becoming ready for re-insemination again
mpsti <- c(1:12) #12 pens for housing of gilts with room for 5 animals within each pen. Placed in section 5

#Placements of the sows at simulation start
if (gTime==1&any(TBes$stable==1&TBes$Type==1)) {
  teama <- TBes$teama[TBes$stable==1&TBes$Type==1] #When sows to be re-inseminated are not included all sows should belong to the same team
  TBes$sec[TBes$stable==1&TBes$Type==1] <- teama 
  for(i in unique(TBes$sec[TBes$stable==1&TBes$Type==1])){
      Index <- TBes$sec[TBes$stable==1&TBes$Type==1]==i
      TBes$sti[TBes$stable==1&TBes$Type==1][Index]<- 1:sum(Index)
   }
  }  
   
  mmovedays <- c(240,Farrowdays + 29) # Age at arrival in the mating unit

  if (any(TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1)){
    if (any(TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999))) {
      secc <- 1:5 #Sows to be re-inseminated are in a separate section untill this have taken place. Following that they will be moved to the same section as the batch they will now be part of
    # here we reset the pen and section to give them new values
  TBes$sti[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)] <- 0
  TBes$sec[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)] <- 0
  
  #Identification of absolutely empty and partly empty sections and pens, and putting pigs into these
  AbsEmptySec <- secc[!secc%in%unique(TBes$sec[TBes$stable==1&TBes$Type==1])]
    
  if (length(AbsEmptySec)>0) {
    if (length(mssti)<length(TBes$Animal[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)])) { #If the no. of animals exceeds stable capacity, these are removed for sale or culled
      influx <- length(TBes$Animal[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)])
      #if there is to many pigs - some of them will be selected for sale/removal:
      msold <- sample(TBes$Animal[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)], influx-length(mssti), rep=F)
      if (length(msold)>0) {
        TBes <- TBes[-which(TBes$Animal%in%msold),] 
      }
    }
    msteamin <- unique(TBes$teama[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)])
    TBes$sec[TBes$teama%in%msteamin&TBes$stable==1&TBes$Type==1&TBes$ri %in% c(0,9999)] <- rep(AbsEmptySec[1], length(TBes$sti[TBes$teama%in%msteamin&TBes$stable==1&TBes$Type==1&TBes$ri %in% c(0,9999)]))
    EmptyMSPens <- mssti[!mssti%in%unique(TBes$sti[TBes$stable==1&TBes$sec==AbsEmptySec[1]])]
    TakenMSPens <- EmptyMSPens[1:length(TBes$Animal[TBes$PACount%in%mmovedays&TBes$Type==1&TBes$stable==1&TBes$ri %in% c(0,9999)])] #Indicate the pen we want to take, not the pens already taken 
    TBes$sti[TBes$PACount%in%mmovedays&TBes$stable==1&TBes$Type==1&TBes$ri %in% c(0,9999)] <- TakenMSPens  
    } else
      {
    OccupiedMSec <- as.data.frame(table(unlist(TBes$sec[TBes$stable==1&TBes$Type==1]))) 
    colnames(OccupiedMSec) <- c("Sec", "Freq")
    OccupiedMSec$Free <- max(mssti)-OccupiedMSec$Freq
    PartlyEmptyMSec <- OccupiedMSec[OccupiedMSec$Free>0,]  
  }     
 }
}

#Placements of the gilts at start
if (gTime==1&any(TBes$stable==1&TBes$Type==3)) { 
  Gteama <- unique(TBes$teama[TBes$stable==1&TBes$Type==3]) 
  ppens <- as.data.frame(cbind(Gteama,seq_along(Gteama))) 
  colnames(ppens) <- c("teama", "sti")
  ppens$sec <- rep(6, nrow(ppens))
  MaGiSecPen <- merge(TBes[TBes$stable==1&TBes$Type==3,c(1:24,27:39)], ppens, by=c("teama"), all.x="TRUE")[, union(names(TBes[TBes$stable==1&TBes$Type==3,c(1:24,27:39)]), names(ppens))]
  TBes <- rbind(TBes[-which(TBes$Animal%in%MaGiSecPen$Animal),], MaGiSecPen)
}

#Update when new gilts are purchased
if (any(TBes$PACount==191&TBes$Type==3&TBes$stable==1) | gTime%in%(h+5) & any(TBes$PACount==191&TBes$Type==3&TBes$stable==1)) {
  TBes$sec[TBes$PACount==191&TBes$Type==3&TBes$stable==1] <- 6
  TBes$sti[TBes$PACount==191&TBes$Type==3&TBes$stable==1] <- 0
  TotallyEmptyMPPens <- mpsti[!mpsti%in%unique(TBes$sti[TBes$stable==1&TBes$sec==6])]
  #giltteamin <- NULL
  Sharedpenspacesneeded <- 0
  giltteamin <- TBes$Animal[TBes$PACount==191&TBes$Type==3&TBes$stable==1]

  if (length(TotallyEmptyMPPens)>0) {
    #giltteamin <- TBes$Animal[TBes$PACount==192&TBes$Type==3&TBes$stable==1]
    spacesneeded <- length(giltteamin)
    TotEmpSpace <- rep(TotallyEmptyMPPens, each=5)
    TakenMPPens <- TotEmpSpace[1:spacesneeded] #Indicate the pens we want to take, not the pens already taken
    Sharedpenspacesneeded <- sum(is.na(TakenMPPens))
    TotEmpTaken <- length(TakenMPPens) - Sharedpenspacesneeded
    TBes$sti[TBes$Animal%in%giltteamin&TBes$stable==1&TBes$Type==3][1:TotEmpTaken] <- na.omit(TakenMPPens)
  } 
 
  if (length(TotallyEmptyMPPens)==0 | Sharedpenspacesneeded>0) {
    #giltteamin <- TBes$Animal[TBes$PACount==192&TBes$Type==3&TBes$stable==1]
    Occupied <- as.data.frame(table(unlist(TBes$sti[TBes$stable==1&TBes$sec==6])))
    colnames(Occupied) <- c("Pen", "Freq")
    Occupied$Free <- ifelse(Occupied$Freq>4, 0, 5-Occupied$Freq)
    PartlyEmptyMPPens <- Occupied[Occupied$Free>0,]      
    placevector <- c(rep(PartlyEmptyMPPens$Pen, times=PartlyEmptyMPPens$Free))
    Sharedpenspacesneeded <- ifelse(Sharedpenspacesneeded==0, length(giltteamin), Sharedpenspacesneeded)
    Takenplace <- placevector[1:Sharedpenspacesneeded]  
    TBes$sti[TBes$Animal%in%giltteamin&TBes$stable==1&TBes$Type==3&TBes$sti==0] <- Takenplace #TBes$sti==0 to avoid overwriting
  }
}

### Placement, when gilts are inserted as sows (not done on mmovedays to allow for culling to happen before)
if (any(TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000)){
  teamnewins <- unique(TBes$teama[TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000]) 
  nsec <- unique(TBes$sec[TBes$stable==1&TBes$Type==1&TBes$teama==teamnewins&TBes$ri==0&TBes$ns==0])
  TBes$sec[TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000] <- nsec
  EmptyNPens <- mssti[!mssti%in%unique(TBes$sti[TBes$stable==1&TBes$sec==nsec])]
  TakenNPens <- EmptyNPens[1:length(TBes$Animal[TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000])] #Indicate the pens we want to take, not the pens already taken 
  TBes$sti[TBes$Animal %in% TBes$Animal[TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000]] <- TakenNPens
  TBes$ns[TBes$Animal %in% TBes$Animal[TBes$PACount%in%(Matdays-3)&TBes$stable==1&TBes$Type==1&TBes$ns==1000]] <- 0
}



#### Alternative strategy for sows to be re-inseminated: identification 3 weeks after insemination and immediate re-semination and movement to the section of the other newly inseminated sows
#As the sows remains with their old team for the first three weeks, only one movement needs to be coded
omlmat <-TBes$Animal[TBes$ri==-1] #Only tagged for one day now, so can be used as sole condition

if (length(omlmat)>0) {
  nhsec <- unique(TBes$sec[TBes$PACount%in%(Matdays)&TBes$stable==1&TBes$Type==1&TBes$sec>0])
  TBes$sec[TBes$Animal %in% omlmat] <- nhsec
  EmptyolPens <- mssti[!mssti%in%unique(TBes$sti[TBes$stable==1&TBes$sec==nhsec])]
  if (length(omlmat)>length(EmptyolPens)) {
    m2sold <- sample(omlmat, (length(omlmat)-length(EmptyolPens)), rep=F)
    TBes$sec[TBes$Animal %in% m2sold] <- 7
    TBes$sti[TBes$Animal %in% m2sold] <- 1
    IndexOmlMat                       <- which(omlmat%in%m2sold)
    omlmat                            <- omlmat[-IndexOmlMat]
  }
   TakenolPens <- EmptyolPens[1:length(omlmat)] #Indicate the pens we want to take, not the pens already taken 
   TBes$sti[TBes$Animal %in% omlmat] <- TakenolPens
}

###############################################################################################
##Housing in the weaner unit (stable=4)
#It is assumed that pigs are sorted according to size, when assigning pen mates (=random mixing of litters)

rwsti <- rep(1:14, each=30) #30 pigs per pen in 14 pens per section
wdaysp1 <- 29 #Age when entering this unit

#Housing of weaners at the very start: 
if (gTime==1&any(TBes$stable==4)) {
  TBes$sti[TBes$stable==4] <- 0
  TBes$sec[TBes$stable==4] <- 0
  penvec <- NULL
  weagevec <- seq(from=30, to=79, by=7)
  TBes$sec[TBes$stable==4] <- match(TBes$PACount[TBes$stable==4], weagevec)
  for (i in 1:8) {
    nstart <- length(TBes[TBes$stable==4&TBes$sec==i,1])
    penvec <- rep(1:ceiling(nstart/30), each=30)
    TBes$sti[TBes$stable==4&TBes$sec==i] <- penvec[1:nstart] #assigmment of pen
  }
}

#Due to the way the herd was initiated the 8 sections might not nescessarily fill up in chronological order
if (any(TBes$PACount%in%wdaysp1&TBes$stable==4)&gTime>1) {                           
  nw <- length(TBes[TBes$PACount%in%wdaysp1&TBes$stable==4,1])      #No. weaned on the day in question
  if (length(rwsti)<nw) {     #If the no. of animals exceeds the units capacity, surplus pigs are removed for sale or culled
    wsold <- sample(TBes$Animal[TBes$PACount%in%wdaysp1&TBes$stable==4], nw-length(rwsti), rep=F)
    if (length(wsold)>0) {
      TBes <- TBes[-which(TBes$Animal%in%wsold),] 
      nw <- length(TBes[TBes$PACount%in%wdaysp1&TBes$stable==4,1])
    }
  }
  sti <- sample(rwsti, nw, replace=FALSE) 
  
  WeaSec <- 1:8
  
  TBes$sti[TBes$PACount%in%wdaysp1&TBes$stable==4] <- 0
  TBes$sec[TBes$PACount%in%wdaysp1&TBes$stable==4] <- 0
  
  EmptyWSec <- WeaSec[!WeaSec%in%unique(TBes$sec[TBes$stable==4])][1]
  
  TBes$sti[TBes$PACount%in%wdaysp1&TBes$stable==4] <- rwsti[1:length(TBes$Animal[TBes$PACount%in%wdaysp1&TBes$stable==4])]
  TBes$sec[TBes$PACount%in%wdaysp1&TBes$stable==4] <- EmptyWSec

}

###############################################################################################
##Housing in the finisher unit (stable=5)

fimoveday <- 80 #PACount, when entering this unit

frsti <- rep(1:24, each=15) #24 pens per section with room for 15 animals in each
redsti <- rep(1:3, each=15) #3 pens in the buffer section 

#Housing of finishers at the very start: 
if (gTime==1&any(TBes$stable==5)) {
  fiagevec <- seq(from=86, to=177, by=7)
  TBes$sec[TBes$stable==5] <- match(TBes$Age[TBes$stable==5], fiagevec)
  for (i in 1:12) {
    nfstart <- length(TBes[TBes$stable==5&TBes$sec==i,1])
    penfvec <- rep(1:ceiling(nfstart/15), each=15)
    TBes$sti[TBes$stable==5&TBes$sec==i] <- penfvec[1:nfstart]
  }
}

#update when animals are moved
if (any(TBes$PACount%in%fimoveday&TBes$stable==5)) { 
  
  TBes$sti[TBes$PACount%in%fimoveday&TBes$stable==5] <- 0
  TBes$sec[TBes$PACount%in%fimoveday&TBes$stable==5] <- 0
    
  nf <- length(TBes[TBes$PACount%in%fimoveday&TBes$stable==5,1]) #No. of weaners moved into the slaughterpig stable on the day in question
  if (length(frsti)<nf) { #If the no. of animals exceeds stable capacity, these are removed for sale or culled
    fisold <- sample(TBes$Animal[TBes$PACount%in%fimoveday&TBes$stable==5], nf-length(frsti), rep=F)
    if (length(fisold)>0) {
      TBes <- TBes[-which(TBes$Animal%in%fisold),] 
    }
    nf <- length(TBes[TBes$PACount%in%fimoveday&TBes$stable==5,1])
  }
  if (nf>150) {
    sti <- sample(frsti, nf, replace=FALSE)
    } else {
    fipenvec <- rep(1:ceiling(nf/15), each=15)
    sti <- fipenvec[1:nf]
    }
  
  FinSec <- 1:14 
  EmptyFinSec <- FinSec[!FinSec%in%unique(TBes$sec[TBes$stable==5])][1]

  TBes$sec[TBes$PACount%in%fimoveday&TBes$stable==5] <- EmptyFinSec
  TBes$sti[TBes$PACount%in%fimoveday&TBes$stable==5&TBes$sec==EmptyFinSec] <- frsti[1:length(TBes$Animal[TBes$PACount%in%fimoveday&TBes$stable==5&TBes$sec==EmptyFinSec])]
  
}

#### Snout contact between pens (two and two)
#Note: curently not used in the model
#Need to be "cut-off" if exceeding the max. no. of pens within the section and/or all pens don't have a partner pen (or else this will just have to be taken into account in the colonization model)
#Definition: Pigs in pens with an odd pen no. can have snout contact to pigs in the pen having their pen no. + 1 / if even: own pen no. -1

is.even <- function(x) x %% 2 == 0
z <- TBes$sti[is.even(TBes$sti)=="TRUE"]
TBes$snoutcont <- ifelse(TBes$sti%in%z,TBes$sti-1, TBes$sti+1)

if (any(TBes$sti==23&TBes$stable==5) | any(TBes$sti==15&TBes$Type==3)) {
  TBes$snoutcont[TBes$sti==21&TBes$stable==5 | TBes$sti==15&TBes$Type==3] <- 0 
}

#####################################
### Deaths/culling of sows
### Probability vectors
Sow_raw_mat_far<-c(0,0.020, 0.016, 0.017, 0.027, 0.018, 0.023, 0.020, 0.034,0.034) # Sow fatality rate between mating and farrowing, based on no. og litters, source: "VSP/SEGES (2013): Udsætningstrategi"
Sow_raw_post_far<-c(0,0.028, 0.026, 0.031, 0.034, 0.031, 0.025, 0.042, 0.028,0.028) #As above, after farrowing (source: "VSP/SEGES (2013): Udsætningstrategi")
Sow_raw_cull_post_far<-c(0,0.04, 0.06, 0.08, 0.11, 0.18, 0.34, 0.45, 1.00,1) #Prop of sows culled after the a given litter is weaned #source: "VSP/SEGES (2013): Udsætningstrategi".(Prop. of sows alive after haven given birth to the relevant no. of litters)

## Day_leav_fst<-c(386, 533, 680, 827, 974, 1121, 1268, 1415) #Day of leaving the farrowing stable after having the given litter no.
Sow_mat_far_daily <-Sow_raw_mat_far/114 #Daily probability of removal between mating and farrowing
Sow_post_far_daily <-Sow_raw_post_far/33 #Daily probability of removal post-farrowing (and pre-mating)

#Removal of sows of high age (nescessary because the litter counter is starting from zero at the the start of model run for sows of higher age)
if (any(TBes$PACount > 1415)) {
  sowsout0 <- TBes$Animal[TBes$PACount>1415]
 if (length(sowsout0)>0) {
  TBes <- TBes[-which(TBes$Animal%in%sowsout0),] 
 }
}

##################################
#Removal of sows after weaning

if (length(notreins)>0) {
  TBes <- TBes[-which(TBes$Animal%in%notreins),] 
}

lit_dam  <- matrix(numeric(0),ncol=3) #Nescessary to avoid carry-over
lit_dam2 <- matrix(numeric(0),ncol=3)
lit_dam3 <- matrix(numeric(0),ncol=3)


#"Strategic culling: selection and removal of sows immediately after weaning 
Postwean_newprob <- Sow_raw_cull_post_far - c(0, 0.12*(1-probreins), 1) #probreins is the probability of re-insemination
#Farrowdays + 29 = Matdays - 4 =Flyttedag fra farestalden til løbestalden
if (any(TBes$PACount %in% (Farrowdays + 29)&TBes$Type==1)) {
Movedsows  <- TBes$Animal[TBes$PACount %in% (Farrowdays + 29)&TBes$Type==1]
lit_dam    <- cbind(Movedsows,TBes$LCount[match(Movedsows,TBes$Animal)])
lit_dam    <- cbind(lit_dam,rbinom(lit_dam[,1], 1, Postwean_newprob[(lit_dam[,2]+1)])) 
if (sum(lit_dam[,3])>0) {
 TBes <- TBes[-which(TBes$Animal%in%lit_dam[lit_dam[,3]==1,1]),] 
}
}

######################################
#Removal of sows during any other time
#Between insemination and farrowing
gestdur <- 147 -28 -5 -1 #One sow cycle is 147 days, minus the 28 days spent on lactation and the 5 days as "dry sow"
Matfardays <- rep(Matdays, each=gestdur) + c(1:gestdur)
matfarsow <- TBes$Animal[TBes$PACount %in% Matfardays & TBes$Type==1]
if (any(TBes$PACount %in% Matfardays & TBes$Type==1)) {
  lit_dam2    <- cbind(matfarsow,TBes$LCount[match(matfarsow,TBes$Animal)])
  lit_dam2    <- cbind(lit_dam2,rbinom(lit_dam2[,1], 1, Sow_mat_far_daily[(lit_dam2[,2]+1)])) 
  if (sum(lit_dam2[,3])>0) {
    TBes <- TBes[-which(TBes$Animal%in%lit_dam2[lit_dam2[,3]==1,1]),] 
 }
}

######################################
#Between farrowing and weaning
wdays    <- rep(Farrowdays, each=28) + 1:28
lactsows <- TBes$Animal[TBes$PACount %in% wdays & TBes$Type==1 &TBes$ns %in% c(0,9999)] #This is a simplification, beacause something can happen to sows with ns %in% c(1,2) as well.
if (any(TBes$PACount %in% wdays & TBes$Type==1)) {
  lit_dam3 <- cbind(lactsows,TBes$LCount[match(lactsows,TBes$Animal)])
  lit_dam3    <- cbind(lit_dam3,rbinom(lit_dam3[,1], 1, Sow_post_far_daily[(lit_dam3[,2]+1)]))
  
  if (length(lit_dam3[lit_dam3[,3]==1,1])>0) {
   piglremov <- TBes$Animal[TBes$Age%in%1:28 & TBes$nsow%in%lit_dam3[lit_dam3[,3]==1,1]]
   pigsowrem <- c(lit_dam3[lit_dam3[,3]==1,1], piglremov)
   TBes <- TBes[-which(TBes$Animal%in%pigsowrem),] 
   }
}

#########################################
#Collecting all removal vectors for sows for the output file
removedsows <-0
removedsows <- sum(lit_dam[,3]==1) + sum(lit_dam2[,3]==1) + sum(lit_dam3[,3]==1)

##########################################
### Mortality for weaners and finishers
#Assumptions re. mortality (source:"Landsgennemsnit i svineproduktionen" (2015 figures)): 3.7% for finishers herd incl. pigs not approved at slaughterhouse and 3.1% for weaner herds

ptb3       <- TBes$Animal[TBes$Type==2&TBes$Age>28]
ptb3Prob   <- ifelse(TBes$PACount[match(ptb3,TBes$Animal)]<80, 0.031/(80-28), (0.037-0.031)/(165-80))
Removeptb3 <- ptb3[rbinom(ptb3, 1, ptb3Prob)==1]
if (length(Removeptb3)>0) {
  TBes <- TBes[-which(TBes$Animal%in%Removeptb3),] 
}

### Mortality for piglets
# Overall mortality per litter for piglets between birth and expected time of weaning at day 28 have been set to: 13.4% (Landsgennemsnit)
# Piglets to be removed are selected independent of which litter they belong to, but maybe they should be selected litterwise?
probPM <- 0.134/28 #mean prob per day of piglets getting culled
daydistr <- c(0.28, 0.15, 0.11, rep(0.05, times=4), rep(0.024, times=7), rep(0.007, times=15)) #Agedependent probability vector #Numbers estimated from figure 5, p. 17 in DJF report on piglet mortality in DK
probPD <- daydistr*0.134 #differentiated prob per day of piglets getting culled

#Using differentiated prob
piglets    <- cbind(TBes$Animal[TBes$stable==3&TBes$Type==2],(TBes$Age[TBes$stable==3&TBes$Type==2]+1))
piglets    <- cbind(piglets,rbinom(piglets[,1], 1, probPD[(piglets[,2])])) 
remove3 <- piglets[piglets[,3]==1,1]


#No. of animals in the pen, the piglet will be removed from:
TBes$Sti <- TBes$stable*10000 + TBes$sec*100 + TBes$sti*1 #create s pen-ID unique across units and sections
Prf <- TBes$Sti[match(remove3, TBes$Animal)] 
AniRem <- NULL
AniRem <- data.frame(table(TBes$Sti[TBes$Sti %in% Prf]))

#AniRem[6,2] <-2 

if (any(AniRem$Freq<3)) {
  sowstobemoved <- TBes$Animal[TBes$Type==1&TBes$Sti%in%AniRem$Var1[AniRem$Freq<3]]
  #tmpvar <- rep(mmovedays, each=5) + c(0:4) #Only days prior to insemination and the insemination day have been included
  secmin <- min(unique(TBes$sec[TBes$stable==1]))
  distm  <- min(TBes$PACount[TBes$stable==1&TBes$Type==1&TBes$sec==secmin])-mmovedays #Når alderen på dyrene i en sektion kendes, kan alderen på dem i de øvrige sektioner forudsiges
  remainder <- unique(distm%%7) #Days the youngest animals are older than mmovedays
  adj <- ifelse(remainder>5,(7-remainder)+7 ,remainder) 
  adj2 <- as.matrix(sapply(sowstobemoved,function(x)mmovedays-TBes$PACount[TBes$Animal%in%x]),drop=F)
  adj3 <- apply(adj2,2,function(x)min(x[x>0])) 
  adj4 <- adj3 + adj
  moveday <- ifelse(remainder>5,gTime+(7-remainder) ,gTime)
  New.team <- unique(TBes$teama[TBes$PACount %in% (mmovedays+remainder)])
  
  #Removal of sows that cannot be moved to the mating section immediatly 
  ToBeRemoved <- sowstobemoved[adj3%%7!=0]
  if(length(ToBeRemoved)>0) TBes <- TBes[-which(TBes$Animal%in%ToBeRemoved),]
    
 if (gTime==moveday) {
    TBes$stable[TBes$Animal%in%sowstobemoved] <- 1
    TBes$sec[TBes$Animal%in%sowstobemoved] <- 0
    TBes$sti[TBes$Animal%in%sowstobemoved] <- 0
    TBes$PACount[TBes$Animal%in%sowstobemoved] <- TBes$PACount[TBes$Animal%in%sowstobemoved]+adj3
    TBes$teama[TBes$Animal%in%sowstobemoved] <- New.team
    TBes$teamo[TBes$Animal%in%sowstobemoved] <-  TBes$teama[TBes$Animal%in%sowstobemoved]
    TBes$sec[TBes$Animal%in%sowstobemoved] <-max(unique(TBes$sec[TBes$Type==1&TBes$stable==1&TBes$teama==New.team]))
    EmpP <- mssti[!mssti%in%unique(TBes$sti[TBes$stable==1&TBes$sec==TBes$sec[TBes$Animal%in%sowstobemoved]])]
    TBes$sti[TBes$Animal%in%sowstobemoved] <- EmpP[1] #Based on that this is unlikely to happen to more than one sow within a round
    TBes$ns[TBes$Animal%in%sowstobemoved] <- 0
    TBes$nsd[TBes$Animal%in%sowstobemoved] <- 0
  }
}
  

#Removal of the selected piglets
if (length(remove3)>0) {
 TBes <- TBes[-which(TBes$Animal%in%remove3),] 
}


####################################################
# Infection model starts here
####################################################

#Currently a burn-in period of 4 years are applied before introduction (defined by "DiseaseStart"). 
#This is nescessay for the number of animls in the herd to stabilize

if(gTime>=(DiseaseStart)){

if (any(TBes$Status==1)) {
  Pers.carriage <- sample(2:4,length(TBes$Animal[TBes$Status==1]), rep=TRUE, prob=c(1-prob.pers,prob.pers*0.99, prob.pers*0.01))
  TBes$Status[TBes$Status==1] <-  Pers.carriage
}

if(gTime==DiseaseStart) { 
  Ini.inf <- sample(TBes$Animal[TBes$Type==3],1) #select type and no. of animals to be infected initially
  #Ini.inf <- sample(TBes$Animal[TBes$stable==4&TBes$Age%in%28:35],1)
  #Pers.carriage <- rbinom(length(Ini.inf), 1, prob.pers) 
  TBes$Infected[TBes$Animal %in% Ini.inf] <- TRUE 
  TBes$Status[TBes$Animal %in% Ini.inf] <- 2 #Define shedder status for the initially infected: 2=intermittent shedder, 3=persistent
  TBes$dcd[TBes$Animal %in% Ini.inf] <- gTime + round(rpert(1, DurCarMin, DurCarMed, DurCarMax)) #Define decolonization day for the animals
  TBes$nsr[TBes$Animal%in%Ini.inf] <- gTime+1 #Used for tagging newly infected animals on their first day of infection
   }
   
TBes$Sti <- TBes$stable*10000 + TBes$sec*100 + TBes$sti*1 # Need to match sorting after loops later on, and be run after movements on the day have been finished
# A new variable i introduced.... make sure that we introduce it at start
TBes$Sec <- TBes$stable*10000 + TBes$sec*100
#A New variable to look up beta, based on stable and type
TBes$BetaLUF <- TBes$stable*10 + TBes$Type * 1 

IndexXMat <- cbind(c(11,13,21,31,32,42,52),c(1:7))

for(i in 1:(dim(IndexXMat)[1])) TBes$BetaLUF[TBes$BetaLUF==IndexXMat[i,1]] <- IndexXMat[i,2]

# Here we create a look up table to represent the within-pen between-pens, between-sections and between-stable transmission rates based on the types
# First we reset all values to 0 to make sure that introductions and removal of pigs did not affect the infections module
# MAKE sure that we introduce these variables at start

WithinPenBeta   <- c(kal11*BetaWPO,kal13*BetaWPO,kal21*BetaWPO,kal31*BetaSO,kal32*BetaSO,kal42*BetaWPW,kal52*BetaWPF) #The betas for the farrowing stable will be overwritten later
BetweenPenBeta  <- c(kal11*BetaBPO, kal13*BetaBPX, kal21*BetaBPX, kal31*BetaBPX, kal32*BetaBPX, kal42*BetaBPW,kal52*BetaBPF)
BetweenSecBeta  <- c(kal11*BetaBSEh,kal13*BetaBSE,kal21*BetaBSE0,kal31*BetaBSEh,kal32*BetaBSEh,kal42*BetaBSE,kal52*BetaBSE)
BetweenStabBeta <- c(kal11*BetaBSTA,kal13*BetaBSTA,kal21*BetaBSTA,kal31*BetaBSTA,kal32*BetaBSTA,kal42*BetaBSTA,kal52*BetaBSTA)

TBes$BetaWPT  <- WithinPenBeta[TBes$BetaLUF]  
TBes$BetaBPT  <- BetweenPenBeta[TBes$BetaLUF]  
TBes$BetaBST  <- BetweenSecBeta[TBes$BetaLUF]  
TBes$BetaBStT <- BetweenStabBeta[TBes$BetaLUF]

# Number of infected animals
TBes$NumWPInf <- ave(TBes$Infected,TBes$Sti,FUN=sum)
TBes$NumBPInf <- ave(TBes$Infected,TBes$Sec,FUN=sum)
TBes$NumBSInf <- ave(TBes$Infected,TBes$stable,FUN=sum) 

#Total number of animals
TBes$NumWPAll <- ave(TBes$Infected,TBes$Sti,FUN=length) #in the pen
TBes$NumBPAll <- ave(TBes$Infected,TBes$Sec,FUN=length) #in the section
TBes$NumBSAll <- ave(TBes$Infected,TBes$stable,FUN=length) #in the stable


PrevelanceSec <- TBes$NumBPInf/TBes$NumBPAll

TBes$nsr      <- ifelse(PrevelanceSec>PrevCuttOff,ProbPershigh,ProbPersLow) #Probability of becoming a persisent carrier for those pre-disposed

# Here we calculate the probability of infection
ProbInfWP  <- 1-exp(-TBes$BetaWPT*TBes$NumWPInf/TBes$NumWPAll)
ProbInfBP  <- ifelse(TBes$sec==-1&TBes$NumWPAll==TBes$NumBPAll,0, 1-exp(-TBes$BetaBPT*(TBes$NumBPInf-TBes$NumWPInf)/(TBes$NumBPAll-TBes$NumWPAll)))
ProbInfBS  <- ifelse(TBes$stable==2,0,1-exp(-TBes$BetaBST*(TBes$NumBSInf-TBes$NumBPInf)/(TBes$NumBSAll-TBes$NumBPAll)))
ProbInfBSt <- 1-exp(-TBes$BetaBStT*(sum(TBes$Infected)-TBes$NumBSInf)/(length(TBes$Animal)-TBes$NumBSAll))
TotProbInf <- 1-((1-ProbInfWP)*(1-ProbInfBP)*(1-ProbInfBS)*(1-ProbInfBSt))
#TotProbInf[is.na(TotProbInf)]   <- 0


if(sum(TBes$Age==0)>0){
  IndexAge1 <- TBes$Age==0 & TBes$Infected[match(TBes$nsow,TBes$Animal)]==1
  TotProbInf[IndexAge1]   <- rpert(1,ProbPNMin,ProbPND, ProbPNMax)

  IndexAge2 <- TBes$Age==0 & TBes$Infected[match(TBes$nsow,TBes$Animal)]==0 & sum(TBes$Infected[TBes$Sec==unique(TBes$Sec[TBes$Age==0])])>0
  TotProbInf[IndexAge2]   <- rpert(1,ProbNPNMin,ProbNPND, ProbNPNMax)

  IndexAge3 <- TBes$Age==0 & TBes$Infected[match(TBes$nsow,TBes$Animal)]==0 & sum(TBes$Infected[TBes$Sec==unique(TBes$Sec[TBes$Age==0])])==0
  TotProbInf[IndexAge3] <- 0
 }

SUSAnimals  <- TBes$Animal[TBes$Infected==0]
NewInfected <- SUSAnimals[rbinom(length(SUSAnimals),1,TotProbInf[TBes$Animal%in%SUSAnimals])==1]
if(length(NewInfected)) {
   TBes$Infected[TBes$Animal%in%NewInfected] <- 1
   TBes$ShedStatus[TBes$Animal%in%NewInfected & TBes$Status==2]           <- 2
   TBes$ShedStatus[TBes$Animal%in%NewInfected & TBes$Status%in%3:4]       <- ifelse(rbinom(n=sum(TBes$Animal%in%NewInfected & TBes$Status%in%3:4),1,prob=TBes$nsr[TBes$Animal%in%NewInfected & TBes$Status%in%3:4]),3,2)
   TBes$dcd[TBes$Animal%in%NewInfected&TBes$ShedStatus==2]                <- gTime + round(rpert(1, DurCarMin, DurCarMed, DurCarMax)) #Decolonization day for intermittent carriers 
   TBes$dcd[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$Status==4] <- gTime + round(rpert(1, (DurCarMin+100), (DurCarMed+100), (DurCarMax+100))) #Decolonization day for persistent carriers, for those where it might be possible (Status=4) 
  }

#Decolonization 
TBes$Infected[TBes$Status%in%c(2,4) & gTime==TBes$dcd] <- FALSE 

}#END of if(gTime>=DiseaseStart)){


#Creation of various variables to be included in output
weanedpigs <- numeric(0)
weanedlitters <- numeric(0)

if (any(TBes$PACount==28)) {
  weanedpigs <- TBes$Animal[TBes$PACount==28]
  weanedlitters <- unique(TBes$litter[TBes$PACount==28])
}

insemipigs <- numeric(0)
reinsemipigs <- numeric(0)
if (any(TBes$PACount%in%Matdays)) {
  insemipigs <- TBes$Animal[TBes$PACount%in%Matdays&TBes$Type==1]
  reinsemipigs <- TBes$Animal[TBes$ri<0]
}

livebornpigs <- numeric(0)
littersborn  <- numeric(0)
if (any(TBes$PACount==0)) {
  livebornpigs <- TBes$Animal[TBes$PACount==0]
  littersborn  <- unique(TBes$litter[TBes$PACount==0])
}

firstparitysows <- numeric(0)

if (any(TBes$PACount%in%(Farrowdays-1)&TBes$LCount==0)) {
  firstparitysows <- TBes$Animal[TBes$PACount%in%(Farrowdays-1)&TBes$LCount==0] #Needs to be counted the day before farrowing, since one is added to LCount at the farrowday
} #Littersborn can be used to determine how many sows gave birth within a given year

#Create output file
Outputherd <- rbind(Outputherd, c(iteration, gTime, length(TBes$Animal[TBes$Type==1]), length(TBes$Animal[TBes$Type==2]), length(Sla), removedsows, length(l), 
                                  length(TBes$Animal[TBes$Type==3]), length(remove3), length(Removeptb3), length(weanedpigs), length(weanedlitters), length(insemipigs), 
                                  length(reinsemipigs),length(livebornpigs), length(littersborn), length(firstparitysows), sum(TBes$Infected), length(TBes$Animal[TBes$ShedStatus==3]),BetaWPW, 
                                  sum(TBes$Infected[TBes$stable==1]), length(TBes$Animal[TBes$ShedStatus==3&TBes$stable==1]), sum(TBes$Infected[TBes$stable==2]), length(TBes$Animal[TBes$ShedStatus==3&TBes$stable==2]),
                                  sum(TBes$Infected[TBes$stable==3]), length(TBes$Animal[TBes$ShedStatus==3&TBes$stable==3]), sum(TBes$Infected[TBes$stable==4]), length(TBes$Animal[TBes$ShedStatus==3&TBes$stable==4]),
                                  sum(TBes$Infected[TBes$stable==5]), length(TBes$Animal[TBes$ShedStatus==3&TBes$stable==5]), sum(TBes$Infected[TBes$Type==1]), length(TBes$Animal[TBes$ShedStatus==3&TBes$Type==1]),
                                  sum(TBes$Infected[TBes$Type==2]), length(TBes$Animal[TBes$ShedStatus==3&TBes$Type==2]),  sum(TBes$Infected[TBes$Type==3]), length(TBes$Animal[TBes$ShedStatus==3&TBes$Type==3]),
                                  sum(TBes$Infected[TBes$Type==2&TBes$PACount<29]), length(TBes$Animal[TBes$ShedStatus==3&TBes$Type==2&TBes$PACount<29]), length(TBes$Animal[TBes$stable==1&TBes$Type==1]),
                                  length(TBes$Animal[TBes$stable==1&TBes$Type==3]), length(TBes$Animal[TBes$stable==2&TBes$Type==1]), length(TBes$Animal[TBes$stable==3&TBes$Type==1]), length(TBes$Animal[TBes$stable==3&TBes$Type==2]), 
                                  length(TBes$Animal[TBes$stable==4&TBes$Type==2]),length(TBes$Animal[TBes$stable==5&TBes$Type==2]), 
                                  length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$Type==1]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$Type==1]),
                                  length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$Type==2]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$Type==2]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$Type==3]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$Type==3]),  
                                  length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==1]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$Status==3&TBes$stable==1]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==2]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$stable==2]),
                                  length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==3]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$Status==3&TBes$stable==3]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==4]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$stable==4]),
                                  length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==5]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$Status==3&TBes$stable==5]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==2&TBes$stable==3&TBes$Type==2]), length(TBes$Animal[TBes$Animal%in%NewInfected&TBes$ShedStatus==3&TBes$stable==3&TBes$Type==2])
))

if(dim(Outputherd)[1]>10000){
  NAME <- paste(runID,"MRSA.txt",sep="-")
  write.table(Outputherd,NAME,append=T,sep=" ",col.names = F,row.names=F)                                                     
  Outputherd <- matrix(numeric(0),ncol=65)
}

if(gTime==Days){
  NAME <- paste(runID,"MRSA.txt",sep="-")
  write.table(Outputherd,NAME,append=T,sep=" ",col.names = F,row.names=F)                                                      
  Outputherd <- matrix(numeric(0),ncol=65)
}

  }#End of gTime "while" loop

}#End of iteration "for" loop

#}# End of Function


###Temperature dependant Forage-rest state model

#amount of time to run modle in hours
Time<-1000

#Size of prey in l
prey<-.2

#prey calorie density(Calories per l)
pval<-1200

#oxycaloric coefficient (mgO2/Kilocal)
OxyCal<-0.000324

#Foraging Sucess(probability of catching prey when foraging)
pSuc<-.1

#temperature of shark(C)
temp<-25

#Body size(kg)
body<-1

#Volume in stomach that causes an individual to have a 50 percent chance to stop foraging/Become Full
FRs50=15

#log rate constant that changes the shape of the steepness of stopping foraging
FRlrc=1

#Volume in stomach that cause an individual to have a 50 percent chance to start foraging/ become hungry
RFs50=3

#log rate constant that changes the shape of the steepness of changing to foraging
RFlrc=1

#temperature at which 50 percent of max gut evacuation rate
rEvd50=20

#rate constant that regulates the rate of gut evacuation increase
rEvlrc=4

#maximum rate of gut evacuation in liters per hour
rEvmax=.5

#Temperature at mean gut percent gut absorbtion
rAbd50=20

#constant that regulates the rate of gub absorbtion decrease
rAblrc=4

#maximum gut absorbtion efficiciency
rAbmax=.9

#minimum gut absorbtion efficiency
rAbmin=.5

#metabolic rate at 0 in mg 02 per kilo per hour
m_Int=16.8272

#coefficienct in front of squared term effecting rate of metabolic increase
m_Slope=.2315

#difference in temperature between foraging and resting state(negative is forage cool,positive is rest cool)
FRTdif=0



#probability of foraging to resting
pFR<-function(x, s50=.7, lrc=.1){
  1/(1+exp((s50-x)/lrc))
}
#probability of resting to foraging
pRF<-function(x, s50=.1, lrc=.1){
  1 -1/(1+exp((s50-x)/lrc))
}

#Gut Evauation Rate()
rEvac<-function(temp,d50=20,lrc=4,max=.1){
  (max)/(1+exp((d50-temp)/lrc))
}

#Gut Absorbtion efficiency()
rAbsorb<-function(temp,d50=20,lrc=4,max=.9,min=.5){
  max+(min-max)/(1+exp((d50-temp)/lrc))
}


#bodyCore
bodyCore<-function(Tb,Ta,k=0.00515,ksd=0.0008){
  Tb+ rnorm(1,k,ksd) * (Ta-Tb)
}

#metabolism
Rmetabolism<-function(temp,m_Int=16.83,m_Slope=.231528){
  (m_Int + m_Slope*temp^2)
}


#Assumptions and Rules of the Model
plot(pFR(seq(0,3,by=.1))~seq(0,3,by=.1),ylab="Probabily of transition from foraging to resting",xlab="Stomach Fullness",type="l")
plot(pRF(seq(0,3,by=.1))~seq(0,3,by=.1),ylab="Probabily of transition from Resting to foragings",xlab="Stomach Fullness",type="l")
plot(rEvac(seq(0,30,by=.1))~seq(0,30,by=.1),ylab="Gut Evacuation Rate", xlab="Temperature")
plot(rAbsorb(seq(0,30,by=.1))~seq(0,30,by=.1),ylab="Percent of evacuated contents absorbed", xlab="Temperature")
plot(Rmetabolism(seq(0,30,by=.1))~seq(0,30,by=.1),ylab="Metabolic Rate", xlab="Temperature")
plot(rep(prey,301)~seq(0,30,by=.1),type="l",ylab="Foraging amount",xlab="Temperature")
plot(rep(pSuc,301)~seq(0,30,by=.1),type="l",ylab="Foraging Sucess",xlab="Temperature")


TempBudgetModel<-function(prey=.1,
                          pval=1200,
                          pSuc=.1,
                          temp=25,
                          body=1,
                          OxyCal=0.000324,
                          Time=1000,
                          FRs50=.7, FRlrc=.1,
                          RFs50=.1, RFlrc=.1,
                          rEvd50=20,rEvlrc=4,rEvmax=.1,
                          rAbd50=20,rAblrc=4,rAbmax=.9,rAbmin=.5,
                          m_Int=16.8272,m_Slope=.2315,
                          FRTdif=0
){
  #Amount of energy accumulated
  energy<-0
  #amount of energy absorbed
  absorbed<-0
  #amount in stomach at start
  stomach<-0
  #State at start
  state<-1
  temp0<-temp
  
  stateRec<-matrix(0,ncol=5,nrow=Time)
  colnames(stateRec)<-c("State","Temp","Stomach","Energy","FSucess")
  
  for(i in 1:Time){
    evacuated<-rEvac(temp,rEvd50,rEvlrc,rEvmax)
    
    if(evacuated>stomach){
      evacuated<-stomach
    }
    stomach<-stomach - evacuated 
    absorbed<-(evacuated*pval) * rAbsorb(temp,rAbd50,rAblrc,rAbmax,rAbmin)
    spent<-(Rmetabolism(temp,m_Int,m_Slope) * body)*OxyCal
    energy<-energy + absorbed - spent
    if(state==2){
      temp<-temp0+FRTdif
      sucess<-sample(c(1,0),size=1,p=c(pSuc,(1-pSuc)))
      if(sucess==1){
        stomach<-stomach+prey
      }
      state<-sample(c(1,2),size=1,p=c(pFR(stomach,FRs50,FRlrc),(1-pFR(stomach,FRs50,FRlrc))))
    }else{
      sucess<-0
      temp<-temp0-FRTdif
      state<-sample(c(1,2),size=1,p=c((1-pRF(stomach,RFs50,RFlrc)),pRF(stomach,RFs50,RFlrc)))
    }
    
    stateRec[i,1]<-state
    stateRec[i,2]<-temp
    stateRec[i,3]<-stomach
    stateRec[i,4]<-energy
    stateRec[i,5]<-sucess
  }
  return(stateRec)
}









out<-TempBudgetModel(prey=.3, pval=1200,
                     pSuc=.2,
                     temp=25,
                     body=1,
                     OxyCal=0.000324,
                     Time=1000,
                     FRs50=.7, FRlrc=.1,
                     RFs50=.02, RFlrc=.1,
                     rEvd50=20,rEvlrc=4,rEvmax=.1,
                     rAbd50=20,rAblrc=4,rAbmax=.8,rAbmin=.8,
                     m_Int=16.8272,m_Slope=.2315,
                     FRTdif=0)




size<-.5
maxRate<-.1
sizes<-NULL
for( i in 0:100){
  
  size <- size - (size * (i/100))
  sizes <- c(sizes,size)
}



y= exp(0.0923 - 0.0252*(0:100))











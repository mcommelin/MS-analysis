# Load and analyse the results generated from Masslinx

#### INITIALIZATION ####

# Packages 
  if(!require(readxl)){install.packages("readxl")}
  if(!require(ggplot2)){install.packages("ggplot2")}
  if(!require(ggplot2)){install.packages("grid")}
  if(!require(tidyverse)){install.packages("tidyverse")}

  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  rm(list=ls()) # cleaning console
  graphics.off() # cleaning plots

# set the compounds list
  Compounds=c("13C-Caffeine","Pendimethalin", "Chlorantranilliprole","Boscalid","Boscalid met.","Pyraclostrobin")

# List of the different batches : excel tabs
  Batchnames=c("07_24","07_27","07_28", "07_31")

# Create a list of Sample list per batch: 
  Samplelist<- vector("list", length(Batchnames)) 
  names(Samplelist)=Batchnames

# (optional) #Create a list of Coumpound list per batch: 
  Compoundlist<- list(Compounds,Compounds,Compounds, Compounds)
  names(Compoundlist)=Batchnames

# Create a Data Frame for raw results:
  Area_DF<-data.frame( Batch=factor(), Compound=factor(), Num=numeric(),      File_Name=factor(), Sample_text=factor(),
                     RT=numeric(),   Area=numeric(),    Sec.Area=numeric(), IR=numeric() )

# Create a Data Frame for Average Ion Ratio for each batch and each compounds:
  IR_av<-data.frame(matrix(0,length(Batchnames),length(Compounds)+1))
  names(IR_av)<-c("Batch",Compounds)
  IR_av$Batch=Batchnames


#### LOADING ####

#Set work directory
  setwd("C:/MassLynx/Default.pro/Data")
# Prepare and load the data into Area_DF

  for(i in 1:length(Batchnames)){   #Loading for each batch, from a different tab
    
    Area_i <-read_excel("LC_Area_Meso.xlsx",i,skip = 6) # Load tab number i into temporary Data Frame Area_i
    if (ncol(Area_i)>6){Area_i[,7: ncol(Area_i)]=NULL }    # Delete columns after Area.sec 

    colnames(Area_i)=c("Num", "File_Name", "Sample_text", "RT", "Area", "Sec.Area") # name of the columns, Rename the columns of the data frame
    
    Area_i$Num=as.numeric(Area_i$Num)   # Replace text in first column by NA
    Area_i=Area_i[!is.na(Area_i$Num),]  # Delete all lines with NA. in first column. Only lines with the samples are left. They are numbered in the first column
    
    Samplelist[[i]]=Area_i$Sample_text[1:max(Area_i$Num)] # Extract the name of the sample from the first compound results
    
    Area_i$Compound=rep(1)       # Create an additional column in the data frame for the name of the compound analised
    
  # Create a column "Compound" with the name of analized comound for each sample
    for(c in 1:length(Compoundlist[[i]])){
      Area_i$Compound[((c-1)*length(Samplelist[[i]])+1):(c*length(Samplelist[[i]]))]=Compoundlist[[i]][c]
    }
    
    Area_i$Batch=Batchnames[i]
    Area_DF=rbind(Area_DF,Area_i)
    
  } # end loading batches

# Create a column Type : Cal, QC, Std, Spike, Blk, Sample :
  Area_DF$Type="Sample"
  Area_DF$Type[grep("Cal", Area_DF$Sample_text)]="Cal"
  Area_DF$Type[grep("QC", Area_DF$Sample_text)]="QC"
  Area_DF$Type[grep("Std", Area_DF$Sample_text)]="Std"
  Area_DF$Type[grep("Blank soil", Area_DF$Sample_text)]="Blk_soil"
  Area_DF$Type[grep("Blank chemical", Area_DF$Sample_text)]="Blk"
  Area_DF$Type[grep("mq H2O", Area_DF$Sample_text)]="H2O"
  Area_DF$Type[grep("pike", Area_DF$Sample_text)]="Spike"

# Replace NA by 0 in RT, Area and sec.Area
  Area_DF$RT=as.numeric(Area_DF$RT)
  Area_DF$Area=as.numeric(Area_DF$Area)
  Area_DF$Sec.Area=as.numeric(Area_DF$Sec.Area)
  
  Area_DF$RT[is.na(Area_DF$RT)]=0
  Area_DF$Area[is.na(Area_DF$Area)]=0
  Area_DF$Sec.Area[is.na(Area_DF$Sec.Area)]=0

### ***Sample Meta data ####
  Sample_Meta=read_excel("LC_Area_Meso.xlsx","Sample_Meta",skip = 0)
  
#### ***13C-Caffeine ####
#  subset(Area_DF, Area_DF$Compound=="13C-Caffeine" & Area_DF$Type=="Sample" & Area_DF$Area==0)

# Delete the area in cal 0 to calculate the concentration in the samples : 
  #Area_DF$Area[ (Area_DF$Compound=="13C-Caffeine" & Area_DF$Sample_text=="Cal 0")]=0

# Remove 13C-Caffeine
  Area_DF=subset(Area_DF, Compound!= "13C-Caffeine")
  Compounds=c("Pendimethalin", "Chlorantranilliprole","Boscalid","Boscalid met.","Pyraclostrobin")
  Compoundlist<- list(Compounds,Compounds,Compounds,Compounds)

#### AREA CORRECTION (optional) #####

# Remove the background signal from the Area
# Substract the Area for blank sample (labeled "H2O")

  Area_DF$Bg_Area=0
  Area_DF$Bg_Sec.Area=0
  for(i in 1:length(Batchnames)){   
    for(c in 1:length(Compoundlist[[i]])){
      Area_DF$Bg_Area[((c-1)*length(Samplelist[[i]])+1):(c*length(Samplelist[[i]]))]=mean(subset( Area_DF,  Compound==Compoundlist[[i]][c] & Batch==Batchnames[i] & Type=="H2O" ,select = Area)[[1]])
      Area_DF$Bg_Sec.Area[((c-1)*length(Samplelist[[i]])+1):(c*length(Samplelist[[i]]))]=mean(subset( Area_DF,  Compound==Compoundlist[[i]][c] & Batch==Batchnames[i] & Type=="H2O" ,select = Sec.Area)[[1]])
    }
  } 

  Area_DF$Area_cor=Area_DF$Area#-Area_DF$Bg_Area
  Area_DF$Sec.Area_cor=Area_DF$Sec.Area#-Area_DF$Bg_Sec.Area
  


#### ION RATIO #### 

# Calculate Ion Ration and replace NA by 0
  Area_DF$IR=Area_DF$Sec.Area_cor/Area_DF$Area_cor
  Area_DF$IR[is.na(Area_DF$IR)]=0
  Area_DF$IR[is.infinite(Area_DF$IR)]=0

# Create a column for the Average Ion ration for each batch each compound 
  Area_DF$IR_av=0
# Create a column for to check the Ion ration range for each sample
  Area_DF$IR_ck=0


### *** Initial calculation of Av_IR ####

# Calculate average Ion Ration of Calibration curve and Standard in matrix for each batch and each compound
  for(i  in 1:length(Batchnames)){
    for(c in 1:length(Compoundlist[[i]])){
      # Average for "Cal" and "Std"
      # Output in both data frame IR_av and column Area_DF$IR_av
      IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]=mean(subset( Area_DF,  Compound==Compoundlist[[i]][c] & Batch==Batchnames[i] & Type %in% c("Cal","Std"),select = IR)[[1]])
      
      Area_DF$IR_av[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]= 
        IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]
    } #end Compound loop
  } #end Batch loop


# replace NA by 0
  IR_av[is.na(IR_av)]=0
  Area_DF$IR_av[is.na(Area_DF$IR_av)]=0


### *** IR check for IR_av<1 (optional) ####  

# #Re-calculate all IR for IR_av>1 : reverse Area and Sec.Area
# Area_DF$IR[Area_DF$IR_av>1]=1/ Area_DF$IR[Area_DF$IR_av>1] # Area_DF$Area[Area_DF$IR_av>1]/Area_DF$Sec.Area[Area_DF$IR_av>1]
# Area_DF$IR[is.na(Area_DF$IR)]=0   # replace NA by 0
# Area_DF$IR[is.infinite(Area_DF$IR)]=0   # replace Inf by 0
# 
# # Re-calculate IR_av with "Cal" and "Std", IR_av<1
#   for(i  in 1:length(Batchnames)){
#     for(c in 1:length(Compoundlist[[i]])){
#       # subset of Batch_i, Compound_c, IR >< 0
#         # SubArea_ic= subset(Area_DF,Area_DF$Batch==Batchnames[i]
#         #                  & Area_DF$Compound==Compoundlist[[i]][c]
#         #                  & Area_DF$IR!=0 )
#       # Average for "Cal" and "Std"
#       # Output in both data frame IR_av and column Area_DF$IR_av
#         IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]] =
#           # ( sum(SubArea_ic$IR[SubArea_ic$Type=="Cal"])+sum(SubArea_ic$IR[SubArea_ic$Type=="Std"]) )/
#           # ( length(SubArea_ic$IR[SubArea_ic$Type=="Cal"])+length(SubArea_ic$IR[SubArea_ic$Type=="Std"]) )
# 
#         Area_DF$IR_av[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]] =
#           IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]
#     } #end Compound loop
#   } #end Batch loop
# 
#   # replace NA by 0
#     IR_av[is.na(IR_av)]=0
#     Area_DF[is.na(Area_DF)]=0


### *** IR check in [0.5IR_av, 1.5IR_av ] ####  

# Re-calculate IR_av only with "Cal" and "Std" in [0.5IR_av, 1.5IR_av ]
  for(i  in 1:length(Batchnames)){
    for(c in 1:length(Compoundlist[[i]])){
     
      # Average for "Cal" and "Std"
      # Output in both data frame IR_av and column Area_DF$IR_av
      IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]=mean(subset( Area_DF,  Compound==Compoundlist[[i]][c] & Batch==Batchnames[i] & Type %in% c("Cal","Std")
                                                                          & Area_DF$IR > Area_DF$IR_av/3 & Area_DF$IR < Area_DF$IR_av*2, select = IR)[[1]])
      
      Area_DF$IR_av[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]] = 
        IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]
    } #end Compound loop
  } #end Batch loop  
  
  # replace NA by 0
  IR_av[is.na(IR_av)]=0


### *** IR check in [0.7IR_av, 1.3IR_av ] ####  

# Check that all IR are in [0.7IR_av, 1.3IR_av ] : 1=True 0=False
  Area_DF$IR_ck[Area_DF$IR<1.3*Area_DF$IR_av &
                  Area_DF$IR>0.7*Area_DF$IR_av]=1 

# Re-calculate IR_av only with "Cal" and "Std" in [0.7IR_av, 1.3IR_av ]
  for(i  in 1:length(Batchnames)){
    for(c in 1:length(Compoundlist[[i]])){
      # Average for "Cal" and "Std"
      # Output in both data frame IR_av and column Area_DF$IR_av
      IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]] = mean(subset( Area_DF,  Compound==Compoundlist[[i]][c] & Batch==Batchnames[i] & Type %in% c("Cal","Std")
                                                                            & Area_DF$IR_ck==1, select = IR)[[1]])
      
      Area_DF$IR_av[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]] = 
        IR_av[IR_av$Batch==Batchnames[i],Compoundlist[[i]][c]]
    } #end Compound loop
  } #end Batch loop  

# replace NA by 0
  IR_av[is.na(IR_av)]=0

# Check that all IR are in [0.7IR_av, 1.3IR_av ] : 1=True 0=False
  sum(Area_DF$IR_ck)
  Area_DF$IR_ck[Area_DF$IR<1.3*Area_DF$IR_av &
                  Area_DF$IR>0.7*Area_DF$IR_av]=1 
  sum(Area_DF$IR_ck)


#### CALIBRATION ####

#**CALIBATION CONCENTRATION [ng/mL] ###
# deviation of back-calculated concentration from true concentration within ?20%   
# Create Data Frame calibration : 
  Cal_DF<-Area_DF[Area_DF$Type=="Cal",]
#Concentration applied in cal [ng/mL]
  Cal_DF$Level=0
  Cal_DF$Level[grep("1", Cal_DF$Sample_text)]=1
  Cal_DF$Level[grep("5", Cal_DF$Sample_text)]=5
  Cal_DF$Level[grep("10", Cal_DF$Sample_text)]=10
  Cal_DF$Level[grep("25", Cal_DF$Sample_text)]=25
  Cal_DF$Level[grep("0.125", Cal_DF$Sample_text)]=0.125
  Cal_DF$Level[grep("0.25", Cal_DF$Sample_text)]=0.25
  Cal_DF$Level[grep("0.5", Cal_DF$Sample_text)]=0.5
  Cal_DF$Level[grep("2.5", Cal_DF$Sample_text)]=2.5


# Calibration concentration [ng/mg]
# Calculate the concentration in vial assuming that the area for "cal 5" represent 5ng/mL
  Cal_DF$Conc=0
  Cal_DF$Slope=0
  Cal_DF$R2=0
  
  Area_DF$Cal_Slope=0
  Area_DF$Cal_R2=0

  for(i  in 1:(length(Batchnames)-1)){
    for(c in 1:length(Compoundlist[[i]])){
      SubCal_ic= subset( Cal_DF, Batch==Batchnames[i] 
                         &  Compound==Compoundlist[[i]][c] )
      A=lm(Level~Area_cor, SubCal_ic)
      A=lm(Level~0 +Area_cor, SubCal_ic)
      
      Cal_DF$Slope[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]]=A[[1]][1]
      Cal_DF$R2[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]]=summary(A)$r.squared
  
    } #end Compound loop
  } #end Batch loop

  Cal_DF$Conc_lm=Cal_DF$Area_cor*Cal_DF$Slope
  Cal_DF$Recovery=Cal_DF$Conc_lm/Cal_DF$Level

# Check that all Recovery are in [0.8, 1.2 ] : 1=True 0=False
  Cal_DF$Recovery_ck=0
  Cal_DF$Recovery_ck[Cal_DF$Recovery<1.2 & Cal_DF$Recovery>0.8]=1

#***Remove cal 25  #### 
  Cal_DF$Slope_25=0
  Cal_DF$R2_25=0

  for(i  in 1:length(Batchnames)){
    for(c in 1:length(Compoundlist[[i]])){
      SubCal_ic= subset( Cal_DF, Batch==Batchnames[i] &
                           Compound==Compoundlist[[i]][c] &
                           Level != 25)
      A=lm(Level~Area_cor, SubCal_ic)
      A=lm(Level~0 +Area_cor, SubCal_ic)
      
      Cal_DF$Slope_25[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]]=A[[1]][1]
      Cal_DF$R2_25[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]]=summary(A)$r.squared
  
    } #end Compound loop
  } #end Batch loop

# *** Plot Calibration ####
 for(i  in 1:length(Batchnames)){
  for(c in 1:length(Compoundlist[[i]])){
    SubCal_ic= subset( Cal_DF, Batch==Batchnames[i] &
                         Compound==Compoundlist[[i]][c] & Level < 26)
    
    A=lm(Level~Area_cor, SubCal_ic)
    #A=lm(Level~0 +Area_cor, SubCal_ic)
    
    R2 <- as.character(as.expression( substitute(
      italic(R)^2~"="~r2,
      list(r2 = format(  summary(A)$r.squared , digits = 3) )
    )))

    PLOT=ggplot( SubCal_ic, aes(x=Area_cor, y=Level, color=Num))+
      geom_point()+
      #geom_abline(intercept =  A[[1]][1], slope =  A[[1]][2], size=1  )+
      geom_abline( slope =  A[[1]][2], intercept = A[[1]][1], size=1  )+
      geom_text(x = 1000, y = 10, hjust = 0 , parse = TRUE, colour="black", show.legend = FALSE, size=5,
                label =  R2  )+
      ggtitle(label = paste(Compoundlist[[i]][c],Batchnames[i], sep = "_"  ))+
      theme(
        # panel.background = element_rect(fill = "white", colour = "white"),
        # panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
        #plot.title =element_text(size=8),
        #axis.title.x = element_blank(), #remove the x label
        #legend.position="none",
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.x  = element_text(size=9),
        axis.text.y =  element_text(size=9))

    print(PLOT)
  } #end Compound loop
} #end Batch loop

  
  Cal_DF$Conc_lm_25=Cal_DF$Area_cor*Cal_DF$Slope_25
  Cal_DF$Recovery_25=Cal_DF$Conc_lm_25/Cal_DF$Level
  
  Cal_DF$Recovery[Cal_DF$Recovery==Inf]=0
  Cal_DF$Recovery_25[Cal_DF$Recovery_25==Inf]=0
  Cal_DF$Recovery[Cal_DF$Recovery==-Inf]=0
  Cal_DF$Recovery_25[Cal_DF$Recovery_25==-Inf]=0
  Cal_DF$Recovery[is.nan(Cal_DF$Recovery) ]=0
  Cal_DF$Recovery_25[is.nan(Cal_DF$Recovery_25) ]=0

# Check that all Recovery are in [0.8, 1.2 ] : 1=True 0=False
  Cal_DF$Recovery_ck_25=0
  Cal_DF$Recovery_ck_25[Cal_DF$Recovery_25<1.2 & Cal_DF$Recovery_25>0.8]=1

#summary

  Recovery_cal=function(x){
    Rec=0
    for (i in 1: length(x)){
      
      Rec=Rec+abs(1-x[i])
    }
    return(Rec)
  }
  
  Cal_Recovery= Cal_DF %>%
    group_by(Batch, Compound) %>%
    summarise(Recovery=Recovery_cal(Recovery), Recovery_25=Recovery_cal(Recovery_25), slope=mean(Slope), slope_sd=sd(Slope), Slope_25=mean(Slope_25))
  
  Cal_Recovery$R25= Cal_Recovery$Recovery-Cal_Recovery$Recovery_25

#### CONCENTRATION IN VIAL [ng/mL] ####
  # # Calculate the concentration in vial assuming that the area for the standard in blanck matrix represent 5ng/mL
  # # Substraction of the pic area for "Cal 0" to correct for noise pic.
  #     
  #   # Create a column for Sample_Concentration
  #   Area_DF$Conc_Vial=0
  #     #// start iteration at 12 because the calculation needs a std before the sample.
  #   for(s in 12:length(Area_DF$Num)){
  #     # Sample_Concentration = (Sample_Area - Cal 0_Area) *5 /mean(neighbouring Std)
  #     Area_DF$Conc_Vial[s]=(Area_DF$Area[s] - Area_DF$Area[s+2-match("Cal 0",Area_DF$Sample_text[(s+1):1])])*5/
  #       mean(c(Area_DF$Area[ s+1-match("Std", Area_DF$Type[s:(s-12)]) ],  
  #              Area_DF$Area[ s+match("Std", Area_DF$Type[(s+1):(s+13)]) ]))
  #   } #end sample loop
  #   #Last one Manual calculation :
  #     Area_DF$Conc_Vial[s]=(Area_DF$Area[s] - Area_DF$Area[s+2-match("Cal 0",Area_DF$Sample_text[(s+1):1])])*5/
  #       mean(c(Area_DF$Area[ s+1-match("Std", Area_DF$Type[s:(s-12)]) ], Area_DF$Area[ s]))
  #     
  #   Area_DF$Conc_Vial[Area_DF$Conc_Vial<0]=0
  #     
  #   std_DF<-Area_DF[Area_DF$Type=="Std",]
  
  Area_DF$Cal_methode=0
  for(i  in 1:length(Batchnames)){
    for(c in 1:length(Compoundlist[[i]])){
      if (Cal_Recovery$Recovery[Cal_Recovery$Batch==Batchnames[i] & Cal_Recovery$Compound==Compoundlist[[i]][c]] - 
          Cal_Recovery$Recovery_25[Cal_Recovery$Batch==Batchnames[i] & Cal_Recovery$Compound==Compoundlist[[i]][c]] <0){
        
        Area_DF$Cal_Slope[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]= unique(Cal_DF$Slope[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]])
        Area_DF$Cal_methode[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]=1
        #  Area_DF$Cal_R2[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]=summary(A)$r.squared
        #  
      } else { 
        Area_DF$Cal_Slope[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]= unique(Cal_DF$Slope_25[Cal_DF$Batch==Batchnames[i] & Cal_DF$Compound==Compoundlist[[i]][c]])
        Area_DF$Cal_methode[Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c]]=25
      }
    }
  }


Area_DF$Conc_Vial=Area_DF$Area_cor*Area_DF$Cal_Slope


#### DILUTION #### 

  Area_DF$Dilution=1
  Area_DF$Dilution[grep("_20",Area_DF$Sample_text)]=20
  Area_DF$Dilution[grep("_50",Area_DF$Sample_text)]=50
  Area_DF$Dilution[grep("_100",Area_DF$Sample_text)]=100
  Area_DF$Dilution[grep("_200",Area_DF$Sample_text)]=200

# A=subset(Area_DF,Area_DF$Dilution>1 & Area_DF$Compound=="13C-Caffeine" ) 

# *** Compare the dilution 

#### CONCENTRATION IN SAMPLE [ng/g] ####
# CORRECTING FACTOR
  Area_DF$Corr_Factor=4.4
  
  Area_DF$Conc_Soil=Area_DF$Conc_Vial*Area_DF$Corr_Factor*Area_DF$Dilution

#### METADATA ####
  Area_DF$Pest=0
  Area_DF$Plastic=0
  Area_DF$Sampling=0
  
  Plastic_lab=unique(Sample_Meta$plastic_lab)
  Pesticide_lab=unique(Sample_Meta$pesticide_lab)
  Sampling_num=unique(Sample_Meta$Sampling_num)
  
  for(pe in 1:length(  Plastic_lab)){
    Area_DF$Pest[Area_DF$Sample_text %in% Sample_Meta$Label[Sample_Meta$pesticide_lab== Pesticide_lab[pe]]]=Pesticide_lab[pe]
  }
  
  for(pl in 1:length(  Plastic_lab)){
    Area_DF$Plastic[Area_DF$Sample_text %in% Sample_Meta$Label[Sample_Meta$plastic_lab== Plastic_lab[pl]]]= Plastic_lab[pl]
  }
  
  for(s in 1:length(  Sampling_num)){
    Area_DF$Sampling[Area_DF$Sample_text %in% Sample_Meta$Label[Sample_Meta$Sampling_num== Sampling_num[s]]]= Sampling_num[s]
  }


#### RECOVERY : SPIKES ####
  # Concentration applied in Qc, Std and Spike
    Area_DF$Level=0
  # Area_DF$Level[grep("QC 0.5", Area_DF$Sample_text)]=0.5  # [ng/g]
  # Area_DF$Level[grep("QC 1", Area_DF$Sample_text)]=1      # [ng/g]
  # Area_DF$Level[grep("QC 2.5", Area_DF$Sample_text)]=2.5  # [ng/g]
  # Area_DF$Level[grep("QC 5", Area_DF$Sample_text)]=5      # [ng/g]
  # Area_DF$Level[grep("QC 10", Area_DF$Sample_text)]=10    # [ng/g]
  Area_DF$Level[grep("Spike", Area_DF$Type)]=5            # [ng/g]
  Area_DF$Level[grep("Std", Area_DF$Type)]=5              # [ng/mL]

# Create a column for recovery
# Area_DF$Recovery=0

# # QC recovery = Area_QC/Level_QC
#   Area_DF$Recovery[Area_DF$Type=="QC"]=Area_DF$Conc_Soil[Area_DF$Type=="QC"]/Area_DF$Level[Area_DF$Type=="QC"]
#   Area_DF$Recovery[Area_DF$Sample_text=="QC 0"]=0
#   
# Spike recovery = (Area_Spike - Area_Sample)/Level_Spike
  Area_DF$Recovery[Area_DF$Type=="Spike"] = (Area_DF$Conc_Soil[Area_DF$Type=="Spike"] - 
                                             Area_DF$Conc_Soil[which(Area_DF$Type=="Spike")-1])/
  Area_DF$Level[Area_DF$Type=="Spike"]

  unique(Area_DF$Sample_text[which(Area_DF$Type=="Spike")])
  unique(Area_DF$Sample_text[which(Area_DF$Type=="Spike")-1])
# Warning Sample_A/std/Sample_Spike 
  Area_DF$Recovery[Area_DF$File_Name=="TQS3_200727_081"] = (Area_DF$Conc_Soil[Area_DF$File_Name=="TQS3_200727_081"] - 
                                                            Area_DF$Conc_Soil[which(Area_DF$File_Name=="TQS3_200727_081")-2])/
  Area_DF$Level[Area_DF$File_Name=="TQS3_200727_081"]



# Warning Sample_A <-> Sample_Spike
  Area_DF$Recovery[Area_DF$File_Name=="TQS3_200724_096"]
  # Area_DF$Recovery[Area_DF$File_Name=="TQS3_200724_096"] = (Area_DF$Conc_Soil[which(Area_DF$File_Name=="TQS3_200724_096")-1] - 
  #                                                             Area_DF$Conc_Soil[which(Area_DF$File_Name=="TQS3_200724_096")])/
  #                                                    Area_DF$Level[Area_DF$File_Name=="TQS3_200724_096"]

# Check that all Recovery are in [0.75, 1.25 ] : 1=True 0=False
  Area_DF$Recovery_ck=0
  Area_DF$Recovery_ck[Area_DF$Recovery<1.25 & Area_DF$Recovery>0.75]=1 
  
  Area_DF$Sample_text[Area_DF$Recovery_ck==1]
  C=Area_DF[Area_DF$Recovery_ck!=1 &Area_DF$Type=="Spike",  c("Sample_text","File_Name","Conc_Vial","Recovery","Pest","Compound","Sampling","Batch")]

# 
# #Load Cal_DF results in Area_DF
# Area_DF$Level[Area_DF$Type=="Cal"]=Cal_DF$Level
# Area_DF$Recovery[Area_DF$Type=="Cal"]=Cal_DF$Recovery
# Area_DF$Recovery_ck[Area_DF$Type=="Cal"]=Cal_DF$Recovery_ck


####SUMMARY####
  Pest_sum = subset(Area_DF, Area_DF$Type=="Sample") %>% group_by(Compound) %>%  # Pest,Plastic, Sampling
    summarise(Min=min(Conc_Vial), Mean=mean(Conc_Vial), Max= max(Conc_Vial), Sd=sd(Conc_Vial))

#*** Rerquired dilution####
  #Identify samples that need to be diluted
  ToDilute=Area_DF[Area_DF$Conc_Vial>25 & Area_DF$Type=="Sample" , c("Sample_text","File_Name","Batch","Compound","Conc_Vial","Pest","Sampling","Batch")]
  length(unique( ToDilute$Sample_text))
  
  
  ToDilute$Dilution=0
  ToDilute$Dilution[ToDilute$Conc_Vial<50]=10
  ToDilute$Dilution[ToDilute$Conc_Vial>50 & ToDilute$Conc_Vial<100 ]=25
  ToDilute$Dilution[ToDilute$Conc_Vial>100 & ToDilute$Conc_Vial<150 ]=50
  ToDilute$Dilution[ToDilute$Conc_Vial>150 ]=100
  D=ToDilute[ToDilute$Sample_text %in% ToDilute$Sample_text[duplicated(ToDilute$Sample_text)], ]
  
  
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_031"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_029"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_083"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_068"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_076"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_078"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200728_040"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_074"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_079"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_066"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_028"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200728_030"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_060"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_057"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_039"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_032"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_108"]=25
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_019"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200724_045"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_099"]=50
   ToDilute$Dilution[ToDilute$File_Name=="TQS3_200727_095"]=25
 
 ToDilute$Dilution= ToDilute$Dilution*2
 
 ToDilute$New_conc=ToDilute$Conc_Vial/ToDilute$Dilution*2
 
 ToDilute$Volume_tot= 500
 ToDilute$Extract=500/ ToDilute$Dilution
 ToDilute$Solvent=500-500/ ToDilute$Dilution
 
 ToDilute_sample=  ToDilute[!duplicated( ToDilute$File_Name), ]
 
 ToDilute_sample=subset(ToDilute_sample, select = c(Sample_text,Batch,Dilution,Extract, Solvent))
 
 #*** check dilutions ####

  ToDilute$Sample_text_Dilution= "s"
  ToDilute$File_Name_Dilution= "f"
  ToDilute$Conc_Vial_Dilution= 0
  
 for (i in 1:nrow(ToDilute)){
   ToDilute$Sample_text_Dilution[i]= paste(ToDilute$Sample_text[i],ToDilute$Dilution[i],sep = "_")
   ToDilute$File_Name_Dilution[i]= Area_DF$File_Name[Area_DF$Sample_text  ==  ToDilute$Sample_text_Dilution[i] & Area_DF$Compound==ToDilute$Compound[i]]
   ToDilute$Conc_Vial_Dilution[i]=Area_DF$Conc_Vial[Area_DF$Sample_text  ==  ToDilute$Sample_text_Dilution[i] & Area_DF$Compound==ToDilute$Compound[i]]
 }
 
  ToDilute$dilute_calc_ratio= ToDilute$Conc_Vial_Dilution/ToDilute$New_conc
  
 #write.csv(ToDilute, "To_be_Diluted3.csv")
 #write.csv(ToDilute_sample, "To_be_Diluted_sample.csv")
 # write.csv(Area_DF, "Area_DF.csv") 

#### LODs #### 
#MAX( lowest QC 2sd area > 50 & QC recovery in [0.75, 1.25 ] & IR in [0.75, 1.25 ] )

# Create a Data Frame for LOD [ng/g] for each batch and each compounds:
LOD<-data.frame(matrix(0,length(Batchnames),length(Compounds)+1))
names(LOD)<-c("Batch",Compounds)
LOD$Batch=Batchnames  

# for QC [ng/g] and  Calibration [ng/mL]
LOD_QC=LOD
LOD_Cal=LOD
for(i  in 1:length(Batchnames)){
  for(c in 1:length(Compoundlist[[i]])){     
    LOD_QC[i,Compoundlist[[i]][c]]= min(Area_DF$Level[Area_DF$Type=="QC" & Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c] &
                                                        Area_DF$Area>75   &  Area_DF$Recovery_ck==1      &   Area_DF$IR_ck==1]     )   
    
    LOD_Cal[i,Compoundlist[[i]][c]]= min(Area_DF$Level[Area_DF$Type=="Cal" &   Area_DF$Batch==Batchnames[i] & Area_DF$Compound==Compoundlist[[i]][c] & 
                                                         Area_DF$Area>75   &   Area_DF$Recovery_ck==1     &   Area_DF$IR_ck==1]     ) 
  } #end Compound loop
} #end Batch loop

LOD_QC[LOD_QC==Inf]=0  #mins that there are no QCs
LOD_Cal[LOD_Cal==Inf]=0
LOD<-LOD[1,2:length(LOD)]
Area_DF$LOD=10


#!?MAY BE NOT THE BEST WAY  
LOD_QC[3:4,]=LOD_Cal[3:4,]
for(c in 1:length(Compounds)){
  LOD[Compounds[c]]=max(c(LOD_QC[,Compounds[c]]))
  
  Area_DF$LOD[Area_DF$Compound==Compounds[c]]= LOD[1,Compounds[c]]
} #end Compound loop


# LOD Spinosyn-A


# Account for the repetition of measure  : 16_07 <-> 10_08 ???

# Round for LOD from Cal : 
LOD[LOD==0.55]=0.5
LOD[LOD==1.1]=1
LOD[LOD==2.2]=2.5
LOD[LOD==4.4]=5
LOD[LOD>10]=10
# =>  LOD in [ng/g]  <=      

####*** Check LOD ####
Area_DF$LOD_ck=0
Area_DF$LOD_ck[(Area_DF$LOD*Area_DF$Dilution)<Area_DF$Conc_Soil]=1

#### RESULTS, selection  ####




# Load Area_DF into Results to delete un-necessary samples and remove B's in results 
# # Results=Area_DF[Area_DF$Type=="Sample",]
Results=Area_DF[(!grepl(" B", Area_DF$Sample_text))&Area_DF$Type=="Sample",]


## *** NON SIGNIFICANT MEASURE ####
#Er -1 : 0<Content<LOD
#ER -2 : Content>LOD but out of linearity
#ER -3 : Dilution required but not available

# Give -1 for Content lower than LOD
#Plastic_positive=Plastic_content
Results2$Conc_Soil_positive=Results$Conc_Soil

for( s in 1:nrow(Results)){
  if( Results$Conc_Soil[s] < Results$LOD[s] * Results$Dilution[s]  & Results$Conc_Soil[s]!=0 ){  #!?!? USE Area_plastic$Content[s] < Area_plastic$LOD[s] * Area_plastic$Dilution[s] | Area_plastic$IR_ck[s]==0 Check with PAUL the condition for IR ION RATIO : "St Mix P a" "07_27" is FALSE
    Results$Conc_Soil_positive[s]=-1
  }
  
  
  
  #if Conc_Vial>26 & if there is in fact a dilution available => asign -2 
  if( Results$Conc_Vial[s]>26 &  length (Results$Sample_text[(grepl(Results$Sample_text[s],Results$Sample_text)) & Results$Compound=="Propamocarb"]) >2 ){
    Results$Conc_Soil_positive[s]=-2
  }
} 



# 1. Remove duplicates : calculate an average of A and B




## *** std with AB table : StAB ####

# table of Bs 
#B=Area_DF[grepl(" B", Area_DF$Sample_text)&Area_DF$Compound=="13C-Caffeine",]
B=Area_DF[grepl(" B", Area_DF$Sample_text),]
# table of As
A=Area_DF[which( grepl(" B", Area_DF$Sample_text) )-1, ]
A[A$Type!="Sample",]=Area_DF[which( Area_DF$File_Name==A$File_Name[A$Type!="Sample"] )-1, ]

# # DupAB= rbind(A,B)
StAB=subset(A, select=-c(RT, Area, Sec.Area, IR, Conc_Vial, Recovery, Recovery_ck, Level))
StAB$Sample_text_B=B$Sample_text
StAB$Conc_Soil_B=B$Conc_Soil
StAB$Mean_Conc_Soil=rowMeans(StAB[,c("Conc_Soil","Conc_Soil_B")])
StAB$Std=apply(StAB[,c("Conc_Soil","Conc_Soil_B")],1, sd, na.rm = TRUE)

# Asign Conc_Soil mean AB to A's //!\\  then Area. Sec>Area, IR... are wrong!
#  # Results$Conc_Soil[grepl(" A", Area_DF$Sample_text),]=StAB$Mean_Conc_Soil
for (i in 1:nrow(StAB)){
  Results$Conc_Soil[Results$File_Name==StAB$File_Name[i] & Results$Compound==StAB$Compound[i]]=StAB$Mean_Conc_Soil[i]
}

#### **** remove plastic ad samples ####
Results=Results[(!grepl("LDPE", Results$Sample_text)) & (!grepl("PAC", Results$Sample_text))&
                  (!grepl("Bio", Results$Sample_text))  & (!grepl("mix", Results$Sample_text))&
                  (!grepl("Mix", Results$Sample_text)) ,]


#### **** remove Dilutions and replicats ####

# 1. asign a sample name (sample_mum) to each measure of plastic_content
# 2. choose a unique value content/Num/compound based on the Conc in vial in Area_Df using "sample file as an index.
#2.1 create the table of dilution 
#2.2 choose the lowest dilution with conc_Vial<25.25ng/mL /\/\/ USe a variable to check the lineartity of Cal 25.
#save the File Name used for each Mum/compound: Plastic_ad_FN


Results_red=data.frame(matrix(0, nrow(Sample_Meta)*length(Compounds),4))  # Results with only one value per sample : \replicats. Dilution
colnames(Results_red)=c("Sample_text", "Compounds", "FN", "Conc_Soil")

# Conc_Soil : Data.Frame of required samples x compounds
Conc_Soil=data.frame(matrix(0,nrow(Sample_Meta),length(Compounds)),row.names=Sample_Meta$Sample_text)
colnames(Conc_Soil)=Compounds

# Warning !!! c(19:36,89:90,37:44) <=> Tricky way to select the sample text !!!
#Plastic_ad = subset(Plastic_content, Plastic_content$File_Name %in% Plastic_list$File_Name[c(19:36,89:90,37:44)], select = -c(File_Name,Batch) ) # plastic ad [ng]

i=1
for (s in 1:nrow(Sample_red)){      # For all required samples
  for( c in 1:length(Compounds)){     # For all compounds
    
    Results_red$Sample_text[i]=Sample_red$Sample_text[s]
    Results_red$Compounds[i]=Compounds[c]
    
    #subset all dilutions for sample s 
    C= subset(Results, grepl(Sample_red$Sample_text[s], Results$Sample_text )& Results$Compound==Compounds[c] )
    
    
    #if the smalest Conc_Vial of the highest dilution is inded smaler than 25.25
    if (min(C$Conc_Vial[C$Dilution==max(C$Dilution)])<25.25){ 
      # Imput the results of the lowest dilution / Conc_Vial<25.25
      Results_red$Conc_Soil[i]= max(C$Conc_Soil_positive[ C$Dilution*C$Corr_Factor== min( C$Dilution[C$Conc_Vial<25.25] * C$Corr_Factor[C$Conc_Vial<25.25] ) ])
      
      Results_red$FN[i]= max(C$File_Name[ C$Dilution*C$Corr_Factor== 
                                            min( C$Dilution[C$Conc_Vial<25.25] * C$Corr_Factor[C$Conc_Vial<25.25] ) ])
      #give all the values if they are egual so take the first one
      # Results_red$Compound[s]=  
      Conc_Soil[ s,Compounds[c] ] = max(C$Conc_Soil_positive[ C$Dilution*C$Corr_Factor== 
                                                                min( C$Dilution[C$Conc_Vial<25.25] * C$Corr_Factor[C$Conc_Vial<25.25] ) ])
      
      
    } else {     
      Results_red$Conc_Soil[i]  -3 
      Conc_Soil[ s,Compounds[c] ] =  -3 
      Results_red$FN[i] =  -3 
    }
    
    i=i+1
  } #end Compound loop
} #end sample loop


#### *** Check data ####


mySummary <- function(x) {
  c(#min = min(x),
    min_pos = min(x[x>0]),
    max = max(x), 
    #mean = mean(x), 
    median_pos  = median(x[x>0]), 
    std_pos  = sd(x[x>0]) )
  
  #var=var(x))
}

A=tapply(Results_red$Conc_Soil, Results_red$Compounds, mySummary)
DATA_SUM=as.data.frame(do.call(rbind, A))
DATA_SUM[DATA_SUM$max==0,c("min_pos","median_pos","std_pos")]=0
DATA_SUM$Compound=row.names(DATA_SUM)
DATA_SUM$LOD=100
for (i in 1:nrow(DATA_SUM)){
  DATA_SUM[i,"LOD"]=LOD[DATA_SUM$Compound[i]][[1]]
}

DATA_SUM[DATA_SUM$min_pos<DATA_SUM$LOD,]


#### *** Export ####



write.csv(Conc_Soil, "Conc_soil_LC_10_07_19.csv")

# options(digits=2)

#### *** Count residues ####
## Conc_soil was a DF of just samples*compounds, Indexing ad index for lines:
# //!\\ Check LOD! Why there are contents < 0.5 ng/mg ?! //!\\ 
Conc_Soil[Conc_Soil<0]=0
Conc_Soil$`13C-Caffeine`=0
Conc_Soil$Num_residues=0
for (s in 1:nrow(Conc_Soil)){
  Conc_Soil$Num_residues[s]=length( Conc_Soil[s,(Conc_Soil[s,Compounds]>0.5) ] ) 
  # length( Conc_Soil[s,(Conc_Soil[s,Compounds]>0.5) ] )                # I don't know why it counts Sum as part of Compounds??!
  # length( Conc_Soil[s,subset(Conc_Soil[s,], select = Compounds)>0.5]) # I don't know why it counts Sum as part of Compounds??!
}

### *** Indexing, Sum####
## as.factor or as.numeric() ??



#Conc_Soil$Sum=rowSums(Conc_Soil)   #It worked but not anymore ??!?
Conc_Soil$Sum=rowSums(Conc_Soil[,Compounds])  

Conc_Soil$Sample_text=Sample_red$Sample_text
Conc_Soil$num=as.numeric( gsub("[^[:digit:]]","",Conc_Soil$Sample_text) ) #remove letters in Sample_text
Conc_Soil$duplicat=gsub("[[:digit:]]","",Conc_Soil$Sample_text) #remove numbers in Sample_text
Conc_Soil$Coop=as.numeric( floor(Conc_Soil$num/10000) ) #1 to 4
Conc_Soil$Farm=as.numeric( floor(Conc_Soil$num/1000) ) #11 to 43
Conc_Soil$Parcel=as.numeric( floor(Conc_Soil$num/100) ) #111 to 433
Conc_Soil$Layer=as.numeric(  Conc_Soil$num%%10 ) # remainder of euclidian division by 10 : 1 (0-10cm) or 0 (10-30) 

## Results_red
Results_red$num=as.numeric( gsub("[^[:digit:]]","",Results_red$Sample_text) ) #remove letters in Sample_text
Results_red$duplicat=gsub("[[:digit:]]","",Results_red$Sample_text) #remove numbers in Sample_text
Results_red$Coop=as.numeric( floor(Results_red$num/10000) ) #1 to 4
Results_red$Farm=as.numeric( floor(Results_red$num/1000) ) #11 to 43
Results_red$Parcel=as.numeric( floor(Results_red$num/100) ) #111 to 433
Results_red$Layer=as.numeric(  Results_red$num%%10 ) # remainder of euclidian division by 10 : 1 (0-10cm) or 0 (10-30) 

## Sample_red
Sample_red$num=as.numeric( gsub("[^[:digit:]]","",Sample_red$Sample_text) ) #remove letters in Sample_text
Sample_red$duplicat=gsub("[[:digit:]]","",Sample_red$Sample_text) #remove numbers in Sample_text
Sample_red$Coop=as.numeric( floor(Sample_red$num/10000) ) #1 to 4
Sample_red$Farm=as.numeric( floor(Sample_red$num/1000) ) #11 to 43
Sample_red$Parcel=as.numeric( floor(Sample_red$num/100) ) #111 to 433
Sample_red$Layer=as.numeric(  Sample_red$num%%10 ) # remainder of euclidian division by 10 : 1 (0-10cm) or 0 (10-30) 



#*** Adding Coumpounds info####

Results_red$Pest_use="-" # INsecticide/ Herbicide/Fungicide
Results_red$Pest_group="-" #Substance group 

for (c in 1:length(Compounds)){
  Results_red$Pest_use[Results_red$Compounds==Compounds_xl$name[c]]=Compounds_xl$usage[c]
  Results_red$Pest_group[Results_red$Compounds==Compounds_xl$name[c]]=Compounds_xl$Substance_group[c]
}

# ### *** mean/ Coop /Layer####
# Coop_Layer= data.frame (matrix(0,8,length(Compounds)+4))
# names(Coop_Layer)<-c("Coop","Layer","Sum_pest","Num_pest", Compounds) #"Sum_pest_ng/g", "Num_pest" mean>0.5ng/g
# for (c in 1:4 ){
#   if (c==1|c==4){l=1
#   Coop_Layer$Coop[(c*2-l)]=c
#   Coop_Layer$Layer[(c*2-l)]=l
#   Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ]= colMeans( subset(Conc_Soil, Conc_Soil$Coop==c & Conc_Soil$Layer==l, select=Compounds ))
#   Coop_Layer$Sum_pest[(c*2-l)]=sum(Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ] )
#   Coop_Layer$Num_pest[(c*2-l)]= length( Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ] [Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ]>0.5]) #"Num_pest" mean>0.5ng/g
#   } else {
#   for (l in 1:0){
#     Coop_Layer$Coop[(c*2-l)]=c
#     Coop_Layer$Layer[(c*2-l)]=l
#     Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ]= colMeans( subset(Conc_Soil, Conc_Soil$Coop==c & Conc_Soil$Layer==l, select=Compounds ))
#     Coop_Layer$Sum_pest[(c*2-l)]=sum(Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ] )
#     Coop_Layer$Num_pest[(c*2-l)]= length( Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ] [Coop_Layer[ (c*2-l),5:ncol(Coop_Layer) ]>0.5]) #"Num_pest" mean>0.5ng/g
#    # Coop_Layer$Num_pest[(c*2-l)]=length(subset(Coop_Layer[(c*2-l),4:ncol(Coop_Layer), >1])
#   }
#   }
# }


#  write.csv(Coop_Layer, "Coop_Layer_LC_23_1.csv")


### *** mean/ parcel /Layer####




#### GRAPH Sector ####
### ***Graph_Sector data.frame ####

# Results_red$sector : Define the sectors

# Results_red$Sector="Other"
# 
# Results_red$Sector[Results_red$Compounds=="Boscalid"]="Boscalid"
# Results_red$Sector[Results_red$Compounds=="Oxyfluorfen"]="Oxyfluorfen"
# Results_red$Sector[Results_red$Compounds=="Pendimethalin"]="Pendimethalin"
# Results_red$Sector[Results_red$Compounds=="Chlorantraniliprole"]="Chlorantraniliprole"

Sector=c("Other","Boscalid","Oxyfluorfen","Pendimethalin","Chlorantraniliprole")
Results_red$Sector=Sector[1]      # initialyzed at "Other"

for (c in 2:(length(Sector))){
  Results_red$Sector[Results_red$Compounds==Sector[c]]=Sector[c]
}

Results_red$Conc_Soil[Results_red$Conc_Soil<0]=0
Results_red$Conc_Soil[Results_red$Compounds=="13C-Caffeine"]=0
Graph_Sector=Results_red[Results_red$Sector!="Other",] #Grapf Sector is asigned all sectores that don't need to be aggregated, summed

# 
# n=nrow(Graph_Sector)
# Graph_Sector[(n+1):(n+nrow(Sample_red)),]<-0 # Graph rows at the end of the data frame to add the sector(s) that need aggregation, sum
# 
# for (s in 1:nrow(Sample_red)) { #We want to add value to every sample
#   
#   Graph_Sector$Sample_text[n+s]=Sample_red$Sample_text[s] 
#   Graph_Sector$Sector[n+s]=Sector[1]
#   Graph_Sector$Conc_Soil[n+s]= sum( Results_red$Conc_Soil[Results_red$Sector==Sector[1] & Results_red$Sample_text==Sample_red$Sample_text[s] ] )
# }


# Aggregate the data in Sector "Other" by adding the conc 
Other=data.frame(matrix(0,nrow(Sample_red),length(Graph_Sector)))
colnames(Other)=colnames(Graph_Sector)                
Other$Sector=Sector[1]

for (s in 1:nrow(Sample_red)) { #We want to add value to every sample
  Other[s,c("Sample_text","Coop", "Farm", "Parcel", "Layer")]=as.numeric(Sample_red[s,c("Sample_text","Coop", "Farm", "Parcel", "Layer")])
  Other$Conc_Soil[s]= sum( Results_red$Conc_Soil[Results_red$Sector==Sector[1] & Results_red$Sample_text==Sample_red$Sample_text[s] ] )
}

# Add the secot "Other" at the end of  Graph_Sector
Graph_Sector=rbind(Graph_Sector,Other)

### ***Order #### 
# Sum
Graph_Sector$Sum=0
Graph_Sector$Sample_text=as.factor(Graph_Sector$Sample_text)

for (s in 1 : nrow(Sample_red)){ #We want to add value to every sample
  
  Graph_Sector$Sum[Graph_Sector$Sample_text==Sample_red$Sample_text[s] ]=sum ( Graph_Sector$Conc_Soil[Graph_Sector$Sample_text==Sample_red$Sample_text[s] ])
  
}

# Order
#Order from small to high Sum of concentration
Graph_Sector= Graph_Sector[with( Graph_Sector, order(Sum)), ]
#Order from Organic to conventional
Graph_Sector$Management="Conventional"
Graph_Sector$Management[Graph_Sector$Coop <3]="Organic"
Graph_Sector= Graph_Sector[with( Graph_Sector, order(Sum)), ]

max(Graph_Sector$Sum[Graph_Sector$Coop<3]) < min(Graph_Sector$Sum[Graph_Sector$Coop>2])

# Rank

Graph_Sector=subset(Graph_Sector, Graph_Sector$Layer==1 & Graph_Sector$Coop %in% c(1,2,3,4))

#  add a rank number for all the sector of a specific sample, Rank according to the 
Graph_Sector$Rank=0

for (s in 1 : (nrow(Graph_Sector)/5)){ #We want to add value to every sample
  
  
  Graph_Sector$Rank[Graph_Sector$Sum==Graph_Sector$Sum[Graph_Sector$Sector=="Other"][s] &  #samples with Sum = Sum^s
                      Graph_Sector$Rank==0 &
                      Graph_Sector$Sample_text==
                      Graph_Sector$Sample_text[Graph_Sector$Sum==Graph_Sector$Sum[Graph_Sector$Sector=="Other"][s] & Graph_Sector$Rank==0 ][1] #Select the first one of the sample text by applying the same conditions 
                    ]=s
}


### ***Graph_Sector PLOT ####



Data=subset(Graph_Sector, select=c(Rank, Conc_Soil,Sector))
Data$Sector=as.factor(Data$Sector)
Data$Conc_Soil=Data$Conc_Soil/1000  # ng/g -> mg/kg
Data$Sector=factor(Data$Sector , levels(Data$Sector)[c(3,1,2,4,5)]) #ordering sectors 

cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#999999", "#CC79A7") #colors

ggplot(Data, aes(x=Rank, y=Conc_Soil, fill=Sector)) + 
  geom_area()+
  xlab("Cartagena agricultural topSoils (0-10cm) (n=108)") +
  ylab("Pesticide Content mg/kg")+
  scale_fill_brewer(palette="Set1")+
  geom_vline(xintercept=mean(c(min(Data$Rank),max(Data$Rank))))+
  annotate("text", x=25, y=0.35, label= "Organic") + 
  annotate("text", x=75, y=0.35, label= "Conventional") 


#### PIE Chart ####

Pie_data=subset(Conc_Soil,Conc_Soil$Layer==0 & Conc_Soil$Coop %in% c(1,2), select=Num_residues)

Pie_residues=c( length(Pie_data[Pie_data$Num_residues>10,1]),
                length(Pie_data[(Pie_data$Num_residues<=10 & Pie_data$Num_residues>5),1] ),
                length(Pie_data[(Pie_data$Num_residues<=5 & Pie_data$Num_residues>1),1] ),
                length(Pie_data[Pie_data$Num_residues==1,1] ),
                length(Pie_data[Pie_data$Num_residues==0,1] )  )
Pie_residues=round(Pie_residues*100/nrow(Pie_data))

Pie_percent=Pie_residues
Pie_percent[1]=paste(Pie_residues[1],"%",sep="")
Pie_percent[2]=paste(Pie_residues[2],"%",sep="")
Pie_percent[3]=paste(Pie_residues[3],"%",sep="")
Pie_percent[4]=paste(Pie_residues[4],"%",sep="")
Pie_percent[5]=paste(Pie_residues[5],"%",sep="")

pie( Pie_residues, labels = Pie_percent, main = "Organic 10-30cm Samples",col = cbbPalette)
legend( c(1,1) , c(">10 residues","6-10 residues" ,"2-5 residues","1 residue","No residues > LOQ"), cex = 1,
        fill = cbbPalette)

# Pie_data[1]=paste("(",(Pie_data[1]),"%)",sep="")
# Pie_data[2]=paste("(",(Pie_data[2]),"%)",sep="")
# Pie_data[3]=paste("(",(Pie_data[3]),"%)",sep="")
# Pie_data[4]=paste("(",(Pie_data[4]),"%)",sep="")
# Pie_data[5]=paste("(",(Pie_data[5]),"%)",sep="")



# Conc_Soil$SectorPie=Conc_Soil$Num_residues
# Conc_Soil$SectorPie[Conc_Soil$Num_residues>10]=">10 residues"
# Conc_Soil$SectorPie[Conc_Soil$Num_residues<=10 & Conc_Soil$Num_residues>5]="6-10 residues"                    
# Conc_Soil$SectorPie[Conc_Soil$Num_residues<=5 & Conc_Soil$Num_residues>1]="2-5 residues"
# Conc_Soil$SectorPie[Conc_Soil$Num_residues==1]="1 residue"                   
# Conc_Soil$SectorPie[Conc_Soil$Num_residues==0]="No residues > LOQ"                   

# # Pie Plot with labels inside : 
# 
# #Pie chart without labels
# pie = ggplot(Conc_Soil, aes(x="", fill=factor(SectorPie))) + 
#   geom_bar(width=1) +
#   coord_polar(theta="y",start=0) +
#   geom_text(aes(y="",label="")) +
#   xlab("") + ylab("") +
#   theme_bw() + 
#   theme(legend.position = "none",
#         panel.grid.major = element_line(color="grey60"),
#         panel.border=element_blank())
# pie
# 
# #Add labels with % in pie chart
# library(grid)
# grid.text(paste(">10 residues"   ,Pie_data[1]  ,sep=" "), x=unit(0.49, "npc"), y=unit(0.91, "npc"), gp=gpar(fontsize=14, col="black"), rot = 00)
# grid.text(paste("6-10 residues"  ,Pie_data[2]  ,sep=" "), x=unit(0.64, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=14, col="black"), rot = 00)
# grid.text(paste("2-5 residues"   ,Pie_data[3]  ,sep=" "), x=unit(0.55, "npc"), y=unit(0.30, "npc"), gp=gpar(fontsize=14, col="black"), rot = 00)
# grid.text(paste("1 residue"        ,Pie_data[4],sep=" "), x=unit(0.35, "npc"), y=unit(0.55, "npc"), gp=gpar(fontsize=14, col="black"), rot = 00)
# grid.text(paste("No residues>LOQ",Pie_data[5],sep=" "), x=unit(0.38, "npc"), y=unit(0.72, "npc"), gp=gpar(fontsize=14, col="black"), rot = 00)
# 
# 
# 
# 

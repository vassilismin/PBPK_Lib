library(deSolve)
library(ggplot2)

# Set the path to load the necessary data 
setwd("C:/Users/vassi/Documents/GitHub/PBPK_Lib/Data")

dose_kg <- 10 # mg/kg rat body
mass <- 250 # g  
dose <- dose_kg*mass/1000 # mg TiO2

# Load raw data from paper Xie et al.2011
df <- openxlsx::read.xlsx("TiO2_iv_rat.xlsx", sheet = 1, 
                          colNames = T, rowNames = T) 
# TiO2 NPs %ID/g of tissue  (Table 1)
excretion <- openxlsx::read.xlsx("Cummulative_Excretion.xlsx", sheet = 1, 
                                 colNames = T, rowNames = F) 
# accumulated excretory rate, expressed as %ID

excretion <- excretion[,c(2:3)]

# Transform to (mg of NPs)/(g of tissue)
df <- (df/100)*dose
excretion <- (excretion/100)*dose

####################################################################
###################             DATA             ###################
####################################################################


### Important!!! each compartment has a specific index vectors Tissue_fractions,
# Regional_flow_fractions, Capillary_fractions and cannot be changed

# The index of each compartment:
#Rest of Body (rob) --> 1
#Heart (ht) --> 2
#Kidneys (ki) --> 3
#Brain (br) --> 4
#Spleen (spl) --> 5
#Lungs (lu) --> 6
#Liver (li) --> 7
#Uterus (ut) --> 8
#Bone (bone) --> 9
#Adipose (ad) --> 10
#Skin (skin) --> 11
#Muscles (mu) --> 12
#Gastrointestinal track (GIT) --> 13


####################
### User's INPUT ###
####################
#### If any of these compartments don not exist in pbpk, just give it the value 
#### NA in compartments vector, example: "Heart" = NA and it will remove it 
#### from the equilibriums and the corresponding V_tis, V_cap, Q will be equal
#### to NA.


compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", 
                      "Brain"= NA, "Spleen"="Spleen", "Lungs"="Lungs",
                      "Liver"="Liver", "Uterus"=NA, "Bone"="Bone", 
                      "Adipose"=NA, "Skin"=NA, "Muscles"=NA, "GIT"="GIT") 
#used as input in function, compartments that are used in pbpk


#####################################
### Function to create Parameters ###
#####################################
create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Bone"="Bone", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles",
                    "GIT"="GIT") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*w^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat physiological parameters.xlsx", sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  rownames(fractions) <- all_comps
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  #Macrophage content as fraction tissue volume for each tissue/organ
  Macrophage_fractions <- fractions[,4] 
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  W_macro <- rep(0,length(comp_names))  #one more for blood compartment
  Q <- rep(0,length(comp_names))
  
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Tissue_fractions[i] <- ifelse(is.na(control), NA, Tissue_fractions[i])
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
    Macrophage_fractions[i] <- ifelse(is.na(control), NA, Macrophage_fractions[i])
    
    ### Calculation of tissue weights  
    W_tis[i] <- w*Tissue_fractions[i]
    
    
    ###Calculation of tissue volumes
    
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
    
    ###Volume of macrophage contents
    W_macro[i] <- W_tis[i]*Macrophage_fractions[i]
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }
  
  #Vm_ven <- 0.01*Vven #macrophage content in veins
  #Vm_art <- 0.01*Vart #0.02*Vart #macrophage content in arteries
  
  ### Calculations for "Soft tissue" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)
  V_tis[1] <- W_tis[1]/d_adipose     
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
  V_cap[1] <- V_tis[1]*Capillary_fractions[1] #Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  W_macro[1] <- W_tis[1]*Macrophage_fractions[1]
  
  parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[],W_macro[]), ncol = 5)
  colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q", "W_macro")
  rownames(parameters) <- all_comps
  
  Vven=0.64*Total_Blood
  Vart=0.15*Total_Blood
  Wm_ven=0.01*Vven
  Wm_art=0.01*Vart
  
  return(c(
    "Q_total"=Q_total, "V_blood"=Total_Blood, "V_ven"=Vven, "V_art"=Vart,
    
    "w_rob"=parameters[1,1], "w_ht"=parameters[2,1], "w_ki"=parameters[3,1], "w_spl"=parameters[5,1], "w_lu"=parameters[6,1], "w_li"=parameters[7,1], "w_bone"=parameters[9,1], "w_git"=parameters[13,1],
    
    "V_tis_rob"=parameters[1,2], "V_tis_ht"=parameters[2,2], "V_tis_ki"=parameters[3,2], "V_tis_spl"=parameters[5,2], "V_tis_lu"=parameters[6,2], "V_tis_li"=parameters[7,2], "V_tis_bone"=parameters[9,2], "V_tis_git"=parameters[13,2], 
    
    "V_cap_rob"=parameters[1,3], "V_cap_ht"=parameters[2,3], "V_cap_ki"=parameters[3,3], "V_cap_spl"=parameters[5,3], "V_cap_lu"=parameters[6,3], "V_cap_li"=parameters[7,3], "V_cap_bone"=parameters[9,3], "V_cap_git"=parameters[13,3],
    
    "Q_rob"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_bone"=parameters[9,4], "Q_git"=parameters[13,4]
    
  ))
}

# Physiological parameters units
# V_blood, V_ven, V_art (ml): Volume of total blood, venous blood and arterial blood
# w_i (g):                    mass of tissue or organ "i"
# V_tis_i (ml):                volume of tissue or organ "i"
# V_cap_i (ml):                volume of capillary blood in tissue "i"
# Q_i, Q_total (ml/h):        regional blood flow of tissue or organ "i"


#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters, dose){
  with( as.list(parameters),{
    M_ht<-0; M_lu<-0; M_li<-0; M_spl<-0; M_ki<-0; M_git<-0; M_bone<-0; M_rob<-0;
    
    M_cap_ht<-0; M_cap_lu<-0; M_cap_li<-0; M_cap_spl<-0; M_cap_ki<-0; M_cap_git<-0; M_cap_bone<-0; M_cap_rob<-0;
    
    M_ven<-dose; M_art<-0
    M_feces<-0; M_urine<-0 
    
    return(c("M_ht" = M_ht, "M_lu" = M_lu, 
             "M_li" = M_li, "M_spl" = M_spl, 
             "M_ki" = M_ki, "M_git" = M_git, 
             "M_bone" = M_bone,"M_rob"=M_rob,
             
             "M_cap_ht" = M_cap_ht, "M_cap_lu" = M_cap_lu, 
             "M_cap_li" = M_cap_li, "M_cap_spl" = M_cap_spl, 
             "M_cap_ki" = M_cap_ki, "M_cap_git" = M_cap_git, 
             "M_cap_bone" = M_cap_bone,"M_cap_rob"=M_cap_rob,
             
             "M_ven" = M_ven, "M_art" = M_art, "M_feces" = M_feces, "M_urine" = M_urine))
    
  })
}

#==============
#3. ODEs System
#==============
ode.func <- function(time, inits, params){
  with(as.list(c(inits, params)),{
    
    
    x_ht <- x_gen
    x_lu <- x_gen
    x_li <- x_gen
    x_spl <- x_gen
    x_ki <- x_gen
    x_git <- x_gen 
    x_bone <- x_gen
    x_rob <- x_gen
    
    P_ht <- P_gen
    P_lu <- P_gen
    P_li <- P_gen
    P_spl <- P_gen
    P_ki <- P_gen
    P_git <- P_gen 
    P_bone <- P_gen
    P_rob <- P_gen
    
    # Concentrations (mg of NPs)/(g of wet tissue)
    C_ht <- M_ht/w_ht
    C_cap_ht <- M_cap_ht/V_cap_ht
    C_lu <- M_lu/w_lu
    C_cap_lu <- M_cap_lu/V_cap_lu
    C_li <- M_li/w_li
    C_cap_li <- M_cap_li/V_cap_li
    C_spl <- M_spl/w_spl
    C_cap_spl <- M_cap_spl/V_cap_spl
    C_ki <- M_ki/w_ki
    C_cap_ki <- M_cap_ki/V_cap_ki
    C_git <- M_git/w_git
    C_cap_git <- M_cap_git/V_cap_git
    C_bone <- M_bone/w_bone
    C_cap_bone <- M_cap_bone/V_cap_bone
    C_rob <- M_rob/w_rob
    C_cap_rob <- M_cap_rob/V_cap_rob
    
    C_ven <- M_ven/V_ven
    C_art <- M_art/V_art
    
    # Heart
    dM_cap_ht <- Q_ht*(C_art - C_cap_ht) - x_ht*Q_ht*(C_cap_ht - C_ht/P_ht)
    dM_ht <- x_ht*Q_ht*(C_cap_ht - C_ht/P_ht) 
    
    # Lungs
    dM_cap_lu <- Q_total*(C_ven - C_cap_lu) - x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
    dM_lu <-  x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
    
    # Liver 
    dM_cap_li <- Q_li*(C_art - C_cap_li) + Q_spl*(C_cap_spl - C_cap_li) - x_li*(Q_li)*(C_cap_li - C_li/P_li)
    dM_li <- x_li*Q_li*(C_cap_li - C_li/P_li) 
    
    # Spleen
    dM_cap_spl <- Q_spl*(C_art - C_cap_spl) - x_spl*Q_spl*(C_cap_spl - C_spl/P_spl)
    dM_spl <- x_spl*Q_spl*(C_cap_spl - C_spl/P_spl) 
    
    # Kidneys
    dM_cap_ki <- Q_ki*(C_art - C_cap_ki) - x_ki*Q_ki*(C_cap_ki - C_ki/P_ki)
    dM_ki <- x_ki*Q_ki*(C_cap_ki - C_ki/P_ki) - CLE_u*M_ki
    
    # GIT - Gastrointestinal Track
    dM_cap_git <- Q_git*(C_art - C_cap_git) - x_git*Q_git*(C_cap_git - C_git/P_git)
    dM_git <- x_git*Q_git*(C_cap_git - C_git/P_git) - CLE_f*M_git
    
    # Bone
    dM_cap_bone <- Q_bone*(C_art - C_cap_bone) - x_bone*Q_bone*(C_cap_bone - C_bone/P_bone)
    dM_bone <- x_bone*Q_bone*(C_cap_bone - C_bone/P_bone) 
    
    
    # RoB - Rest of Body
    dM_cap_rob <- Q_rob*(C_art - C_cap_rob) - x_rob*Q_rob*(C_cap_rob - C_rob/P_rob)
    dM_rob <- x_rob*Q_rob*(C_cap_rob - C_rob/P_rob) 
    
    # Urine
    dM_urine <- CLE_u*M_ki
    
    # Feces
    dM_feces <- CLE_f*M_git
    
    # Venous Blood
    dM_ven <- Q_ht*C_cap_ht + (Q_li + Q_spl)*C_cap_li + Q_ki*C_cap_ki + Q_git*C_cap_git +
      Q_bone*C_cap_bone + Q_rob*C_cap_rob - Q_total*C_ven
    
    # Arterial Blood
    dM_art <-  Q_total*C_cap_lu - Q_total*C_art
    
    Blood_total <- M_ven + M_art
    Blood <- Blood_total/(V_blood)
    
    list(c("dM_ht" = dM_ht, "dM_lu" = dM_lu, 
           "dM_li" = dM_li, "dM_spl" = dM_spl, 
           "dM_ki" = dM_ki, "dM_git" = dM_git, 
           "dM_bone" = dM_bone,"dM_rob"=dM_rob,
           
           "dM_cap_ht" = dM_cap_ht, "dM_cap_lu" = dM_cap_lu, 
           "dM_cap_li" = dM_cap_li, "dM_cap_spl" = dM_cap_spl, 
           "dM_cap_ki" = dM_cap_ki, "dM_cap_git" = dM_cap_git, 
           "dM_cap_bone" = dM_cap_bone,"dM_cap_rob"=dM_cap_rob,
           
           "dM_ven" = dM_ven, "dM_art" = dM_art, "dM_feces" = dM_feces, 
           "dM_urine" = dM_urine),
         
         "Blood"=Blood,
         "Heart"=C_ht, "Lungs"=C_lu, "Liver"=C_li, "Spleen"=C_spl,
         "Kidneys"=C_ki, "Git"=C_git, "Bone"=C_bone, "Rob"=C_rob,
         "Feces"=M_feces, "Urine"=M_urine)
  })
}


#---------------------
# Custom AUC function
#---------------------
AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
}

#---------------------------------------------------
# SENSITIVITY ANALYSIS FUNCTION for PBPK Models
#---------------------------------------------------

PBPK_sensitivity <- function(model, parms, ranges, targets, method,
                             ode_settings, heatmap = FALSE){
  
  if(is.null(model)){
    stop("The ODEs of the PBPK model must be provided as a function compatible
         with the deSolve library")
  }
  if(is.null(parms)){
    stop("A vector with the names and the values of the parameters for the 
         analysis must be provided")
  }
  if(is.null(method)){
    stop("The method of the sensitivity analysis must be defined
         (\"Local\" or \"Global\")")
  }
  
  if(method=="Local"){
    # Local Sensitivity Analysis
    
    # Store the initial values of the parms in a vector
    parms_0 <- parms
    
    # The variation dp of each parameter will be equal to "ranges" value 
    if(!is.numeric(ranges) | length(ranges) != 1 | ranges<=0 | ranges>1){
      # "ranges" must be a single numeric value in (0,1]
      stop("For local sensitivity analysis \"ranges\"  should be 
           a single numeric value in (0,1]")
    }else{dp <- ranges}
    
    # Get the number of the parameters and targets
    N_parms <- length(parms_0) # The total number of parameters to be analysed
    N_targets <- length(targets)
    
    # Take all objects from the ode_settings list 
    constant_params <- ode_settings[[1]] # Take the constant parameters
                                         # of the model
    inits <- ode_settings[[2]] # Initial conditions of the ODEs
    sample_time <- ode_settings[[3]] # Time points of solution
    solver <- ifelse(ode_settings[[4]] == "default", "bdf", ode_settings[[4]])
    # Select the solver to use later
    
    # Calculate the ODEs for the initial parameters
    solution_0 <- ode(times=sample_time, func=model, y=inits, 
                      parms=c(constant_params, parms_0), 
                      method=solver, rtol=1e-5, atol=1e-5)
    
    if(sum(!(targets %in% colnames(solution_0))) != 0){
      stop("As \"targets\" should be provided the name(s) of one or more
      of the outputs (compertments) of the PBPK model")
    }
    
    # Calculate AUC of each target-compartment for the initial parameters
    AUC_0 <- c()
    for (i in 1:N_targets) {
      AUC_0[i] <- AUC(solution_0[,"time"],solution_0[,targets[i]])
    }
    
    # Initialize a vector to store the sensitivity indexes
    SI <- matrix(NA, nrow = N_parms, ncol = N_targets)
    for (i in 1:N_parms) {
      parms <- parms_0
      parms[i] <- parms[i]*(1 + dp)
      # Merge the new parms vector with the constant params of the model
      params <- c(constant_params, parms)
      
      # Solve the ode system for given initial values and time
      solution <- ode(times=sample_time, func=model, y=inits, parms=params, 
                      method=solver, rtol=1e-5, atol=1e-5)
      
      for (j in 1:N_targets) {
        # Calculate AUC for the target compartment j
        AUC_j <- AUC(solution[,"time"],solution[,targets[j]])
        
        # Calculate sensitivity index of parameter i 
        # Relative Sensitivity = (dAUC/AUC)/(dp/p)
        SI[i,j] <- ((AUC_j-AUC_0[j])/AUC_0[j])/(dp/ parms_0[i])
      }
    }
    rownames(SI) <- names(parms)
    colnames(SI) <- targets
    
  }else{
    print("Global Sensitivity not ready!")
  }
  
  #----------
  # Heatmaps
  #----------
  if(heatmap){
    
    if(length(targets)==1){
      stop("Provide more than 1 targets in order to create a heatmap
           or turn \"heatmap\" to \"FALSE\".")
    }
    
    # Transform SI matrix to long format dataframe
    data_to_plot <- melt(SI)
    colnames(data_to_plot) <- c("Parameters","Targets", "SI")
    
    heatmap_plot <- ggplot(data_to_plot, aes(x = Targets, y = Parameters,
                                             fill = SI)) +
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1)+ # Add white borders to tiles
      
      geom_text(aes(label = format(round(SI, digits=3))),
                color = "white", size = 4) + # Show the SI value
      
      coord_fixed() + #Transforms the tiles of the heatmap to squares 
      
      guides(fill = guide_colourbar(barwidth = 0.5,
                                    barheight = 20)) +
      
      labs(title = "Sensitivity Coefficients",
           y = "Parameters",
           x = "Targets") +
      
      theme(plot.title =element_text(hjust = 0.5, size=30, face="bold"),
            axis.title.y =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.y=element_text(size=18),
            axis.title.x =element_text(hjust = 0.5, size=20, face="bold"),
            axis.text.x=element_text(size=18, 
                                     angle = ifelse(length(targets)>5, 90, 0)),
            legend.title=element_text(hjust = 0.01, size=20), 
            legend.text=element_text(size=18))
    Plots <- list(heatmap_plot)
  }
  
  #-------------------------------------
  # Create a list to return the results 
  #-------------------------------------
  
  data_list <- list("Method"=method, "Targets"=targets,
                    "Change of parameters" = dp,
                    "Parameters" = parms_0, 
                    "Normalized Sensitivity Coefficients"=data.frame(SI), 
                    "Heatmap" = ifelse(heatmap, Plots, "FALSE"))
  
  return(data_list)
  
  
}

#-----------------
# Prepare input
#-----------------

parms <- c("x_gen" = 0.1, # random value - unitless
           "P_gen" = 5, # random value - unitless
           "CLE_f" = 0.36,
           "CLE_u" = 0.0026)

params<-create.params(compartments,mass)
inits <- create.inits(params, dose)
sample_time <- c(0,1,3,7, 15, 30)*24 # hours

ode_settings <- list(params=params, 
                     inits=inits,
                     sample_time=sample_time,
                     solver= "default")


targets <- c("Lungs", "Blood", "Liver")
targets <- c("Blood", "Heart", "Lungs", "Liver", "Spleen", "Kidneys",
             "Git", "Bone", "Rob", "Feces", "Urine")



test <- PBPK_sensitivity(model = ode.func, parms=parms, ranges = 0.1,
                         targets = targets,
                         method = "Local",ode_settings = ode_settings, 
                         heatmap = TRUE) 

test



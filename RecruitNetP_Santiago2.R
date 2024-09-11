# download_RN

download_RN <- function(path = getwd(), destfile = "RN.zip", unzip = TRUE) {
  
  utils::download.file("https://zenodo.org/record/6567608/files/Data_S1.zip?download=1",
                       destfile = file.path(path, destfile), mode = "wb")
  
  if (isTRUE(unzip)) {
    utils::unzip(file.path(path, destfile), exdir = path)
    file.remove(file.path(path, destfile))
  }
  
}

# -----------------------------------------

# comm_subset

comm_subset <- function(int_data, site = NULL) {
  
  if(length(unique(int_data$Study_site))>1 & is.null(site)==TRUE) stop("You must enter a site name.")
  # Formatting
  ifelse(is.null(site)==TRUE, 
         df <- int_data, 
         df <- int_data[int_data$Study_site %in% site, ])
  return(df)
}

#-----------------------------------------------

# comm_summary

comm_summary <- function(int_data = NULL) {
  
  dfList <- list()
  
  for (i in 1:length(unique(int_data$Study_site))) {
    data <- comm_subset(int_data, unique(int_data$Study_site)[i])
    
    # Plots data
    n_plots <- length(unique(data$Plot))
    plotDim <- data$PlotdimX[1] * data$PlotdimY[1]
    areaSampled <- plotDim * n_plots
    
    # Nodes data
    isOpen <- ifelse(is.element('Open', data$Standardized_Canopy) == TRUE, "Yes", "No")
    n_nodes <- length(unique(c(
      data$Standardized_Canopy, data$Standardized_Recruit
    )))
    n_sp <- ifelse(isOpen == "Yes", n_nodes - 1, n_nodes)
    num_woody <- length(which(unique(data.frame(
      c(data$LifeHabit_Canopy, data$LifeHabit_Recruit),
      c(data$Standardized_Canopy, data$Standardized_Recruit)
    ), margin = 1) == "W"))
    num_herbs <- length(which(unique(data.frame(
      c(data$LifeHabit_Canopy, data$LifeHabit_Recruit),
      c(data$Standardized_Canopy, data$Standardized_Recruit)
    ), margin = 1) == "H"))
    num_others <- n_sp - num_woody - num_herbs
    
    # Function output
    df <- data.frame(
      c(
        data$Study_site[1],
        data$Country[1],
        data$Latitude[1],
        data$Longitude[1],
        data$Sampling_date[1],
        data$Site_responsible[1],
        data$Biome[1],
        data$Vegetation_type[1],
        data$Community[1],
        data$Successional_stage[1],
        data$Disturbance[1],
        data$Sampling_method[1],
        n_plots,
        plotDim,
        areaSampled,
        n_sp,
        isOpen,
        num_woody,
        num_herbs,
        num_others
      )
    )
    colnames(df) <- c("Value")
    rownames(df) <- c(
      "Local Community",
      "Country",
      "Latitude",
      "Longitude",
      "Year of sampling",
      "Site responsible",
      "Biome",
      "Vegetation",
      "Plant Community",
      "Successional stage",
      "Disturbance",
      "Sampling method",
      "Number of plots",
      "Plot area (m2)",
      "Area sampled (m2)",
      "Number of plant species",
      "Contains Open node",
      "Number of woody species",
      "Number of herb species",
      "Number of other growth habits"
    )
    
    dfList[[i]] <-  t(df)
    
  }
  
  dfAll <-  as.data.frame(do.call(rbind, dfList))
  
  dfAll <- utils::type.convert(dfAll, as.is = TRUE)
  
  return(dfAll)
  
}

# ---------------------------------------------------

# remove_no_cover

remove_no_cover <- function(int_data=NULL, cover_data=NULL) {
  
  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_data$Canopy))
  RN_list <- sort(unique(c(int_data$Canopy, int_data$Recruit)))
  lack_cover <- setdiff(RN_list, cover_list)
  
  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    df <- int_data
  } else {
    df <- int_data[-which(int_data$Recruit %in% lack_cover), ]
  }
  
  if (length(lack_cover) == 0) {
    df
  } else {
    if(length(which(df$Canopy %in% lack_cover))==0) {
      df
    } else {
      df <- df[-which(df$Canopy %in% lack_cover), ]
    }
  }
  return(df)
}

# -------------------------------------------------------

# aggr_RN

aggr_RN <- function(int_data) {
  
  # Sum the number of recruits per interaction across plots
  RN <- aggregate(Frequency ~ Canopy*Recruit, data = int_data, FUN = sum)
  colnames(RN) <- c("Canopy", "Recruit", "Fcr")
  RN$Icr <- aggregate(Frequency ~ Canopy*Recruit, data=int_data, FUN = NROW)[[3]]
  RN$Pcr <- ifelse(RN$Icr==0,0,1)
  RN$Canopy <- gsub("[[:space:]]", "_", RN$Canopy)
  RN$Recruit <- gsub("[[:space:]]", "_", RN$Recruit)
  
  # Incorporate the unobserved interactions
  species_list <- unique(c(RN$Canopy, RN$Recruit))
  df <- expand.grid(Canopy = species_list, Recruit = species_list)
  RN <- merge(df, RN, all = TRUE)
  RN[is.na(RN)] <- 0
  RN$Canopy <- as.character(RN$Canopy)
  RN$Recruit <- as.character(RN$Recruit)
  # RN <- RN[which(RN$Recruit!="Open"),] # Remove Open from the Recruit species
  return(RN)
  
}

# -------------------------------------------------------

# aggr_cover

aggr_cover <- function(cover_data) {
  # Calculations
  cover_data$Ac <- (cover_data$Cover / 100) * cover_data$Sampled_distance_or_area
  cover <- stats::aggregate(Ac ~ Canopy, data = cover_data, sum)
  return(cover)
}

# -------------------------------------------------------

# comm_to_RN

comm_to_RN <- function(int_data, cover_data) {
  
  # Aggregate data sets across plots.
  int_df <- aggr_RN(int_data)
  cover_df <- aggr_cover(cover_data)
  
  # Find species present in RN but that lack data on cover.
  cover_list <- sort(unique(cover_df$Canopy))
  RN_list <- sort(unique(c(int_df$Canopy, int_df$Recruit)))
  lack_cover <- setdiff(RN_list, cover_list)
  
  # Remove species lacking cover from RN
  if (length(lack_cover) == 0) {
    RNc <- int_df
  } else {
    RNc <- int_df[-which(int_df$Recruit %in% lack_cover), ]
  }
  
  if (length(lack_cover) == 0) {
    RNc
  } else {
    if(length(which(RNc$Canopy %in% lack_cover))==0) {
      RNc
    } else {
      RNc <- RNc[-which(RNc$Canopy %in% lack_cover), ]
    }
  }
  
  # Add variables with the cover of the canopy (Ac) and recruit (Ar) species
  RNc$Ac <- RNc$Canopy
  RNc$Ar <- RNc$Recruit
  for (i in 1:dim(RNc)[1]) {
    RNc$Ac[i] <- as.numeric(replace(
      RNc$Ac[i],
      match(RN_list, RNc$Canopy[i]),
      cover_df$Ac[match(RNc$Canopy[i], cover_df$Canopy)]
    ))
  }
  
  for (i in 1:dim(RNc)[1]) {
    RNc$Ar[i] <- as.numeric(replace(
      RNc$Ar[i],
      match(RN_list, RNc$Recruit[i]),
      cover_df$Ac[match(RNc$Recruit[i], cover_df$Canopy)]
    ))
  }
  
  RNc <- utils::type.convert(RNc, as.is = TRUE)
  RNc <- RNc[order(RNc$Canopy, RNc$Recruit),]
  
  return(RNc)
  
}

# -------------------------------------------------------

# RN_to_matrix()

RN_to_matrix <- function(int_data=NULL, weight = "Fcr"){
  
  # Check column names
  if ("Canopy" %in% names(int_data) == FALSE) warning("ERROR: your data lacks a column named: Canopy")
  if ("Recruit" %in% names(int_data) == FALSE) warning("ERROR: your data lacks a column named: Recruit")
  
  data <- int_data
  # Formatting 
  list_Canopy <- sort(unique(data$Canopy))
  list_Recruit <- sort(unique(data$Recruit))
  Num_canopy <- length(list_Canopy)
  Num_recruit <- length(list_Recruit)
  RNmat <- data[[weight]]
  dim(RNmat) <- c(Num_recruit, Num_canopy)
  colnames(RNmat) <- list_Canopy
  rownames(RNmat) <- list_Recruit
  return(RNmat)
  
}

# -------------------------------------------------------

# link_completeness()

link_completeness <- function(int_data = NULL,
                              type = c("incidence", "abundance")) {
  
  stopifnot(
    c("Plot",
      "Canopy",
      "Recruit",
      "Frequency"
    ) %in% names(int_data))
  
  type <- match.arg(type)
  
  data_raw <- int_data
  data_RN <- aggr_RN(data_raw)
  
  # Completeness based on incidence data.
  
  if (type == "incidence") {
    netRaw <- data.frame(cbind(data_raw$Plot, paste(data_raw$Canopy, data_raw$Recruit)))
    colnames(netRaw) <- c("Plot", "Pair")
    nPlots <- length(unique(netRaw$Plot))
    
    # Check points.
    
    if (nPlots == 1)
      stop(
        "ERROR: your data is not structured in multiple plots. Incidence approach cannot be used. Try the abundance approach."
      )
    if (nPlots < 10)
      warning(
        "WARNING: your are using the incidence approach with very few plots. Consider using the abundance approach if appropriate."
      )
    
    # Combine the lists of canopy and recruit species to obtain the total list of canopy-recruit pairs (links) sampled.
    
    a1 <- split(netRaw[-1], f = netRaw[1])
    a2 <- lapply(a1, unique)
    a3 <- unlist(unlist(a2, recursive = FALSE, use.names = FALSE))
    
    # Table showing the incidence of each canopy-recruit pair in the study site
    
    a4 <- table(a3)
    linkIncidence <- as.data.frame(a4)
    colnames(linkIncidence) <- c("Pair", "Incidence")
    
    # Incidence list to be passed to iNEXT
    
    data_iNEXT <- c(nPlots, sort(linkIncidence$Incidence, decreasing = TRUE))
    
    # Call to iNEXT to obtain completeness values
    
    out <- iNEXT::iNEXT(
      data_iNEXT,
      q = c(0, 1),
      datatype = "incidence_freq",
      se = FALSE,
      size = nPlots
    )
    Lobs <- out$AsyEst[1, 1]
    Lest <- out$AsyEst[1, 2]
    Lest_LCL <- out$AsyEst[1, 4]
    Lest_UCL <- out$AsyEst[1, 5]
    Cq0_L <- Lobs / Lest
    Cq1_L <- out$DataInfo[1, 5]
    df <- data.frame(c(Lobs, Lest, Cq0_L, Cq1_L))
    colnames(df) <- c("Incidence based estimate")
    rownames(df) <- c("Lobs",
                      "Lest",
                      "Completeness Links (q=0)",
                      "Coverage Links (q=1)")
  }
  
  # Completeness based on abundance or frequency of recruits.
  
  if (type == "abundance") {
    # Call to iNEXT to obtain completeness values
    
    out <- iNEXT::iNEXT(data_RN$Fcr[which(data_RN$Fcr > 0)], q = 0, datatype = "abundance")
    Lobs <- out$AsyEst[1, 1]
    Lest <- out$AsyEst[1, 2]
    Cq0_L <- Lobs / Lest
    Cq1_L <- out$DataInfo[1, 4]
    df <- data.frame(c(Lobs, Lest, Cq0_L, Cq1_L))
    colnames(df) <- c("Abundance based estimate")
    rownames(df) <- c("Lobs",
                      "Lest",
                      "Completeness of links (q=0)",
                      "Coverage of links (q=1)")
  }
  
  return(df)
}

# -------------------------------------------------------

# pre_associndex()

pre_associndex <- function(int_data = NULL) {
  # Formatting
  mydata <- int_data[,c("Canopy", "Recruit", "Fcr","Ac")]
  
  # Frequency in open
  Fcr_matrix_df <- as.data.frame(RN_to_matrix(mydata))
  Open_matrix <- rep(Fcr_matrix_df$Open, length(Fcr_matrix_df))
  dim(Open_matrix) <- c(length(Fcr_matrix_df), length(Fcr_matrix_df))
  mydata$Fro <- unlist(as.vector(Open_matrix)) 
  
  # Cover of Open
  cov_matrix_df <- as.data.frame(RN_to_matrix(mydata, "Ac"))
  Open_cov_matrix <- rep(cov_matrix_df$Open, length(cov_matrix_df))
  dim(Open_cov_matrix) <- c(length(cov_matrix_df), length(cov_matrix_df))
  mydata$Ao <- unlist(as.vector(Open_cov_matrix)) 
  
  # Remove absent interactions
  mydata <- mydata[mydata$Fcr+mydata$Fro > 0,]
  # Remove "Open" from the list of canopy species
  mydata <- mydata[-which(mydata$Canopy=="Open"),]
  return(mydata)
  
}

# -------------------------------------------------------

# associndex()

associndex <- function(int_data = NULL,
                               threshold_density = 100) {
  
  thr <- threshold_density
  
  # Assemble the data
  db_inter <- pre_associndex(int_data)
  
  # Incorporate density of recruitment (recruits/m2) under each canopy species and in open.
  db_inter$Dcr <- db_inter$Fcr/db_inter$Ac
  db_inter$Dro <- db_inter$Fro/db_inter$Ao
  
  # Retain the interactions with estimated density below the threshold.
  db_inter <- db_inter[which(db_inter$Dcr<thr & db_inter$Dro<thr), ]
  
  #Obtain the maximum recruitment density for each recruit under the canopy species or in open.
  db_inter$Max_Recr_Density <- pmax(db_inter$Dcr,db_inter$Dro)
  
  db_inter <- utils::type.convert(db_inter, as.is = TRUE)
  
  max_rd <- stats::aggregate(Max_Recr_Density ~ Recruit, data = db_inter, FUN = "max")
  
  # Add a variable max_Recr to each pair indicating the maximum recruitment density of the recruit species in the study site
  Recr_list <- sort(unique(c(db_inter$Recruit)))
  Dens_list <- sort(unique(max_rd$Recruit))
  lack_dens <- setdiff(Recr_list, Dens_list)
  
  db_inter$max_Recr <- db_inter$Recruit
  for (i in 1:(dim(db_inter)[1])) {
    db_inter$max_Recr[i] <- replace(
      db_inter$max_Recr[i],
      match(Recr_list, db_inter$max_Recr[i]),
      max_rd$Max_Recr_Density[match(db_inter$max_Recr[i], max_rd$Recruit)]
    )
  }
  
  db_inter <- utils::type.convert(db_inter, as.is = TRUE)
  
  # Calculate indices Ns, NintC, NintA and RII
  db_inter$Ns <- (db_inter$Dcr - db_inter$Dro)/db_inter$max_Recr
  db_inter$NintC <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dcr + db_inter$Dro)+abs(db_inter$Dcr-db_inter$Dro))
  db_inter$NintA <- 2*(db_inter$Dcr - db_inter$Dro)/((db_inter$Dro) + abs(db_inter$Dcr-db_inter$Dro))
  db_inter$RII <- (db_inter$Dcr - db_inter$Dro)/(db_inter$Dcr + db_inter$Dro)
  
  removed <- names(db_inter) %in% c("Frequency", "Max_Recr_Density")
  db_inter <- db_inter[!removed]
  return(db_inter)
  
}

# -------------------------------------------------------

# int_significance()

int_significance <- function(int_data){
  
  df <- pre_associndex(int_data)
  n_tests <- dim(df)[1]
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)
  
  # Testability through Binomial test
  
  df$Ftot <- df$Fcr+df$Fro
  
  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p
  
  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability
  
  # Binomial (or Chi square) Test Significance
  
  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- binom.test(df$Fcr[i], df$Fcr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- chisq.test(c(df$Fcr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance
  
  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fcr[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fcr[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type
  
  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")
  
  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fcr[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Enhancing",
                         Effect_int[i] <- "Depressing")
           )
    )
  }
  
  df$Effect_int <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]
  return(df)
}

# -------------------------------------------------------

# canopy_effect_on_rec()

canopy_effect_on_rec <- function(int_data){
  
  df <- pre_associndex(int_data)
  sp_Fc <- stats::aggregate(Fcr ~ Canopy, data = df, FUN = sum)
  sp_Ac <- stats::aggregate(Ac ~ Canopy, data = df, FUN = max)
  sp_Fro <- stats::aggregate(Fro ~ Canopy, data = df, FUN = sum)
  sp_Ao <- stats::aggregate(Ao ~ Canopy, data = df, FUN = max)
  n_tests <- dim(sp_Fc)[1]
  df <- data.frame(c(sp_Fc, sp_Ac, sp_Fro, sp_Ao))
  myvars <- names(df) %in% c("Canopy.1", "Canopy.2", "Canopy.3")
  df <- df[!myvars]
  colnames(df) <- c("Canopy", "Fc", "Ac", "Fro", "Ao")
  df$exp_p <- df$Ac/(df$Ac+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)
  
  # Testability through Binomial test
  
  df$Ftot <- df$Fc+df$Fro
  
  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p
  
  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability
  
  # Binomial (or Chi square) Test Significance
  
  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fc[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fc[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Significance[i] <- binom.test(df$Fc[i], df$Fc[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- chisq.test(c(df$Fc[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance
  
  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fc[i]+df$Fro[i])*(df$Ac[i]/(df$Ac[i]+df$Ao[i]))<=5 | (df$Fc[i]+df$Fro[i])*(df$Ao[i]/(df$Ac[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type
  
  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")
  
  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fc[i]/df$Ac[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Facilitative",
                         Effect_int[i] <- "Depressive")
           )
    )
  }
  
  df$Canopy_effect <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]
  return(df)
}

# -------------------------------------------------------

# veg_effect_on_rec()

veg_effect_on_rec <- function(int_data){
  
  df <- pre_associndex(int_data)
  sp_Fr <- stats::aggregate(Fcr ~ Recruit, data = df, FUN = sum)
  sp_Av <- stats::aggregate(Ac ~ Recruit, data = df, FUN = sum)
  sp_Fro <- stats::aggregate(Fro ~ Recruit, data = df, FUN = max)
  sp_Ao <- stats::aggregate(Ao ~ Recruit, data = df, FUN = max)
  n_tests <- dim(sp_Fr)[1]
  df <- data.frame(c(sp_Fr, sp_Av, sp_Fro, sp_Ao))
  myvars <- names(df) %in% c("Recruit.1", "Recruit.2", "Recruit.3")
  df <- df[!myvars]
  colnames(df) <- c("Recruit", "Fr", "Av", "Fro", "Ao")
  df$exp_p <- df$Av/(df$Av+df$Ao) # Expected probability of success (i.e. of recruiting under canopy)
  
  # Testability through Binomial test
  
  df$Ftot <- df$Fr+df$Fro
  
  extreme_p <- c()
  for(i in 1:n_tests){
    extreme_p[i] <- min(df$exp_p[i], 1-df$exp_p[i])
  }
  df$extreme_p <- extreme_p
  
  testability <- c()
  for(i in 1:n_tests) {
    testability[i] <- binom.test(df$Ftot[i], df$Ftot[i], df$extreme_p[i], alternative ="two.sided")$p.value
  }
  df$testability <- testability
  
  # Binomial (or Chi square) Test Significance
  
  Significance <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fr[i]+df$Fro[i])*(df$Av[i]/(df$Av[i]+df$Ao[i]))<=5 | (df$Fr[i]+df$Fro[i])*(df$Ao[i]/(df$Av[i]+df$Ao[i]))<=5),
           Significance[i] <- binom.test(df$Fr[i], df$Fr[i]+df$Fro[i], df$exp_p[i], alternative ="two.sided")$p.value,
           Significance[i] <- chisq.test(c(df$Fr[i], df$Fro[i]), p = c(df$exp_p[i], 1-df$exp_p[i]))$p.value
    )
  }
  df$Significance <- Significance
  
  Test_type <- c()
  for(i in 1:n_tests) {
    ifelse(((df$Fr[i]+df$Fro[i])*(df$Av[i]/(df$Av[i]+df$Ao[i]))<=5 | (df$Fr[i]+df$Fro[i])*(df$Ao[i]/(df$Av[i]+df$Ao[i]))<=5),
           Test_type[i] <- "Binomial",
           Test_type[i] <- "Chi-square"
    )
  }
  df$Test_type <- Test_type
  
  if(length(unique(df$Test_type))>1) warning("Different tests were used for different canopy-recruit pairs. Check column Test_type")
  
  Effect_int <- c()
  for(i in 1:n_tests) {
    ifelse((df$testability[i]>0.05),
           Effect_int[i] <- "Not testable",
           ifelse(df$Significance[i] > 0.05,
                  Effect_int[i] <- "Neutral",
                  ifelse((df$Fr[i]/df$Av[i])>(df$Fro[i]/df$Ao[i]),
                         Effect_int[i] <- "Facilitated",
                         Effect_int[i] <- "Depressed")
           )
    )
  }
  
  df$Veg_effect <- Effect_int
  drops <- c("exp_p", "Ftot", "extreme_p")
  df <- df[ , !(names(df) %in% drops)]
  return(df)
}

# -------------------------------------------------------

# node_degrees()

node_degrees <- function(int_data) {
  
  RNc <- int_data
  
  matrix_Fcr <- RN_to_matrix(RNc, weight = "Fcr")
  matrix_Pcr <- RN_to_matrix(RNc, weight = "Pcr")
  matrix_Ac <- RN_to_matrix(RNc, weight = "Ac")
  
  # Degree properties
  canopy_service_width <- colSums(matrix_Pcr)
  canopy_contribution <- colSums(matrix_Fcr)
  recruit_niche_width <- rowSums(matrix_Pcr)
  recruit_bank_abundance <- rowSums(matrix_Fcr)
  node_abund <- matrix_Ac[1,]
  
  # Specialization
  effective_canopy_service <- exp(-1*rowSums((t(matrix_Fcr)/canopy_contribution)*log(t(matrix_Fcr)/canopy_contribution), na.rm=TRUE))
  effective_recruit_niche <- exp(-1*rowSums((matrix_Fcr/recruit_bank_abundance)*log(matrix_Fcr/recruit_bank_abundance), na.rm=TRUE))
  
  node_deg <- data.frame(cbind(names(canopy_service_width), node_abund, canopy_service_width, canopy_contribution, effective_canopy_service, recruit_niche_width, recruit_bank_abundance, effective_recruit_niche))
  colnames(node_deg) <- c("Node", "Ac", "canopy_service_width", "canopy_contribution", "effective_canopy_service", "recruitment_niche_width", "recruit_bank_abundance", "effective_recruitment_niche")
  node_deg <- utils::type.convert(node_deg, as.is = TRUE)
  rownames(node_deg) <- NULL
  return(node_deg)
}

# -------------------------------------------------------

# partial_RNs()

partial_RNs <- function(int_data, k) {
  
  if (!"Plot" %in% names(int_data)) stop("ERROR: your interactions data lacks a column named Plots. This function requires data assembled in plots.")
  
  require(igraph)
  # Prepare the data
  int_data$Pcr <- ifelse(int_data$Frequency==0,0,1)
  netRaw <- data.frame(cbind(int_data$Plot, int_data$Canopy, int_data$Recruit, int_data$Pcr))
  colnames(netRaw) <- c("Plot", "Canopy", "Recruit", "Pcr")
  netRaw <- transform(netRaw, Pcr = as.numeric(Pcr))
  nPlots <- length(unique(netRaw$Plot))
  
  # Make the network of each plot
  plot_RNs <- c()
  for(i in 1:nPlots) {
    plot_RNs[[i]] <- igraph::graph_from_data_frame(netRaw[netRaw$Plot == i, 2:4], directed = TRUE)
  }
  
  # Combine the networks of n plots k times.
  union_RNs <- list()
  for(i in 1:nPlots) {
    union_RNs[[i]] <- replicate(k, do.call(igraph::union, sample(plot_RNs, i)))
  }
  
  # Add the incidence of each interaction across plots as an edge weight property to each partial network.
  for(i in 1:nPlots){
    for(j in 1:k){
      Icr <- rowSums(do.call(cbind.data.frame, igraph::edge.attributes(union_RNs[[i]][[j]])), na.rm=TRUE)
      union_RNs[[i]][[j]]<-igraph::set_edge_attr(union_RNs[[i]][[j]], "Icr", value = Icr)
    }}
  
  return(union_RNs)
}

# -------------------------------------------------------

# cum_values()

cum_values <- function(int_data, property, k = 100, title){
  require(ggplot2)
  part_RNs <- partial_RNs(int_data, k)
  nSteps <- length(part_RNs) 
  borrar <- unlist(part_RNs, recursive = FALSE)
  df <- data.frame(unlist(lapply(borrar, property)))
  colnames(df) <- c("Value")
  df$sampleSize <- sort(rep(c(1:nSteps),k))
  plot_cumm_value <- ggplot(df, aes(x=as.factor(sampleSize), y=Value)) +
    geom_jitter(colour="turquoise3", alpha=0.5, height = 0, width=0.1) +
    geom_point(stat="summary", fun="mean") +
    geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width=0.3) +
    labs(x="Sample Size (Num. Plots)", y="Value (mean + 95%CI)") +
    ggtitle(title)
  return(plot_cumm_value)
}

# -------------------------------------------------------

# RN_dims()

RN_dims <- function(int_data){
  df <- int_data
  n_nodes <- length(unique(c(df$Canopy, df$Recruit)))
  n_links <- sum(df$Pcr)
  connectance <- n_links/(n_nodes^2 - n_nodes)
  
  out <- data.frame(c(n_nodes, n_links, connectance))
  colnames(out) <- c("Value")
  rownames(out) <- c("Num. Nodes", "Num. Links", "Connectance")
  return(out)
}

# -------------------------------------------------------

# node_topol()

node_topol <- function(int_data) {
  RN_igraph <- igraph::graph_from_adjacency_matrix(t(RN_to_matrix(int_data, weight = "Pcr")), mode = "directed")
  eigen_cent <- igraph::eigen_centrality(RN_igraph, directed=TRUE, scale=FALSE, options = list(which="LR"))$vector
  out_neigh <- igraph::neighborhood_size(RN_igraph, order=gorder(RN_igraph), mode="out", mindist=1)
  in_neigh <- igraph::neighborhood_size(RN_igraph, order=gorder(RN_igraph), mode="in", mindist=1)
  df <- data.frame(eigen_cent, out_neigh,in_neigh)
  df[, 1] <- round(df[, 1], digits = 4)
  colnames(df) <- c("Eigenvector centrality", "Extended canopy service", "Extended recruitment niche")
  return(df)
}

# -------------------------------------------------------

# funtopol()

funtopol <- function(int_data){
  
  if (!"Open" %in% int_data$Canopy) stop("ERROR: your data does not contain a node named Open or it is spelled differently.")
  
  int_data <- int_data[which(int_data$Fcr!=0), c("Canopy", "Recruit")]
  g <- igraph::graph_from_data_frame(int_data, directed = TRUE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = FALSE)
  NEdges <- igraph::gsize(g)
  NNodes <- igraph::gorder(g)
  CDirected <- NEdges/(NNodes*(NNodes - 1))
  SCCs <- igraph::components(g, mode = "strong")
  numSCCs <- SCCs$no
  numNTSCCs <- sum(SCCs$csize > 1)
  coreSize <- max(SCCs$csize)
  SCC_memb <- SCCs$membership
  SCC_memb <- as.data.frame(SCC_memb)
  SCC_subgraphs <- igraph::decompose(g, mode = "strong") # Makes a subgraph of each SCC
  IDcore <- match(coreSize, SCCs$csize) # locates the position of the core in the list of SCCs
  MembersCore <- igraph::V(SCC_subgraphs[[IDcore]])$name # List of the species in the core
  IDOpen <- SCC_memb$SCC_memb[match("Open", row.names(SCC_memb))] # Locate the position of the "open" node in the list of SCCs
  outReachFromOpen <- names(igraph::subcomponent(g, "Open", "out")[-1]) # List of nodes reachable from the "open"
  outReachFromCore <- vector("list", coreSize) # List of nodes reachable from core nodes
  for (i in 1:coreSize) {
    outReachFromCore[[i]] <- igraph::subcomponent(g, MembersCore[i], mode = "out")
  }
  a <- unlist(outReachFromCore)
  a <- unique(names(a))
  MembersSatellites <- setdiff(a, MembersCore)
  MembersTransients <- setdiff(igraph::V(g)$name,c(MembersCore,MembersSatellites))
  MembersTransients <- MembersTransients[!MembersTransients == "Open"]
  MembersStrictTransients <- setdiff(MembersTransients, outReachFromOpen)
  MembersDdTransients <- setdiff(MembersTransients, MembersStrictTransients)
  numSat <- length(MembersSatellites)
  numTransAll <- length(MembersTransients)
  numDdTrans <- length(MembersDdTransients)
  numStrictTrans <- length(MembersStrictTransients)
  propCore <- coreSize/(NNodes - 1)
  propSat <- numSat/(NNodes - 1)
  propTrans <- numTransAll/(NNodes - 1)
  propStrTrans <- numStrictTrans/(NNodes - 1)
  propDdTrans <- numDdTrans/(NNodes - 1)
  persistence <- propCore + propSat
  
  # Function output
  
  df <- data.frame(
    c(NNodes,
      NEdges,
      CDirected,
      numNTSCCs,
      coreSize,
      propCore,
      numSat,
      propSat,
      numDdTrans,
      propDdTrans,
      numStrictTrans,
      propStrTrans,
      persistence)
  )
  colnames(df) <- c("Value")
  rownames(df) <- c(
    "Num. nodes",
    "Num. edges",
    "Connectance",
    "Num. non-trivial SCCs",
    "Num. core species",
    "Prop. core species",
    "Num. satellite species",
    "Prop. satellite species",
    "Num. disturbance-dependent transients",
    "Prop. disturbance-dependent transients",
    "Num. strict transients",
    "Prop. strict transients",
    "Qualitative Persistence")
  
  classif <- list(
    MembersSatellites,
    MembersCore,
    MembersStrictTransients,
    MembersDdTransients)
  classif <- stats::setNames(classif,
                             c("Satellites",
                               "Core",
                               "Strict_transients",
                               "Disturbance_dependent_transients")
  )
  
  df0 <- classif
  df_Sat <- data.frame(df0$Satellites, rep("Satellite", length(df0$Satellites)))
  colnames(df_Sat) <- c("id", "group")
  df_Core <- data.frame(df0$Core, rep("Core", length(df0$Core)))
  colnames(df_Core) <- c("id", "group")
  df_Str <- data.frame(df0$Strict_transients, rep("Strict_transients", length(df0$Strict_transients)))
  colnames(df_Str) <- c("id", "group")
  df_Ddtr <- data.frame(df0$Disturbance_dependent_transients, rep("Disturbance_dependent_transients", length(df0$Disturbance_dependent_transients)))
  colnames(df_Ddtr) <- c("id", "group")
  df0 <- rbind(df_Sat, df_Core, df_Str, df_Ddtr)
  df0 <- df0[order(df0$id),]
  
  outputs <- list("Descriptors" = df, "Functional_classification" = df0)
  
  return(outputs)
  
}

# -------------------------------------------------------

# visu_funtopol()

visu_funtopol <- function(int_data){
  
  if (!"Open" %in% int_data$Canopy) stop("ERROR: your data does not contain a node named Open or it is spelled differently.")
  
  nodes_list <- funtopol(int_data)$Functional_classification
  open_df <- c("Open", "Open")
  nodes_list <- rbind(nodes_list, open_df)
  nodes_list$label <- nodes_list$id
  int_data <- int_data[which(int_data$Fcr!=0), c("Canopy", "Recruit")]
  g <- igraph::graph_from_data_frame(int_data, directed = TRUE)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = FALSE)
  g<- igraph::as_data_frame(g, what = "both")
  edges_list <- g$edges
  
  # nodes data.frame for legend
  lnodes <- data.frame(label = c("Open", "Core", "Satellite", "Strict transient", "Disturbance-dependent transient"),
                       shape = c( "dot"), color = c("#F0E442", "#009E73", "#0072B2", "#D55E00", "#CC79A7"),
                       title = "Functional types", id = 1:5)
  
  # Network visualization and export to html
  require(visNetwork)
  network <- visNetwork(nodes_list, edges_list) %>%
    visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
    visEdges(arrows ="to") %>%                          
    visGroups(groupname = "Open", color = "#F0E442") %>%  
    visGroups(groupname = "Core", color = "#009E73") %>%
    visGroups(groupname = "Satellite", color = "#0072B2") %>%
    visGroups(groupname = "Sctrict_transient", color = "#D55E00") %>%
    visGroups(groupname = "Disturbance_dependent_transients", color = "#CC79A7") %>%
    visOptions(nodesIdSelection = TRUE) %>%
    visNetwork::visLegend(addNodes = lnodes, useGroups = FALSE)
  
  visSave(network, file = "network.html") # Save the html version of the network
  
  return(network)
}

# -------------------------------------------------------

# RN_heatmap

RN_heatmap <- function(int_data, weight_var = c("Fcr", "Dcr", "Icr", "Pcr"), scale_top = 1) {
  require(ggplot2)
  # manually set node order
  canopy_order <- unique(int_data$Canopy)
  canopy_order <- canopy_order[!canopy_order %in% c('Open')]
  canopy_order <- c("Open", canopy_order)
  int_data$Canopy2 <- factor(int_data$Canopy, levels = canopy_order)
  recruit_order <- sort(unique(int_data$Canopy), decreasing = TRUE)
  recruit_order <- recruit_order[!recruit_order %in% c('Open')]
  recruit_order <- c(recruit_order, "Open")
  int_data$Recruit2 <- factor(int_data$Recruit, levels = recruit_order)
  
  # Add recruitment density as another weighting variable
  int_data$Dcr <- int_data$Fcr/int_data$Ac
  
  # Make weight variable
  int_data$weight <- int_data[weight_var]
  
  # Lowest (non-zero) and highest values of the weighting variable
  highest_W <- max(int_data$weight)
  lowest_W <- min(int_data$weight[int_data$weight>0])
  
  # Plot the heatmap
  ggplot(int_data, aes(Canopy2, Recruit2, fill= Dcr)) + 
    geom_tile(colour="gray", size=0.25, aes(height = 1)) +
    scale_fill_gradientn(colours = c("#F5F5F5", "#E69F00","#0072B2"), values = c(0,lowest_W, scale_top*highest_W)) +
    #    scale_fill_gradient(low="white", high="turquoise3")+
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
}

# -------------------------------------------------------

# RN_heatmap

RN_heatmap_non_sq <- function(int_data, weight_var = c("Fcr", "Dcr", "Icr", "Pcr"), scale_top = 1) {
  require(ggplot2)
  # manually set node order
  canopy_order <- unique(int_data$Canopy)
  int_data$Canopy2 <- factor(int_data$Canopy, levels = canopy_order)
  recruit_order <- sort(unique(int_data$Recruit), decreasing = TRUE)
    int_data$Recruit2 <- factor(int_data$Recruit, levels = recruit_order)
  
  # Add recruitment density as another weighting variable
  int_data$Dcr <- int_data$Fcr/int_data$Ac
  
  # Make weight variable
  int_data$weight <- int_data[weight_var]
  
  # Lowest (non-zero) and highest values of the weighting variable
  highest_W <- max(int_data$weight)
  lowest_W <- min(int_data$weight[int_data$weight>0])
  
  # Plot the heatmap
  ggplot(int_data, aes(Canopy2, Recruit2, fill= Dcr)) + 
    geom_tile(colour="gray", size=0.25, aes(height = 1)) +
    scale_fill_gradientn(colours = c("#F5F5F5", "#E69F00","#0072B2"), values = c(0,lowest_W, scale_top*highest_W)) +
    #    scale_fill_gradient(low="white", high="turquoise3")+
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0))
}

# -------------------------------------------------------

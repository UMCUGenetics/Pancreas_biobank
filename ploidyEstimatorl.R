

readFile <- function(filename){
  #  filename <- list.files(pattern = filename, path = "/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/data", recursive = TRUE, full.names = T)
  first_line  <- readLines(filename, n=1)
  sample_data <- read.table(filename, header=FALSE, sep= "\t")
  names(sample_data) <- as.vector(unlist(strsplit(first_line[length(first_line)],"\t")))
  return(sample_data)
}

addMajorMinorAllele <- function(sample_data){
  sample_data$copyNumber[sample_data$copyNumber < 0] <- 0
  sample_data["unboundMinorAllele"] <- (1 - sample_data$actualBAF) * sample_data$copyNumber
  sample_data$unboundMinorAllele[sample_data$unboundMinorAllele < 0.25 ] <- 0
  sample_data["majorAllele"] <- sample_data$copyNumber - sample_data$unboundMinorAllele
  sample_data$majorAllele[sample_data$majorAllele < 0 ] <- 0
  return(sample_data)
}
chromArmPloidy <- function(sample_data, filename){
  majorPloidy_list <- c()
  minorPloidy_list <- c()
  lapply(c(1:22), function(c){
    if(c %!in% c(13,14,15,21,22)){
      subset_data <- subset(sample_data, sample_data$`#chromosome` == c)
      loc_centromeer <- which(subset_data$segmentEndSupport == "CENTROMERE")
      
      
      #Calculates everything for P arm
      subset_p  <- subset_data[1:loc_centromeer,]
      df_subset <- majorGroup(subset_p)
      
      ploidyFrameP <- createPloidyFrame(df_subset, c, "p")
      return_listP <- aneuPloidy(ploidyFrameP, c, "p", filename)
      
      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listP[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listP[2]
      
      
      #Calculates everything for Q arm
      subset_q   <- subset_data[(loc_centromeer+1):length(subset_data$start),]
      df_subset2 <- majorGroup(subset_q)
      
      ploidyFrameQ <- createPloidyFrame(df_subset2,c, "q")
      return_listQ <- aneuPloidy(ploidyFrameQ, c, "q", filename)
      
      majorPloidy_list[length(majorPloidy_list)+1] <<- return_listQ[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_listQ[2]
      
    }else{
      
      #Calculates everything for the autosomes that are treated as one armed
      subset_data <- subset(sample_data, sample_data$`#chromosome` == c)
      df_subset   <- majorGroup(subset_data)
      
      ploidyFrame <- createPloidyFrame(df_subset, c, "q")
      return_list <- aneuPloidy(ploidyFrame, c, "q", filename)
      
      majorPloidy_list[length(majorPloidy_list)+1] <<- return_list[1]
      minorPloidy_list[length(minorPloidy_list)+1] <<- return_list[2]
    }
  })
  return(list(majorPloidy_list, minorPloidy_list))
}


optimalization <- function(ploidy_list, filename, tp){
  ploidy_list  <- as.numeric(ploidy_list)
  common_value <- names(sort(summary(as.factor(ploidy_list)), decreasing=TRUE, na.last = TRUE))[1:2]
  genome_level <- common_value[1]
  ploidy_list  <- append(ploidy_list, values = genome_level)
  chromArms    <- getChromArms()
  tmp_df <- as.data.frame(cbind(chromArms, ploidy_list), stringsAsFactors = FALSE)
  
  doubt_df <- doubt_df[(doubt_df$file == filename),]
  doubt_df <- doubt_df[(doubt_df$tp == tp),]
  if(length(doubt_df$chromArm) == 0){
    return(ploidy_list)
  }else{
    for(i in 1:nrow(doubt_df)){
      row <- as.list(doubt_df[i,])
      if(row$highest == genome_level){
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- genome_level
      }else if(row$second == genome_level){
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- genome_level
      }else{
        tmp_df[which(tmp_df$chromArms == row$chromArm), "ploidy_list"] <- row$highest
      }
      tmp_df$ploidy_list <- as.numeric(tmp_df$ploidy_list)
    }
    ploidy_list <- tmp_df$ploidy_list
  }
  return(ploidy_list)
}


majorGroup <- function(data_frame){
  data_frame["majorGroup"] <- cut(data_frame$majorAllele,
                                  breaks= c(-1.5,-0.5,0.5,c(1*(1.5:20))), right = FALSE,
                                  include.lowest = TRUE, labels = c(-1:19))
  
  data_frame["group"] <- rleid(data_frame$majorGroup)
  data_frame$majorAllele[data_frame$majorAllele <= 0 ] <- NA
  data_frame["mapWeight"] <- rep(1, length(data_frame[,1]))
  data_frame$mapWeight[is.na(data_frame$mapWeight)] <- 1
  
  return(data_frame)
}

createPloidyFrame <- function(arm_subset, chromosome, arm){
  grouplist <- unique(arm_subset$group)
  
  tmp_list  <- lapply(1:length(unique(arm_subset$group)), function(i){
    group_set <- subset(arm_subset, arm_subset$group == grouplist[i])
    start <- group_set$start[1]
    end   <- group_set$end[length(group_set$end)]
    
    segment_size <- (group_set$end - group_set$start) + 1
    weight <-  segment_size * group_set$mapWeight
    
    averageMAP <- group_set$majorAllele * weight
    averageMAP <- sum(averageMAP) / sum(weight)
    
    averageMinor <- group_set$unboundMinorAllele * weight
    averageMinor <- sum(averageMinor) / sum(weight)
    
    averageCN <- group_set$copyNumber * weight
    averageCN <- sum(averageCN) / sum(weight)
    
    if(!is.na(averageMinor) && averageMinor <= 0.25){
      loh <- TRUE
    }else{
      loh <- FALSE
    }
    tmp_vector <- c(-1:19)
    return(c(paste0(chromosome,arm), as.integer(start), as.integer(end), as.numeric(averageMAP), as.numeric(averageMinor),
             as.numeric(averageCN), as.integer(tmp_vector[group_set$majorGroup[1]]), as.logical(loh)))
  })
  
  ploidyFrame <- as.data.frame(do.call(rbind, tmp_list), stringsAsFactors = FALSE)
  colnames(ploidyFrame) <- c("chromArm", "start", "end", "averageMajorAP", "averageMinorAP", "averageCN", "majorAPcategory", "lohStatus")
  ploidyFrame["minorAPcategory"] <- round(as.numeric(ploidyFrame$averageMinorAP))
  
  return(ploidyFrame)
}



aneuPloidy <- function(ploidyFrame, chromosome, arm, filename){
  ploidyFrame["segmentSize"] <- (as.numeric(ploidyFrame$end) - as.numeric(ploidyFrame$start))
  
  major_list <-  aggregate(segmentSize~majorAPcategory, data = ploidyFrame, FUN = sum)
  major_list["percentage"]   <- (major_list$segmentSize / sum(major_list$segmentSize)) * 100
  
  minor_list <-  aggregate(segmentSize~minorAPcategory, data = ploidyFrame, FUN = sum)
  minor_list["percentage"]   <- (minor_list$segmentSize / sum(minor_list$segmentSize)) * 100
  
  tmpList  <- list(major_list, minor_list)
  category <- c("major", "minor")
  
  counter <- 1
  ploidyList <- lapply(tmpList, function(subList){
    
    column <- paste0(category[counter], "APcategory")
    if(length(order(subList$percentage)) > 1){
      if(max(subList$percentage) >= 50  && (subList[which.max(subList$percentage), "percentage"] - sort(subList$percentage, decreasing = TRUE)[2]) > 10){
        ploidyLevel <- subList[which.max(subList$percentage), column]
      }else{
        ploidyLevel <- NA
        doubt_df[nrow(doubt_df)+1,] <<- c(filename, paste0(chromosome, arm),  subList[which.max(subList$percentage), column],
                                          subList[which(subList$percentage == (sort(subList$percentage,decreasing=TRUE)[2])), column], category[counter])
      }
    }else{
      ploidyLevel <- subList[which.max(subList$percentage), column]
    }
    counter <<- counter + 1
    return(ploidyLevel)
  })
  
  return(list(as.numeric(unlist(ploidyList[1])), as.numeric(as.character(unlist(ploidyList[2])))))
}

getChromArms <- function(){
  chromArms = c()
  for(i in c(1:22)){
    if(i %!in% c(13,14,15,21,22)){
      chromArms <- c(chromArms, c(paste0(i,"p"), paste0(i,"q")))
    }else{
      chromArms <- c(chromArms, c(paste0(i,"q")))
    }
  }
  chromArms <- append(chromArms, "genome_level")
  return(chromArms)
}









#Imports
library("data.table")

#Global settings, functions and variables
options(stringsAsFactors = FALSE)
options(scipen = 999)


'%!in%' <- function(x,y)!('%in%'(x,y))
doubt_df <- data.frame(file = character(),chromArm = character(), highest= integer(), second= integer(), tp = character(),
                       stringsAsFactors = FALSE)


ploidyEstimator <- function(path){
  files <- list.files(path = path,  pattern = "\\.purple.cnv$", recursive = TRUE, full.names = TRUE)
  
  chromArms <- getChromArms()
  df_majorPloidy <- data.frame(chromArms)
  df_minorPloidy <- data.frame(chromArms)
  df_copyNumber  <- data.frame(chromArms)
  
  return_list <- lapply(files, function(file){
    #file="/Users/avanhoeck/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data//FR11124015.purple.cnv"
    file_name_CN <- unlist(strsplit(file, "CN_data//"))[2]
    sample_data <- readFile(file)
    sample_data <- addMajorMinorAllele(sample_data)
    return_list <- chromArmPloidy(sample_data, file)
    
    majorPloidyList <- optimalization(unlist(return_list[1]), file, "major")
    minorPloidyList <- optimalization(unlist(return_list[2]), file, "minor")
    
    df_majorPloidy[file_name_CN] <<- as.numeric(majorPloidyList)
    df_minorPloidy[file_name_CN] <<- as.numeric(minorPloidyList)
    df_copyNumber[file_name_CN]  <<- as.numeric(as.character(majorPloidyList)) + as.numeric(as.character(minorPloidyList))
    
  })
  
  write.table(df_majorPloidy, file = "Results_MajorAP.tsv"    , sep="\t", col.names = TRUE, row.names = FALSE)
  write.table(df_minorPloidy, file = "Results_MinorAP.tsv"    , sep="\t", col.names = TRUE, row.names = FALSE)
  write.table(df_copyNumber , file = "Results_CopyNumber.tsv" , sep="\t", col.names = TRUE, row.names = FALSE)
}

setwd("~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/")
ploidyEstimator(path = "~/surfdrive/projects/SU2C/Colorectal-Pancreas/Else_analyses/CN_data/")











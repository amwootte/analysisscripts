#!/bin/Rscript
#Calc.Generic.Skill
#Calculates a skill score for an input dataset that is ds-truth
#We might be adding more skill scores as time goes on

CalcGenericSkill <- function(input_pred, var_pred, input_ref, var_ref, 
                             input_perf, var_perf='na', err_stat='mae', 
                             output='not needed yet', args=list()){
  start.time <- proc.time()
  message(paste("Script started at:", Sys.time()))
  
  pred.data <- build.fudge.object(input_pred,paste0(var_pred,"_",err_stat), 
                                  dim=c('spatial', 'temporal', 'ensemble'))
  ref.data  <- build.fudge.object(input_ref, paste0(var_ref, "_",err_stat))
  #Similarity checks on pred.data and ref.data should go in here
  
  if (is.na(suppressWarnings(as.numeric(input_perf)))){
    #This is a character string, not a number
    perf.data  <- build.fudge.object(input_perf, paste0(var_perf, "_",err_stat))
    perf.string <- paste("located in", normalizePath(input_perf))
  }else{
    #fill array of same dimensions as input with input_perf number
    perf.data <- list(data=list(as.numeric(input_perf)))
    perf.string <- paste("of",input_perf)
  }
  
  #Determine output metadata, based in part off of the err_stat
  err_prop <- switch(err_stat,
                     'mae' = list('Mean Absolute Error'), 
                     'rmse' = list('Root Mean Square Error'), 
                     'mse' = list('Mean Square Error'),
                     stop(paste("Error in CalcGeneralSkill: error statistic", 
                                err_stat, 'was not one of supported choices!')))
  
  out.data <- list()
  #Calculate skill score
  out.data$skill_score <- calc.wilks.skill(pred.data$data[[1]], ref.data$data[[1]], perf.data$data[[1]])
  #It's not technically neccessary to make out.data a list of arrays in this case - 
  #there is only a single variable. If you write more than one variable's worth of
  #output to a single file, however, you need to write each variable as a separate
  #component of the list.
  
  #Write to file
  elapsed.time <- proc.time()-start.time
  message(paste("Calculation took", elapsed.time[3], "seconds to complete"))
  
  #Setting parameters for writing to file
  units.all <- "1"
  precs <- "float"
  var.longname <- ref.data$metadata[[1]]$long_name
  stat.longname <- err_prop[[1]]
  var.longnames <- paste("Skill score based on", var.longname)  
  
  #Metadata
  output.comment <- sub("  ", " ", paste("Prediction skill for the file", normalizePath(input_pred),
                                         "using reference", normalizePath(input_ref), 
                                         "and perfect prediction", perf.string,
                                         "with CalcGenericSkill 1.0", sep=" ")
  )
  append.globals <- list(history=paste(Sys.time(), get.git.branch(base.path),
                                       "Rscript CalcGenericSkill.R", 
                                       get.input.args(args$in.args)), 
                         comment=output.comment, 
                         institution='NOAA-GFDL', 
                         title=paste("Prediction skill measurment for", basename(input_pred)),
                         git_branch = paste("Run with", get.git.branch(base.path), 
                                            "at", Sys.time()
                         ), 
                         skill_type = "generic skill score")
  print(paste("writing to file", output))
  new.out.file <- WriteNC(output, out.data,
                          names(out.data), 
                          dim.list=pred.data$metadata[[1]]$dim,
                          var.data=pred.data$metadata[[1]]$vars,
                          units="skill units",
                          longname=var.longnames,
                          prec="float", 
                          clone_globals=F, global_clone_file=input_pred,
                          global_append = append.globals,
                          verbose=args$verbose)
  
  end.time <- proc.time()-start.time
  message(paste("Script took", end.time[3], 'seconds to run.'))
  message(paste("File", paste(output, collapse=""), "written.", sep=" "))  
}

############## HELPER FUNCTIONS ################################################
calc.wilks.skill <- function(pred, ref_pred, perf=0){
  #Perf should be of same dimensions as data and ref
  skill <- ( (pred-ref_pred)/(perf-ref_pred) ) * 100
  #replace non-numeric and Inf values
  return(skill)
}

read.timediffstats.file <- function(file, err_stat){
  test.nc <- nc_open(file)
  var.names <- names(test.nc$var)
  err.var.names <- var.names[grepl(paste0("_",err_stat), var.names)]
  out.data <- list(data=list(),metadata=list())
  for (v in err.var.names){
    #replace this with a call to the fudgeio primitives?
    out.data$data[[v]] <- ncvar_get(test.nc, v)
  }
  nc_close(test.nc)
  return(out.data)
}
################################################################################

base.path <- Sys.getenv("S8_PATH")
if(base.path==""){
  #If S8_tools not set, assume that you are working from the correct directory
  base.path <- getwd()
}
out <- source(paste0(base.path,'/s8_utils/LoadRPackages.R'))
out <- sapply(list.files(pattern="[.]R$", path=paste(base.path,'/FudgeIO_reboot/',sep=''), full.names=TRUE), source);
out <- sapply(list.files(pattern="[.]R$", path=paste0(base.path,'/s8_utils/'), full.names=TRUE), source, .GlobalEnv);

#Parsing function
ParseSkillScoreArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c("--input_pred"), action="store", default='na',
                help=paste("Input file produced by CalcTimeDiffStats1 containing",
                           "error statistics for the prediction being evaluated.")),
    make_option(c("--var_pred"), action="store", default='na',
                help=paste("Original variable of the prediction; if", 
                           "${var_pred}_${err_stat} is not found in input_pred",
                           "the script will throw an error.")),
    make_option(c("--input_ref"), action="store", default='na',
                help=paste("Input file produced by CalcTimeDiffStats1 containing",
                           "error statistics for the reference dataset used", 
                           "in the evaluation")),
    make_option(c("--var_ref"), action="store", default='na',
                help=paste("Original variable of the reference dataset; if", 
                           "${var_ref}_${err_stat} is not found in input_ref",
                           "the script will throw an error.")),
    make_option(c("--input_perf"), action="store", default='0',
                help=paste("Either an input file produced by CalcTimeDiffStats1",
                           "containing error statistics for a perfect prediction", 
                           "or a single numeric value representing the error", 
                           "of the prediction (i.e. 0). Defaults to 0.")),
    make_option(c("--var_perf"), action="store", default='na',
                help=paste("If --input_pref is a perfect prediction file,", 
                           "the original variable of the prediction file.", 
                           "${var_pred}_${err_stat} is not found in input_pred",
                           "the script will throw an error; if $input_pred", 
                           "is a single numeric value, this variable is never", 
                           "used.")),    
    make_option(c("-o", "--output"), action="store", default="",
                dest='output',
                help=paste("Output from the script; will be in the form of filename.nc")),    
    make_option(c("--err_stat"), action='store', default='mae', 
                help=paste("Which error statistic in the diffstats files to use", 
                           "for calculation of the skill score. Can be one of",
                           "the following:",
                           "\n\t1) mae, Mean Absolute Error",
                           "\n\t2) rmse, Root Mean Square Error",
                           "\n\t3) mse, Mean Square Error (MSE not present in earlier versions of CalcTimeDiffStats1)",
                           "\n\tIf not specified, defaults to mae.")),
    make_option(c("--prof"), action='store', default=0, 
                help=paste("Profiling information for the function as a whole.", 
                           "Prints no information for 0, total memory consumption",
                           "for 1, and total memory plus subfunction timing for 2.",
                           "Defaults to 0.")),
    make_option(c("--verbose"), action="store_true", default=FALSE,
                help=paste("Whether to print extra status messages to facilitate debugging.", 
                           "Defaults to FALSE (no verbose messages added). "))
  )
  description = paste('Given two files produced by CalcTimeDiffStats1,', 
                      'a variable of interest in each file, an error statistic,', 
                      'and either a single value or another set of files and vars',
                      'computes the skill score for the predictions',
                      'in the input_pred file. Output will have the same',
                      'i,j,ens,t dimensions as its input files (all inputs',
                      'require the same i,j,ens,t) and a single output',
                      'var, skill_score.')
  epilouge = "Please note: flags may be specified in any order, and '=' not required to specify strings."
  usage = paste("usage: %prog --input_pred input_pred.nc --var_pred var_pred --input_ref input_ref.nc",
                "--var_rev var_ref --input_perf 0 [--var_perf var_perf]",
                "--err_stat mae|rmse|mse|se --output output.nc",
                "[--verbose] [-h --help]")
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge)))
}
#Main method of script and argument parsing
args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)
if (arg.len > 0){
  parsed.args <- ParseSkillScoreArgs(args)
  parsed.args <- CheckStatArgs(parsed.args)
}else{
  #Assume that you are running from R, and want to test results
  parsed.args <- list()
  parsed.args$input_pred <- gcp.tmpdir("/net3/kd/PROJECTS/DOWNSCALING/DATA/CExplr/diagnostics/tasmax/CalcTimeDiffStats1/tasmax_day_CalcTimeDiffStats1_EDQM_cm3_8p5-minus-c360hires_cm3_8p5_CONUS_12mask.nc")
  parsed.args$var_pred <- "tasmax"
  parsed.args$input_ref <- gcp.tmpdir(file.path("/net3/kd/PROJECTS/DOWNSCALING/DATA/CExplr/diagnostics/tasmax/CalcTimeDiffStats1/",
                                                "tasmax_day_CalcTimeDiffStats1_c360coars_cm3_8p5-minus-c360hires_cm3_8p5_CONUS_12mask.nc"))
  parsed.args$var_ref <- "tasmax"
  parsed.args$input_perf <- 0
  parsed.args$var_perf <- 'NA'
  parsed.args$err_stat <- "mae"
  parsed.args$output = "/home/cew/Code/testing/skill_analysis.nc"
  parsed.args$verbose <- TRUE
  parsed.args$prof <- 3
}
#Include args for metadata
parsed.args$args <- parsed.args

#Any var checks needed?

profile(parsed.args$prof, CalcGenericSkill, 
        parsed.args$input_pred, parsed.args$var_pred, 
        parsed.args$input_ref, parsed.args$var_ref, 
        parsed.args$input_perf, parsed.args$var_perf, 
        parsed.args$err_stat, 
        parsed.args$output, args=parsed.args$args)
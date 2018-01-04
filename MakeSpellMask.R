#! /usr/local/R-2.15.2/bin/Rscript --vanilla --default-packages=utils
#' MakeSpellMask.R
#' given a file, threshold, and comparison, counts the number of hits along the 
#' time and ensemble axes at each i,j point that match the threshold 
#' and comparison (i.e ==, > <=) for at least the given # of days.

MakeSpellMask <- function(input, var, cond, thold, out.file, 
                          units='na', spell_len=7, verbose=FALSE, args=list()){
  start.time <- proc.time()
  message(paste("Script started at:", Sys.time()))
  
  #Read in data and convert units
  message(paste("Reading in the input file", input))
  data <- build.fudge.object(input, var, dim=c('spatial', 'temporal', 'ensemble'))
  data.dim<- dim(data$data[[var]])
  data.units <- data$metadata[[var]]$units
  if (units != "na" && units != data.units){
    message("Converting units")
    orig.thold <- thold
    thold <- convert.units(thold, unit1=units, unit2=data.units)
    if (args$eps != "na"){
      args$eps <- convert.units(args$eps, unit1=units, unit2=data.units)
    }
  }else{
    orig.thold <- thold
    units <- data.units
  }
  gc()
  data.prec <-attr(data$clim.in, "prec")  
  
  var.names <- paste0(var, '_spell')
  cap.names <- paste("Spells of points", cond, orig.thold, units, 
                     "at least", spell_len, 'tlevels long')
  var.longnames <- paste(cap.names) 
  
  ###Okay, control of the looping structure needs to take place down here
  message("Calculating spell masking")
  thold.time <- proc.time()
  out.data <- array(NA, dim=data.dim)
  #The advantage to forcing all inputs into 6 dimensions is that you don't have
  #to worry about ensemble dimensions
  #dim order: x, y, z, ens, n, t
  for (i in 1:data.dim[1]){
    for(j in 1:data.dim[2]){
      for(n in 1:data.dim[4]){
        tmp.data <- as.numeric(count_rle_thold(data$data[[var]][i,j,,n,,], 
                                    thold, cond, spell_len, 
                                    args$cond.args))
        tmp.data[tmp.data==0] <- NA
        out.data[i,j,,n,,] <-  tmp.data
      }
    }
  }
  print(paste("Calculation took", proc.time()[3]-thold.time[3], "seconds."))
all.dimlist <- data$metadata[[1]]$dim
all.varlist <- data$metadata[[1]]$vars
dim.calc.string <- "Mask of t points at each ens,x,y location"

#Metadata
output.comment <- sub("  ", " ", 
                      paste(dim.calc.string, "in", 
                            basename(input), 
                            "where", var, 'was', cond, thold, 
                            "for at least", spell_len, "timelevels",
                            'at', Sys.time(), sep=" "))
dataset.name <- get.dataset.name(input)
append.globals <- list(history=paste(Sys.time(), get.git.branch(base.path),
                                     "Rscript MakeSpellMask.R", 
                                     get.input.args(args$in.args)), 
                       comment=output.comment, 
                       institution='NOAA-GFDL', 
                       title=sub("  ", " ", 
                                 paste("Count of", var, cond, orig.thold, units, 
                                       "for at least", spell_len, "days",
                                       dataset.name)),
                       git_branch = paste("Run with", get.git.branch(base.path), 
                                          "at", Sys.time()
                       ), 
                       skill_type = "generic skill score")

    if(verbose){message(paste("Starting write to to file", out.file))}

print(paste("writing to file", out.file))
new.out.file <- WriteNC(out.file, out.data,
                        var.names, 
                        dim.list=all.dimlist,
                        var.data=all.varlist,
                        units="count",
                        longname=var.longnames,
                        prec="float", 
                        clone_globals=F, global_clone_file=input_pred,
                        global_append = append.globals,
                        verbose=args$verbose)
    message(paste('Variables', paste(var.names, collapse=", "), "written to file"))
  end.time <- proc.time() - start.time
  message(paste("Script took", end.time[3], 'seconds to run.'))
  message(paste("File(s)", paste(out.file, collapse=""), "written.", sep=" "))
}
######## INTERNAL HELPER FUNCTIONS ########
count_thold <- function(data, thold, cond, cond.args){
  #'Helper function for tallying number of counts
  if(any(!is.na(data))){
    if(!is.null(cond.args)){
      return(sum(do.call(cond, list(data, thold, cond.args)), na.rm=TRUE))
    }else{ return(sum(do.call(cond, list(data, thold)), na.rm=TRUE))}
  }else{return(NA)}
}

count_rle_thold <- function(data_vec, thold, cond, spell_len = 1, cond.args=NULL){
  #'Helper function for tallying the spells from the 
  #'threshold passed in
  if(any(!is.na(data_vec))){
    if(!is.null(cond.args)){
      data.cond <- do.call(cond, list(data_vec, thold, cond.args))
    }else{ 
      data.cond <- do.call(cond, list(data_vec, thold))
    }
    #Calculate any runs that pass greater than the spell length
    cond.rle <- rle(data.cond)
    inv.spell.cond <- (cond.rle$values & cond.rle$lengths < spell_len)
    cond.rle$values[inv.spell.cond] <- F
    out.data <- inverse.rle(cond.rle)
    return(out.data)
  }else{return(NA)}  
}
#######INPUT PARSING#######
ParseSpellMaskArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c('-i', "--input"), action="store", default='na',
                dest='input',
                help=paste("The input netCDF file to observe.")),
    make_option(c('-v', "--var"), action="store", default='na',
                dest='var',
                help=paste("The variable of interest in the input file.", 
                           "If not present, an error will be thrown.")),
    make_option(c("-o", "--output"), action="store", default="",
                dest='output',
                help=paste("Output from the script;", 
                           "will throw an error if .nc is not present in the filename.")),    
    make_option(c("--threshold"), action="store", default='na',
                help=paste("The threshold that will be applied to the input", 
                           "in a manner determined by $condition.",
                           "Can be any number greater than zero.",
                           "If not present, an error will be thrown.")),
    make_option(c("--condition"), action="store", default='na',
                help=paste("The condition that determines how the threshold will", 
                           "be applied to the input data. May be one of LT, LE,", 
                           "GT, GE, EQ (for less than, less than or equal, etc.", 
                           "or the corresponding symbols (< <= > >= =).", 
                           "For the sixth option, AQ or ~= (approximately", 
                           "equal to), the parameter eps also must be set")),
    make_option(c("--eps"), action="store", default='na',
                help=paste("If condition is AQ (or ~=), the +/- epsilon to which", 
                           "the data will be matched. Has no effect for the other conditions,", 
                           "and is in the same units as threshold.")),
    make_option(c("--units"), action='store', default='na',
                help=paste("The units in which the threshold and epsilon are be calculated,", 
                           "presented in a udunits-parseable form.", 
                           "(For more information on accepted units, please consult:", 
                           "cat /home/esd/Code/section_8_tools/udunits_list.txt)",
                           "If not specified, defaults to the units of the input file.")),
    make_option(c("--spell-length"), action='store', dest='spell_length', default=7, 
                help=paste("The minimum length of a spell that you want to search", 
                           "for in the input dataset - a minimum number of time levels", 
                           "for which the threshold and condiation are met.", 
                           "Defaults to 7.")),
    make_option(c("--prof"), action='store', default=0, 
                help=paste("Profiling information for the function as a whole.", 
                           "Prints no information for 0, total memory consumption",
                           "for 1, and total memory plus subfunction timing for 2.",
                           "Defaults to 0.")),
    make_option(c("--verbose"), action="store_true", default=FALSE,
                help=paste("Whether to print extra status messages to facilitate debugging.", 
                           "Defaults to FALSE."))
  )
  description = paste('Given an input file, a condition, and a threshold, applies', 
                      "the threshold to the input file in the manner determined",
                      "by condition, and counts the number of t points at",
                      "each i,j,ens location (or the number of t,ens at each", 
                      "i,j if pool_ens is set) that meet the condition. Output", 
                      "will have the i,j,ens or i,j dimensions of the input file.")
  epilouge = paste("Please note: flags may be specified in any order, and '='", 
                    "not required to specify strings.")
  usage = paste("usage: %prog --input input.nc --var var", 
                "--condition GT | GE | LT | LE | EQ | AQ", 
                "--threshold threshold -O out_base --eps eps [--units units] --spell-length", 
                "[--verbose] [-h --help]")
  arg.list <- allow.neg.commandline.opts(arg.list, '--threshold')
  return(parse_args(OptionParser(option_list=option_list, usage=usage, 
                                 description = description, epilogue=epilouge), 
                    args=arg.list))
}

###### MAIN METHOD ###### 
base.path <- Sys.getenv("S8_PATH")
if(base.path==""){
  #If S8_tools not set, assume that you are working from the correct directory
  base.path <- getwd()
}
#Note: if you use this strucutre, it is VITAL that you do not include ANY scripts in the directories that
#you are sourcing. Those scripts will run, and you will be frustrated.
out <- source(paste0(base.path,'/s8_utils/LoadRPackages.R'))
out <- sapply(list.files(pattern="[.]R$", path=paste(base.path,'/FudgeIO_reboot/',sep=''), full.names=TRUE), source);
out <- sapply(list.files(pattern="[.]R$", path=paste0(base.path,'/s8_utils/'), full.names=TRUE), source, .GlobalEnv);

args <- commandArgs(trailingOnly = TRUE)
arg.len <- length(args)
if(arg.len != 0){
  parsed.args <- ParseSpellMaskArgs(args)
  parsed.args <- CheckStatArgs(parsed.args)
}else{
  args <- list()
  parsed.args <- list()
  parsed.args$input <- "/home/cew/Code/testing/test_minfiles/9_pt_tasmax.nc"
  parsed.args$var <- 'tasmax'
  parsed.args$threshold <- 280
  parsed.args$units <- 'K'
  parsed.args$condition <- 'LT' #AQ
  parsed.args$eps <- 0.5
  parsed.args$spell_length <- 7
  parsed.args$output <- "/home/cew/Code/testing/spell_test_i.nc"
  parsed.args$verbose <- TRUE
  
  parsed.args <- CheckStatArgs(parsed.args)
  parsed.args$in.args <- args
  parsed.args$prof <- 0

}
###Add method-specific checks here

#Check on presence of var
if(parsed.args$var=='na'){
  stop(paste("Error in MakeSpellMask: --var not specified"))
}

#parse on condition 
cond <- switch(parsed.args$condition, 
               "GT" = '>', 
               "LT" = "<", 
               "GE" = ">=", 
               "LE" = '<=', 
               "EQ" = "==", "=" = "==", 
               "AQ" = 'approx.eq', '~=' = 'approx.eq')
if(!any(cond==c('==', "<", ">", "<=", ">=", 'approx.eq'))){
  stop(paste("Error in MakeSpellMask: condition", 
             parsed.args$condition, "was not of an accepted type!"))
}else{
  parsed.args$condition <- cond
  if(parsed.args$condition=='approx.eq'){
    if(parsed.args$eps=='na'){
      stop(paste("Error in MakeSpellMask: condition AQ or ~=", 
                 "requires eps to be set!"))
    }else{
      parsed.args$cond.args <- as.numeric(parsed.args$eps)
    }
  }else{
    parsed.args$cond.args=NULL
  }
}

#If passing in command-line args, this needs to be checked
thold.string <- gsub("[()]", "", parsed.args$threshold)
parsed.args$threshold <- suppressWarnings(as.numeric(thold.string))
if(is.na(parsed.args$threshold)){
  stop("Error in MakeSpellMask! Threshold could not be interpreted as numeric!")
}

#And make sure spell_length is a positive number
if (parsed.args$spell_length < 1){
  stop(paste("Error in MakeSpellMask: spell_length should be >= 0!"))
}


#Function call
profile(parsed.args$prof, MakeSpellMask, parsed.args$input, parsed.args$var,
              parsed.args$condition, parsed.args$threshold,
              parsed.args$output, parsed.args$units,
        parsed.args$spell_length,
              parsed.args$verbose, parsed.args)

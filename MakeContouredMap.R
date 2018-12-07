MakeContouredMap <- function(inputfiles,stat_used,var,colorchoice,spatmask.file='na',units='na',bin_limit,outfilepath,outfile_base="testgraphic",outfile_by_filename=FALSE,use_fixed_scale=FALSE,fixed_scale='na',type="raw",useseason=FALSE,season=NA){
  
  ##################
  # Input check
  
  if(outfile_by_filename==FALSE & outfile_base=='na'){
    stop("Error in MakeContouredMap: the outfile_base argument must be set if outfile_by_filename is FALSE!")
  }

  if(useseason==TRUE){
     if(season=="DJF") seasonidx=1
     if(season=="MAM") seasonidx=2
     if(season=="JJA") seasonidx=3
     if(season=="SON") seasonidx=4
  }
  
  ###################
  # Gathering inputs
  
  message(paste("Reading in the input files"))
  if(length(inputfiles)==1){
    ##############
    # For only reading 1 file
    
    data.1 <- build.fudge.object(inputfiles, var, dim=c('spatial'),verbose=TRUE)
    #Check for identical dimensions
    check.result <- CheckIdenticalDims(data.1$dim, data.1$dim)
    if(check.result!=""){
      stop(paste("Error in MakeContouredMap: File",attr(data.1, 'filename'), "and file", attr(data.1, 'filename'),
                 "have", paste(check.result), "dimensions that do not have the same values.")) }  
    lat = data.1$metadata[[var]]$vars$lat_bnds[[16]][[2]]$vals
    lon = data.1$metadata[[var]]$vars$lon_bnds[[16]][[2]]$vals
    data.units <- data.1$metadata[[var]]$units
    
    if (units != "na" && units != data.units){
      message("Converting units")
      if(spatmask.file!="na"){
        message("applying spatial mask")
        tmp1 <- apply.any.mask(data = eval(parse(text=paste("data.",i,"$data[[var]]",sep=""))), 
                               mask = spatmask$masks[[1]], 
                               dim.apply = 'spatial')
        data.1$data[[var]] <- convert.units(tmp1, unit1=data.units, unit2=units)
      } else {
        data.1$data[[var]] <- convert.units(eval(parse(text="data.1$data[[var]]")), unit1=data.units, unit2=units)
      }
      
    }else{
      if(spatmask.file!="na"){
        message("applying spatial mask")
        data.1$data[[var]] <- apply.any.mask(data = eval(parse(text=paste("data.",i,"$data[[var]]",sep=""))), 
                       mask = spatmask$masks[[1]], 
                       dim.apply = 'spatial')
      }
      units <- data.units
    }

    if(use_fixed_scale==TRUE){
    if(length(which(data.1$data[[var]]>fixed_scale[2]))>0){data.1$data[[var]][which(data.1$data[[var]]>fixed_scale[2])]=fixed_scale[2]}
    if(length(which(data.1$data[[var]]<fixed_scale[1]))>0){data.1$data[[var]][which(data.1$data[[var]]<fixed_scale[1])]=fixed_scale[1]} 
    }

      if(useseason==FALSE){
        vals = data.1$data[[var]]
      } else {
        vals=  data.1$data[[var]][,,1,1,1,seasonidx]
      }

    } else {
      
      
      ############
      # for reading multiple files
      
      message("Reading multiple input files")
    message("Working with ",length(inputfiles)," total")
    for(i in 1:length(inputfiles)){
      message("working on file ",inputfiles[i])
      assign(paste("data.",i,sep=""),build.fudge.object(inputfiles[i], var, dim=c('spatial'))) 
    }
    #Check for identical dimensions
    check.result <- CheckIdenticalDims(data.1$dim, data.2$dim)
    if(check.result!=""){
      stop(paste("Error in MakeContouredMap: File",attr(data.1, 'filename'), "and file", attr(data.2, 'filename'),
                 "have", paste(check.result), "dimensions that do not have the same values.")) 
    }
    
    
    lat = data.1$metadata[[var]]$vars$lat_bnds[[16]][[2]]$vals
    lon = data.1$metadata[[var]]$vars$lon_bnds[[16]][[2]]$vals
    vals = c()
    
    for(i in 1:length(inputfiles)){
      data.units <- data.1$metadata[[var]]$units
      #message("Value of spatmask.file= ",spatmask.file)
      if (units != "na" && units != data.units){
        message("Converting units")
        
        if(spatmask.file!='na'){# spatial mask and unit conversion
          message("applying spatial mask")
          if(i==1) spatmask <- ReadMaskNC(nc_open(spatmask.file))
          tmp1 = apply.any.mask(data = eval(parse(text=paste("data.",i,"$data[[var]]",sep=""))), 
                                mask = spatmask$masks[[1]], 
                                dim.apply = 'spatial')
          tmp = convert.units(tmp1, unit1=data.units, unit2=units)
        } else { # no spatial mask with unit conversion
          tmp = convert.units(eval(parse(text=paste("data.",i,"$data[[var]]",sep=""))), unit1=data.units, unit2=units)
        } 
 	if(use_fixed_scale==TRUE){
    	if(length(which(tmp>fixed_scale[2]))>0){tmp[which(tmp>fixed_scale[2])]=fixed_scale[2]}
    	if(length(which(tmp<fixed_scale[1]))>0){tmp[which(tmp<fixed_scale[1])]=fixed_scale[1]} 
    	}

        assign(paste("dvals",i,sep=""),tmp)  
      
        }else{
        
        if(spatmask.file!='na'){ # spatial mask and no unit conversion
          if(i==1) spatmask <- ReadMaskNC(nc_open(spatmask.file))
          message("applying spatial mask")
          tmp = apply.any.mask(data = eval(parse(text=paste("data.",i,"$data[[var]]",sep=""))), 
                                mask = spatmask$masks[[1]], 
                                dim.apply = 'spatial')  
        } else {  # no mask and no unit conversion
	  tmp = eval(parse(text=paste("data.",i,"$data[[var]]",sep="")))          
        } 
        
        if(use_fixed_scale==TRUE){
    	if(length(which(tmp>fixed_scale[2]))>0){tmp[which(tmp>fixed_scale[2])]=fixed_scale[2]}
    	if(length(which(tmp<fixed_scale[1]))>0){tmp[which(tmp<fixed_scale[1])]=fixed_scale[1]} 
    	}
        assign(paste("dvals",i,sep=""),tmp)

        units <- data.units
        
      }
      if(useseason==FALSE){
        vals = c(vals,eval(parse(text=paste("dvals",i,sep=""))))
      } else {
        vals = c(vals,eval(parse(text=paste("dvals",i,"[,,1,1,1,",seasonidx,"]",sep=""))))
      }

    }
     
  }
  
  ##########################
  # Determine Color bars
  message("Creating color ramp")
  if(spatmask.file!="na"){
    vals = vals[-which(is.na(vals==TRUE))]
  }
  
  colorbar_info=colorramp(inputdata=vals,colorchoice=colorchoice,Blimit=bin_limit,use_fixed_scale=use_fixed_scale,fixed_scale=fixed_scale,type=type)

  ##########################
  # set up plotting margins
  
 #set the dimensions of the output image (plot + mar), if needed
  mar.default <- c(1, 5, 4, 4)  + 0.1
  mar.set <- mar.default
  message("writing plots")
  if(length(inputfiles)==1){
    split1 = do.call("c",strsplit(eval(parse(text=paste("data.",1,"$metadata[[var]]$filename",sep=""))),"/",fixed=TRUE))
    split_filename = do.call("c",strsplit(split1[length(split1)],"_",fixed=TRUE))

    if(split_filename[3]=="GFDL-HIRAM-C360-COARSENED"){
     DStech = "c360coars"
    } else {
     split2 = do.call("c",strsplit(split_filename[3],"-",fixed=TRUE))    
     DStech = split2[2]
    }
    if(split_filename[3]=="GFDL-HIRAM-C360") DStech="c360hires"
    varname = split_filename[1]
    
    if(split_filename[4]=="amip") epoch="hist"
    if(split_filename[4]!="amip"){
      if(split_filename[4]=="sst2030") scen = "4p5"
      if(split_filename[4]=="sst2090") scen = "8p5"
      
      if(split_filename[5]=="r1to2i1p1" | split_filename[5]=="r1to3i1p1") GCMv = "esm"
      if(split_filename[5]=="r1to2i1p2" | split_filename[5]=="r1to3i1p2") GCMv = "cm3"
	if(split_filename[5]=="DSCVh") GCMv="DSCVh"
      epoch = paste(GCMv,scen,sep="_")
    } 
    message("Creating plot titles")
    
    ptitle = paste("Difference Stats for ",varname,sep="")
    xlab="LONGITUDE"
    ylab="LATITUDE"
    if(stat_used!="raw") extratext = paste(varname,"_",stat_used," ",DStech,"_-minus-_c360hires ",epoch,sep="")
	if(stat_used=="raw") extratext = paste(varname," ",DStech,"_",epoch,sep="")
      
    if(outfile_by_filename==FALSE){
      outfilename = paste(outfile_base,".png",sep="")
    } else {
      if(stat_used!="raw") outfilename = paste(varname,"_",DStech,"_",epoch,"_-minus-c360hires-_",epoch,"_",stat_used,".png",sep="")
	if(stat_used=="raw") outfilename = paste(varname,"_",DStech,"_",epoch,"_-minus-c360hires-_",epoch,"_",stat_used,".png",sep="")
    }
    
if(useseason==FALSE){
sfc = list(x=lon-360,y=lat,z=eval(parse(text=paste("dvals",1,"[,,1,1,1,1]",sep=""))))
} else {
sfc = list(x=lon-360,y=lat,z=eval(parse(text=paste("dvals",1,"[,,1,1,1,",seasonidx,"]",sep=""))))
}    

    png(paste(outfilepath,outfilename,sep=""),width=600,height=480,units="px",pointsize=12)
    surface(sfc,type="I",main=ptitle,xlab=xlab,ylab=ylab,zlim=colorbar_info[[1]],breaks=colorbar_info[[2]],col=colorbar_info[[3]])#
    map("world",add=TRUE)
    mtext(extratext,side=1,line=-1.5)
    dev.off()
    
    
  } else {
    
    message("Printing multiple plots")
    for(i in 1:length(inputfiles)){
      
      message("Working on plot ",i)
      split1 = do.call("c",strsplit(eval(parse(text=paste("data.",i,"$metadata[[var]]$filename",sep=""))),"/",fixed=TRUE))
      split_filename = do.call("c",strsplit(split1[length(split1)],"_",fixed=TRUE))
    if(split_filename[3]=="GFDL-HIRAM-C360-COARSENED"){
     DStech = "c360coars"
    } else {
     split2 = do.call("c",strsplit(split_filename[3],"-",fixed=TRUE))    
     DStech = split2[2]
    }
    if(split_filename[3]=="GFDL-HIRAM-C360") DStech="c360hires"
    varname = split_filename[1] # for the filenames bits, perfect model regular, 3,4,5 for split_filename, sometimes 4,5,6, 
				# for MAE ratio of the CLIMDEX variables varname:idx=1, DStech:idx=3,epoch:idx=4,GCMv:idx=5
      
      if(split_filename[4]=="amip") epoch="hist"
      if(split_filename[4]!="amip"){
        if(split_filename[4]=="sst2030") scen = "4p5"
        if(split_filename[4]=="sst2090") scen = "8p5"
        
        if(split_filename[5]=="r1to2i1p1" | split_filename[5]=="r1to3i1p1") GCMv = "esm"
        if(split_filename[5]=="r1to2i1p2" | split_filename[5]=="r1to3i1p2") GCMv = "cm3"
	if(split_filename[5]=="DSCVh") GCMv="DSCVh"
        epoch = paste(GCMv,scen,sep="_")
      } 
      
      ptitle = paste("Difference Stats for ",varname,sep="")
      xlab="LONGITUDE"
      ylab="LATITUDE"
      if(stat_used!="raw") extratext = paste(varname,"_",stat_used," ",DStech,"_-minus-_c360hires ",epoch,sep="")
	if(stat_used=="raw") extratext = paste(varname," ",DStech,"_",epoch,sep="")
      
      if(outfile_by_filename==FALSE){
        outfilename = paste(outfile_base,"_",i,".png",sep="")
      } else {
        if(stat_used!="raw") outfilename = paste(varname,"_",DStech,"_",epoch,"_-minus-c360hires-_",epoch,"_",stat_used,".png",sep="")
	if(stat_used=="raw") outfilename = paste(varname,"_",DStech,"_",epoch,"_-minus-c360hires-_",epoch,"_",stat_used,".png",sep="")
      }
      if(useseason==FALSE){
	sfc = list(x=lon-360,y=lat,z=eval(parse(text=paste("dvals",i,"[,,1,1,1,1]",sep=""))))
      } else {
	sfc = list(x=lon-360,y=lat,z=eval(parse(text=paste("dvals",i,"[,,1,1,1,",seasonidx,"]",sep=""))))
      }    

      png(paste(outfilepath,outfilename,sep=""),width=600,height=480,units="px",pointsize=12)
      surface(sfc,type="I",main=ptitle,xlab=xlab,ylab=ylab,zlim=colorbar_info[[1]],breaks=colorbar_info[[2]],col=colorbar_info[[3]])#
      map("world",add=TRUE)
      mtext(extratext,side=1,line=-1.5)
      dev.off()
	message(outfilename," Written!")
      
    }
  }
 message("Finished making plots!")
}


#########################
# Helper function

colorramp = function(inputdata,colorchoice,Blimit,type = "difference",use_fixed_scale = FALSE, fixed_scale=c(-100,100)){
  
  if(use_fixed_scale==FALSE){
  message("Not using a fixed scale")
  datarange = range(inputdata,na.rm=TRUE)
  datarange[1]=floor(datarange[1])
  if(datarange[1] %% 2 != 0) datarange[1]=datarange[1]-1
  datarange[2]=ceiling(datarange[2])
  if(datarange[2] %% 2 != 0) datarange[2]=datarange[2]+1
  } else {
  message("Using a fixed scale")
  #tmp = strsplit(fixed_scale,",")
  #datarange = c(as.numeric(tmp[[1]][1]),as.numeric(tmp[[1]][2]))
  datarange = fixed_scale
  }

  if(datarange[1]>=0 & type=="difference"){centerpoint = 0; startpoint=0; datarange[1]=0; message("type=difference")}
  if(datarange[1]<0 & type=="difference"){centerpoint = 0; startpoint=datarange[1]; message("type=difference"); if(datarange[2]<0){datarange[2]=0}}
  
  if(type=="ratio"){centerpoint=1; startpoint=datarange[1]; message("type=ratio")}
  if(type=="raw"){centerpoint=datarange[1]; startpoint=datarange[1]; message("type=raw")}
  
  breakcheck = 1
  breaklist = c(0.001,0.002,0.0025,0.005,0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.3,0.5,1,2,3,4,5,10,20,25,30,50,100,200,250,300,500,1000)
  
  actualbins = diff(datarange)/breaklist
  actidx = which(actualbins<Blimit)
  dataact = actualbins[actidx]-floor(actualbins[actidx])
  
  if(any(dataact==0)==TRUE){
    message("exact match for bins")
    dataidx = which(dataact==0)
    breakcheck=actidx[dataidx[1]]
  } else {
    message("no exact match going through while loop")
    checkpoint = any(dataact==0)
    counter=1
    while(checkpoint==FALSE){
      datarange[1] = floor(datarange[1]/(10^counter))*10^counter
      datarange[2] = ceiling(datarange[2]/(10^counter))*10^counter
      actualbins = diff(datarange)/breaklist
      actidx = which(actualbins<Blimit)
      dataact = actualbins[actidx]-floor(actualbins[actidx])
      dataidx = which(dataact==0)
      
      if(length(dataidx)>=1){
        breakcheck=actidx[dataidx[1]]
        checkpoint = any(dataact==0)
        break
      } else {
        counter=counter+1
        checkpoint = any(dataact==0)
      }
      
    }
  }
  
  if(datarange[2]==0 & colorchoice=="redtoblue") colorchoice="redtowhite"
  if(startpoint==0 & centerpoint==0 & colorchoice=="bluetored") colorchoice="whitetored"
  if(startpoint==0 & centerpoint==0 & colorchoice=="browntogreen") colorchoice="whitetogreen"
  
  zlimdiff = datarange
  breaksdiff = c(seq(datarange[1],datarange[2],by=breaklist[breakcheck]))
  
  if(any(breaksdiff==centerpoint)==FALSE & zlimdiff[1]<centerpoint & zlimdiff[2]>centerpoint){
    idx = which(abs(breaksdiff)==min(abs(breaksdiff)))
    if(length(idx)==1) breaksdiff[idx]=centerpoint
    if(length(idx)>1){
      breaksdiff = c(breaksdiff[1:(idx[1]-1)],centerpoint,breaksdiff[(idx[2]+1):length(breaksdiff)])
    }
  } 

message("zlimdiff: ",zlimdiff)
message("centerpoint: ",centerpoint)
message("startpoint: ",startpoint)
  
  if(startpoint==centerpoint){
message("startpoint matches centerpoint")
    if(colorchoice == "whitetored") colorbardiff = colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-1)
    if(colorchoice == "yellowtored") colorbardiff = colorRampPalette(c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"))(length(breaksdiff)-1)
    
    if(colorchoice == "whitetogreen") colorbardiff = colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-1)
  } else {
    if(datarange[2]>centerpoint){
     message("datarange[2] > centerpoint")
     message("colorchoice = ",colorchoice)
      zeroidx = which(breaksdiff==centerpoint)
      if(colorchoice == "bluetored"){
        colorbardiff = c(colorRampPalette(c("#053061","#2166ac","#4393c3","#92c5de","#d1e5f0","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#fddbc7","#f4a582","#d6604d","#b2182b","#67001f"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "redtoblue"){
        colorbardiff = c(colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "browntogreen"){
        colorbardiff = c(colorRampPalette(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))(length(breaksdiff)-zeroidx))
      } 
      if(colorchoice == "greentobrown"){
        colorbardiff = c(colorRampPalette(c("#003c30","#01665e","#35978f","#80cdc1","#c7eae5","#f5f5f5"))(zeroidx-1),colorRampPalette(c("#f5f5f5","#f6e8c3","#dfc27d","#bf812d","#8c510a","#543005"))(length(breaksdiff)-zeroidx))
      } 
      
    } else {
      message("odd if")
      if(colorchoice == "redtowhite"){
        colorbardiff = colorRampPalette(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f5f5f5"))(length(breaksdiff)-1)
      } 
      
    }
  }
  
  output = list(zlimdiff,breaksdiff,colorbardiff)
  output
  
}

################################################################################

base.path <- Sys.getenv("S8_PATH")
if(base.path==""){
  #If S8_tools not set, assume that you are working from the correct directory
  base.path <- getwd()
}
#Note: if you use this strucutre, it is VITAL htat you do not include ANY scripts in the directories that
#you are sourcing. Those scripts will run, and you will be frustrated.
out <- source(paste0(base.path,'/s8_utils/LoadRPackages.R'))
out <- sapply(list.files(pattern="[.]R$", path=paste(base.path,'/FudgeIO_reboot/',sep=''), full.names=TRUE), source);
out <- sapply(list.files(pattern="[.]R$", path=paste0(base.path,'/s8_utils/'), full.names=TRUE), source, .GlobalEnv);
#source("/home/cew/Code/section_8_tools/s8_utils/gfdl_bplt.R")
library(graphics)
library(fields)
library(maps)

###############################################################################

###Main method of the script
args <- as.list(commandArgs(trailingOnly = TRUE))
ParseSnakeArgs <- function(arg.list){
  option_list <- list(
    #Input and output options first: input, output, timestamp, overlay
    make_option(c("-i", "--input"), action="store", default='na',
                dest='input',
                help=paste("Input to the script in the form of two or more files (-i 'input1.nc, input2.nc');", 
                                "all included in one list, but remember more files take more power to run")),
    make_option(c("-o", "--out_base"), action="store", default='na',
                help=paste("Base for the image outfile with no extension.", 
                           "There are pngs produced, if outbase is provided, ",
                           "then numbers are appended to outbase for instances, ",
                           "with multiple input files. Should not be used with ",
                           "-f")),
    make_option(c("-p","--out_path"), action="store", default='',
                help=paste("Path for writing the output images, ", 
                           "if no value given, will write in working directory.")),
    make_option(c("-f","--file_base"), action='store', default=TRUE, 
                help=paste("Should the filenames of the input files be used in the output filenames?", 
                           "If TRUE the function will use the input filenames to create the output filenames. ",
                           "If FALSE, -o is required. ",
                           "Will throw an error if -f is TRUE and -o is provided.")),
    make_option(c("-b",'--bin_limit'), action='store', default=20,
                help=paste("Maximum number of bins for the contour plots.", 
                           ' Must be a minimum of 5')),
    make_option(c("--units"), action='store', default='na',
                help=paste("The units in which you want the output presented,", 
                           "presented in a udunits-parseable form.", 
                           "(For more information on accepted units, please consult:", 
                           "cat /home/esd/Code/section_8_tools/udunits_list.txt)",
                           "If not specified, defaults to the units of the input file. ",
                           "If --units is not provided, graphics are made in the input units. ",
                           "If --units is provided, the data are converted to the provided units")),
    make_option(c("--spatmask"), action="store", default='na',
                help=paste("Whether to mask the input data by one or more spatial masks.", 
                           "Is applied to all input and will throw an error if dims disagree.")),
    make_option(c("-v","--var"), action="store", default='',
                help=paste("Variable to be plotted. Must be the same between all files, ", 
                           "otherwise an error will be thrown")),
    make_option(c("-s","--stat_used"), action="store", default='mae',
                help=paste("Statistic which will be plotted. If not", 
                           "set, there will be an error")), 
    make_option(c("-c","--color_choice"), action="store", default='yellowtored',
                help=paste("Color ramp choice for the contours being plotted", 
                           "Defaults to 'yellowtored', other options include, ",
                           "'whitetored', 'whitetogreen', 'bluetored', 'redtoblue', ",
                           "'browntogreen', and 'greentobrown'. All are color-blind friendly.")),
    make_option(c("-x","--use_fixed_scale"), action="store", default=FALSE,
                help=paste("Should there be a fixed scale used?", 
                           "Defaults to FALSE (dynamic scale based on input data). ",
                           "Will throw an error if --use_fixed_scale is TRUE, ",
                           "and no fixed scale is provided.")),
    make_option(c("-g","--fixed_scale"), action="store", default='-100,100',
                help=paste("Fixed scale for plotting, listed as two comma separated values min and max (-g -100,100). ", 
                           "Will throw an error if --use_fixed_scale is TRUE, ",
                           "and no fixed scale is provided.")),
    make_option(c("-t","--type"), action="store", default="raw",
                help=paste("What type of data are you providing? This changes how the color scaling responds a bit. Options include 'raw','difference','ratio'. ", 
                           "For option 'raw' (default)- this is for things like mae and raw variable data. ",
                           "For option 'difference'- this is for things like bias or projected change. ",
			   "For option 'ratio'- this is for things like ratios of mae, or similar stats.")),
    make_option(c("-z","--useseason"), action="store", default="FALSE",
                help=paste("Are you using something where a seasonmask was applied? FALSE - no season mask, TRUE - season mask used. ", 
                           "Will break if -useseason is TRUE and -season isn't provided. ")),
make_option(c("-S","--season"), action="store", default="DJF",
                help=paste("If -useseason = TRUE, which season do you want? ", 
                           "Currently uses only 'DJF', 'MAM', 'JJA', and 'SON' as options."))
  )
  description = paste('Given one or more input netcdf files with results from CalcTimeDiffStats1,',
                      "creates contour plots for a variable of choice. ",
                      "The internal function will create a color ramp which matches for all ",
                      "plots if multiple files are used.")
  epilouge = "Please note: flags may be specified in any order, and '=' not required to specify strings."
  arg.list.mod <- allow.neg.commandline.opts(arg.list, c('--ylim'))
  usage = paste("usage: %prog -i input -o out_base -p out_filepath -v var", 
                "[-f TRUE] [--stat_used mae]", 
                "[--bin_limit 20] [--units mm/day]",
                "[--stat_used statistic] [--spatmask] [--color_choice color_choice]",
                "[--used_fixed_scale FALSE] [--fixed_scale -100,100] [--type raw] [--verbose] [-h --help]")
  method.parser <- OptionParser(option_list=option_list, usage=usage, 
                                description = description, epilogue=epilouge)
  return(parse_args(object=method.parser, unlist(arg.list.mod)))
}

if (length(args) > 0){
  input.params <- ParseSnakeArgs(args)  
  input.params$plot.file <- T  
}else{
  input.params <- list()
  input.params$out_base="testgraphic"
  input.params$input = "/home/Adrienne.Wootten/testing/outputs/prannstats_day_PMgCRprp1-BCQM-A00r1e2X01K00_amip_r1to2i1p1_US48_19790101-20081231_ensmbl_0mask.nc,/home/Adrienne.Wootten/testing/outputs/prannstats_day_PMgCRprp1-BCQM-A00r2e2X01K00_amip_r1to2i1p1_US48_19790101-20081231_ensmbl_0mask.nc,/home/Adrienne.Wootten/testing/outputs/prannstats_day_PMgCRprp1-BCQM-A01r1e3X01K00_sst2090_r1to3i1p1_US48_20860101-20951231_ensmbl_0mask.nc,/home/Adrienne.Wootten/testing/outputs/prannstats_day_PMgCRprp1-BCQM-A02r1e3X01K00_sst2090_r1to3i1p2_US48_20860101-20951231_ensmbl_0mask.nc"
  input.params$bin_limit = 20
  input.params$file_base = FALSE
  input.params$out_path = "/home/Adrienne.Wootten/testing/graphics/"
  input.params$spatmask = "/net3/kd/PROJECTS/DOWNSCALING/DATA/MASKS/geomasks/US48/CONUS/CONUS_masks.nc"
  input.params$color_choice = "yellowtored" 
  input.params$var = "pr_jan2dec_count_mae"
  input.params$stat_used = 'mae'
  input.params$units ='mm/day'
  input.params$use_fixed_scale = FALSE
  input.params$fixed_scale = 'na'
  input.params$type="raw"
  input.params$useseason=FALSE
  input.params$season = 'na'
}

#############
###Add method-specific checks here



#Check on presence of var
if(input.params$fixed_scale=='na' & input.params$use_fixed_scale==TRUE){
  stop(paste("Error in MakeContouredMap.R: --fixed_scale not specified and --use_fixed_scale is TRUE"))
}

#print(parsed.args)
if(input.params$fixed_scale!='na'){
fscale = unlist(strsplit(input.params$fixed_scale,",",fixed=TRUE))

if(1 %in% grep("n",fscale)){
  fscale[1] = paste("-",substr(fscale[1],2,nchar(fscale[1])),sep="")
}
if(2 %in% grep("n",fscale)){
  fscale[2] = paste("-",substr(fscale[2],2,nchar(fscale[2])),sep="")
}


fscale = as.numeric(fscale)
} else {
fscale = input.params$fixed_scale
}

###########################
# parsing inputs

#Parse the inputs
dset.vec <- unlist(strsplit(input.params$input, ","))
print(input.params$var)
###########################

#print(input.params)

MakeContouredMap(inputfiles=dset.vec,stat_used=input.params$stat_used,var=input.params$var,
                 outfile_by_filename=input.params$file_base,outfilepath=input.params$out_path,
                 colorchoice=input.params$color_choice, outfile_base=input.params$out_base, 
                 bin_limit=input.params$bin_limit, units=input.params$units,spatmask.file=input.params$spatmask,
		 fixed_scale=fscale,use_fixed_scale=input.params$use_fixed_scale,type=input.params$type,useseason = input.params$useseason,season=input.params$season)

#MakeContouredMap <- function(inputfiles,stat,var,colorchoice,units,bin_limit,outfile_base="testgraphic",outfile_by_filename=FALSE){
  

#gcms = c("MIROC5","MRI-CGCM3","MPI-ESM-LR","MPI-ESM-MR")
gcms = "MRI-CGCM3"
years = 1950:2005

#tasmax_bcc-csm1-1-m_r1i1p1_historical_1961.nc
#tasmax_bcc-csm1-1-1-m_r1i1p1_historical_1961.nc

for(g in 1:length(gcms)){
for(y in 1:length(years)){
  command = paste("ncks -O -h --mk_rec_dmn time tasmax_",gcms[g],"_r1i1p1_historical_",years[y],".nc tasmax_",gcms[g],"_r1i1p1_historical_",years[y],".nc",sep="")
  system(command,wait=TRUE)
  message("Finished conversion for year ",years[y])
}
  message("Finished conversation for ",gcms[g])
}

#####

#ncrcat tasmax_bcc-csm1-1-m_r1i1p1_historical_1950.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1951.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1952.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1953.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1954.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1955.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1956.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1957.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1958.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1959.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1960.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1961.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1962.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1963.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1964.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1965.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1966.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1967.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1968.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1969.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1970.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1971.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1972.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1973.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1974.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1975.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1976.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1977.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1978.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1979.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1980.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1981.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1982.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1983.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1984.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1985.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1986.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1987.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1988.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1989.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1990.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1991.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1992.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1993.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1994.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1995.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1996.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1997.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1998.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_1999.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2000.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2001.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2002.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2003.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2004.nc tasmax_bcc-csm1-1-m_r1i1p1_historical_2005.nc tasmax_bcc-csm1-1-m_r1i1p1_historical.nc

gcms = "MRI-CGCM3"
years = 1950:2005


commandstart ="ncrcat "

for(g in 1:length(gcms)){
  for(y in 1:length(years)){
    if(y==1){
      command = paste(commandstart,"tasmax_",gcms[g],"_r1i1p1_historical_",years[y],".nc ",sep="")
    } else {
      command = paste(command,"tasmax_",gcms[g],"_r1i1p1_historical_",years[y],".nc ",sep="")
    }
  }
  command = paste(command,"tasmax_",gcms[g],"_r1i1p1_historical.nc",sep="")
  system(command,wait=TRUE)
  message("Finished concatenation for ",gcms[g])
}






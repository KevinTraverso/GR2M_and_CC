'EXTRACCION DE LA INFORMACION HISTORICA DE MODELOS GCM 
PARA SU USO EN CAMBIO CLIMATICO
CASO DE APLICACIÓN: CUENCA VILCANOTA

@autor: Kevin Arnold Traverso
mail: arnold.traverso@gmail.com
verificado: Waldo Lavado
mail: waldo.lavado@gmail.com'

# CARGAR/INSTALAR LOS PAQUETES NECESARIOS ---------------------------------

if(!require(raster)) install.packages("raster")
if(!require(ncdf4))  install.packages("ncdf4")
if(!require(rgdal))  install.packages("rgdal")
if(!require(zoo))    install.packages("zoo")
if(!require(dplyr))  install.packages("dplyr")

# wd <- "D:/2019/2_SENAMHI/Application_GR2M_and_CC_Kevin/"; setwd(wd); getwd() # Definir Directorio
rm(list = ls())

# INGRESAR EL MODELO Y LA CUENCA ------------------------------------------

shp <- shapefile("./Shp_basin/Basin.shp")
# cARPETA DONDE ESTAN LOS MODELOS HISTORICOS
files <- list.files(path = './GCM_model/Historical/', 
                    pattern = '\\.nc$',
                    full.names = TRUE)

# FUNCIONES (cargar)
generate_mask_geom <- function(cov,geom){
  specialcov = cov
  specialcov[] = 1:ncell(cov)
  Position_rowcol <- function(i){
    quad1 <- unlist(raster::extract(specialcov, geom[i,], small=TRUE))
  }
  position<-lapply(1:length(geom), Position_rowcol)
  return(position)
}
mask_fast_extract <- function(cov, positionP, fun = mean, na.rm = TRUE){ 
  matrix_R <- t(as.matrix(cov)) 
  sapply(1:length(positionP), function(i){
    Value <- matrix_R[positionP[[i]]]
    fun(Value, na.rm=TRUE)
  })
}

# EXTRAER LOS DATOS HISTORICOS ADAPTADOS PARA PISCO ----------------------

for (i in 1:length(files)) {
  file <- files[i]
  ncin <- nc_open(file)
  print(ncin)
  institution <- ncatt_get(ncin, 0, "model_id")
  mod<-institution$value
  expm <- ncatt_get(ncin, 0, "experiment_id")
  expmn <- expm$value
  prc <- ncin$var$pr$name
  tas <- ncin$var$tas$name
  pr<-raster(file)
  pr<-brick(file)
  pr<-rotate(pr)
  res(pr)
  ex<-extent(x = c(-90,-62), y=c(-20,-0))
  ext_r<-raster(ex,crs='+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
  pr_z<-crop(pr,ext_r)
  gmask_gcm <- generate_mask_geom(cov = pr_z[[1]], geom = shp)
  serie <- mapply(function(i) mask_fast_extract(pr_z[[i]], gmask_gcm),
                     i=1:nlayers(pr_z))
  serie <- serie*30
  Date.Ini = "1950-01-01"
  Date.Fin = "2005-12-01"
  Date <- seq(as.Date(Date.Ini), as.Date(Date.Fin), by="month")
  serie2 <- zoo(serie, order.by = Date)
  serie2 <- window(x = serie2, start = "1981-01-01", end = "2005-12-01")
  # Guardar en Rdata
  save(serie2, file = paste0(mod,'_',expmn,'_',prc,tas,'.Rdata'))
  # Guardar en csv
  # write.csv(serie2, paste0(mod, '_',expmn, '_', prc, tas ,'.csv'))
  print("Full processing !! the files have been generated !!")
}

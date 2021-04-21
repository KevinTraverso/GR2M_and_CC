'EXTRAER INFORMACION DE TEMPERATURAS

@autor: Kevin Arnold Traverso
mail: tkevinarnold@gmail.com
Verificado: Waldo Lavado
mail: waldo.lavado@gmail.com'

# CARGAR/INSTALAR LOS PAQUETES NECESARIOS ---------------------------------

if(!require(ncdf4)) install.packages("ncdf4")
if(!require(raster)) install.packages("raster")
if(!require(zoo)) install.packages("zoo")
if(!require(lubridate)) install.packages("lubridate")
if(!require(dplyr)) install.packages("dplyr")

# wd <- "D:/2019/2_SENAMHI/Application_GR2M_and_CC_Kevin/"; setwd(wd); getwd() # Definir Directorio
rm(list = ls()) # Liberar Enviroment

# INGRESO DE VARIABLES ----------------------------------------------------

# Ingresar el shape de la cuenca (la proyeccion debe estar en UTM)
shp <- shapefile("./Shp_basin/Basin.shp")
summary(shp) # Informacion del Shapefile
# Cargar la informacion de PISCO precipitaciones y evapotranspiracion
PISCO_tx <- raster::brick("./PISCO-TEM/PISCOmtx_v1.1.nc")
PISCO_tn <- raster::brick("./PISCO-TEM/PISCOmtn_v1.1.nc")

# RECOPIACION DE INFORMACIÓN ----------------------------------------------

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
    fun(Value,na.rm=T)
  })
}
gmask_pisco <- generate_mask_geom(cov = PISCO_tx[[1]], geom = shp)
gmask_pisco <- generate_mask_geom(cov = PISCO_tn[[1]], geom = shp)
# Extraer los datos de precipitacion mensual
serie_tx <- mapply(function(i) mask_fast_extract(PISCO_tx[[i]], gmask_pisco),
                   i=1:nlayers(PISCO_tx))
serie_tn <- mapply(function(i) mask_fast_extract(PISCO_tn[[i]], gmask_pisco),
                   i=1:nlayers(PISCO_tn))
# promediando los datos
serie_ts <- rowMeans(cbind(serie_tx, serie_tn))
# Guardando los datos
Date_ini = as.Date("1981-01-01")
Date_fin = as.Date("2016-12-31")
BasinTEMP <- data.frame("DateR" = as.POSIXlt(seq(from = Date_ini,
                                                to = Date_fin,
                                                by = "month"),
                                            format = "%Y-%m-%d"), 
                       "Tx" = as.numeric(serie_tx),
                       "Tn" = as.numeric(serie_tn),
                       "Ts" = as.numeric(serie_ts))
BasinTEMP <- BasinTEMP[1:300,] # Guardando solo hasta 2005
save(BasinTEMP, file = "BasinTEMP.Rdata")

'SIMULACIÓN DEL MODELO GR2M A PARTIR DE LA INFORMACIÓN PISCO
CASO DE APLICACIÓN: CUENCA VILCANOTA

@autor: Kevin Arnold Traverso
mail: tkevinarnold@gmail.com
Verificado: Waldo Lavado
mail: waldo.lavado@gmail.com'

# CARGAR/INSTALAR LOS PAQUETES NECESARIOS ---------------------------------

if(!require(ncdf4)) install.packages("ncdf4")
if(!require(raster)) install.packages("raster")
if(!require(airGR)) install.packages("airGR")
if(!require(zoo)) install.packages("zoo")
if(!require(hydroTSM)) install.packages("hydroTSM")
if(!require(lubridate)) install.packages("lubridate")
if(!require(dplyr)) install.packages("dplyr")

# wd <- "D:/2019/2_SENAMHI/Application_GR2M_and_CC_Kevin/"; setwd(wd); getwd() # Definir Directorio
rm(list = ls()) # Liberar Enviroment

# INGRESO DE VARIABLES ----------------------------------------------------

# Ingresar el shape de la cuenca (la proyeccion debe estar en UTM)
shp <- shapefile("./Shp_basin/Basin.shp")
summary(shp) # Informacion del Shapefile
# Cargar la informacion de PISCO precipitaciones y evapotranspiracion
PISCO_pp <- raster::brick("./PISCO-PP/PISCOpm.nc")
PISCO_et <- raster::brick("./PISCO-ETO/PISCOmpe_oudin_v1.1.nc")
# visualizar Pisco
spplot(PISCO_et$X1981.01.16, main = "Evapotraspiración")
spplot(PISCO_pp$X0, main = "Precipitación")

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
gmask_pisco <- generate_mask_geom(cov = PISCO_pp[[1]], geom = shp)
# Extraer los datos de precipitacion mensual
serie_pp <- mapply(function(i) mask_fast_extract(PISCO_pp[[i]], gmask_pisco),
                   i=1:nlayers(PISCO_pp))
plot(serie_pp, type = "l", col = "blue",
     main = "Serie de Precipitación - PISCOp")
# Extraer los datos de Evapotranspiracion
serie_et <- mapply(function(i) mask_fast_extract(PISCO_et[[i]], gmask_pisco),
                   i=1:nlayers(PISCO_et))
plot(serie_et, type = "l", col = "red",
     main = "Serie de Evapotranspiracion - PISCOe")
# Ingresar la data de caudales [m3/s]
Qobs_d <- read.zoo("./Caudales/Caudales_km105.csv", header = TRUE, sep = ",")
plot(Qobs_d, type= "l", col = "blue", main = "Serie de Caudales diarios")
# Transformar los caudales a mensuales
Qobs_m <- daily2monthly(x = Qobs_d, FUN = mean, na.rm = TRUE)
plot(Qobs_m, type = "l", col = "blue", main = "Serie de caudales menusuales")
# Transformar los cairales de m3/s a mm
Q_mm <- Qobs_m*days_in_month(Qobs_m)*24*3600/(area(shp)/1000)
plot(Q_mm, col = "blue", main = "Serie de caudales en mm")

# INFORMACION DE LA CUENCA ------------------------------------------------

Date_ini = as.Date("1981-01-01")
Date_fin = as.Date("2016-12-31")
BasinObs <- data.frame("DateR" = as.POSIXlt(seq(from = Date_ini,
                                                     to = Date_fin,
                                                     by = "month"),
                                            format = "%Y-%m-%d"), 
                       "P" = as.numeric(serie_pp), 
                       "E" = as.numeric(serie_et),
                       "Qm3s" = as.numeric(Qobs_m),
                       "Qmm" = as.numeric(Q_mm))
# Guargar la data generada de la cuenca (1981-2016)
save(BasinObs, file = "BasinOBS.Rdata")

# Aplicacion para el cambio climatico
# para la aplicacion de cambio climatico primero se hará
# la seleccion de los datos de la cuenca hasta diciembre del 2005
# necesarios para el analisis de downscaling
BasinCC <- BasinObs[1:300,] # Seleccion de datos hasta el 2005
save(BasinCC, file = "BasinCC.Rdata") # Grardando nueva base

# INGRESO DE DATOS AL MODELO GR2M -----------------------------------------

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M, 
                                 DatesR = BasinObs$DateR,
                                 Precip = BasinObs$P, 
                                 PotEvap = BasinObs$E)
head(BasinObs)
# Definir el periodo de Calentamiento
Ind_Wup <- seq(which(format(BasinObs$DateR, format = "%Y-%m")=="1981-01"),
               which(format(BasinObs$DateR, format = "%Y-%m")=="1983-12"))
# Definir el periodo de Calibracion
Ind_Run <- seq(which(format(BasinObs$DateR, format = "%Y-%m")=="1984-01"),
               which(format(BasinObs$DateR, format = "%Y-%m")=="2011-12"))
# Definir el periodo de Validacion
Ind_Runv <- seq(which(format(BasinObs$DateR, format = "%Y-%m")=="2011-01"),
                which(format(BasinObs$DateR, format = "%Y-%m")=="2016-12"))

# PRIMERA SIMULACION ------------------------------------------------------

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M, 
                               InputsModel = InputsModel,
                               IndPeriod_WarmUp = Ind_Wup,
                               IndPeriod_Run = Ind_Run)
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, 
                               InputsModel = InputsModel,
                               RunOptions = RunOptions, 
                               Obs = BasinObs$Qmm[Ind_Run])
Param <- c(265.072, 1.040) # Parametros aleatorios(base)
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, 
                              RunOptions = RunOptions, Param = Param)
plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run]) # Ver Resultados
OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

# CALIBRACION -------------------------------------------------------------

CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR2M, 
                                   FUN_CALIB = Calibration_Michel)
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
                                   RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, 
                                   CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR2M)
Param <- OutputsCalib$ParamFinalR # Transferencia de nuevos parametros
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
                              RunOptions = RunOptions,
                              Param = Param)
plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run])
OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

# Guardando los parametros encontrados del modelo
save(Param, file = "Parametros_Calibrados.Rdata")

# CALIBRACION POR CRITERIO DE EFICIENCIA ----------------------------------

## Criterio de Eficiencia: Nash-Sutcliffe Efficiency
# InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
#                                RunOptions = RunOptions, Obs = BasinObs$Qmm[Ind_Run])
# OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)
# # Criterio de Eficiencia: Klling-Gupta Efficiency
# InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_KGE, InputsModel = InputsModel,
#                                RunOptions = RunOptions, Obs = BasinObs$Qmm[Ind_Run])
# OutputsCrit <- ErrorCrit_KGE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

# Ingresar el criterio
## ErrorCrit_KGE, ErrorCrit_KGE2, ErrorCrit_NSE, ErrorCrit_RMSE
# OutputsCalib <- Calibration_Michel(InputsModel = InputsModel,
#                                    RunOptions = RunOptions,
#                                    InputsCrit = InputsCrit, 
#                                    CalibOptions = CalibOptions,
#                                    FUN_MOD = RunModel_GR2M, 
#                                    FUN_CRIT = ErrorCrit_KGE) # Cambiar criterio
# Param <- OutputsCalib$ParamFinalR # Ingreso de nuevos Parametros
# OutputsModel <- RunModel_GR2M(InputsModel = InputsModel,
#                               RunOptions = RunOptions,
#                               Param = Param)
# plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Run])
# OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

# VALIDACIÓN --------------------------------------------------------------

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                               InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Runv)
Param <- OutputsCalib$ParamFinalR # Ingreso de Parametros calibrados
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, 
                              RunOptions = RunOptions, Param = Param)
plot(OutputsModel, Qobs = BasinObs$Qmm[Ind_Runv])
OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = OutputsModel)

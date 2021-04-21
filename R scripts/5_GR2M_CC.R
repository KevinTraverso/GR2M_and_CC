'SIMULACIÓN DEL MODELO GR2M A PARTIR DE LA INFORMACIÓN GCM CAMBIO CLIMATICO
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
centroide <- centroides <- rgeos::gCentroid(shp, byid = TRUE)
centroide2 <- spTransform(centroide, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cent <- centroide2
cent # Verificar el dato de la latitud
# INFORMACION DE LA CUENCA ------------------------------------------------

load(file = "GCM_DOWNSCALED_PR.Rdata")
precip_cc <- gcmdw
plot(precip_cc)
load(file = "GCM_DOWNSCALED_TS.Rdata")
temp_cc <- gcmdw

Date_ini = as.Date("1981-01-15")
Date_fin = as.Date("2100-12-15")
BasinCC <- data.frame("DateR" = as.POSIXlt(seq(from = Date_ini,
                                               to = Date_fin,
                                               by = "month"),
                                           format = "%Y-%m-%d"), 
                      "P_CC" = as.numeric(precip_cc), 
                      "T_CC" = as.numeric(temp_cc))

JD <- seq.Date(Date_ini, Date_fin, by = "month")
JD <- format(JD, "%j")
JD <- as.numeric(JD)
# Calculo de la Evapotranspiracion Por Oudin
PE_CC <- PE_Oudin(JD, BasinCC$T_CC, -13.82, # Ingresar el dato de Latitud
                  LatUnit = "deg",
         TimeStepIn = "daily", TimeStepOut = "daily")
# Actualizamos la base de la cuenca en CC
BasinCC <- data.frame(BasinCC, "PE_CC" = PE_CC)
head(BasinCC)
# Guargar la data generada de la cuenca (1981-2016)
save(BasinCC, file = "BasinCC.Rdata")

# Aplicacion para el cambio climatico
# para la aplicacion de cambio climatico primero se hará
# la seleccion de los datos de la cuenca desde 2006 hasta 2100
# INGRESO DE DATOS AL MODELO GR2M -----------------------------------------

InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M, 
                                 DatesR = BasinCC$DateR,
                                 Precip = BasinCC$P_CC, 
                                 PotEvap = BasinCC$PE_CC)
# Definir el periodo de Simulacion
Ind_Run <- seq(which(format(BasinCC$DateR, format = "%Y-%m")=="2006-01"),
               which(format(BasinCC$DateR, format = "%Y-%m")=="2100-12"))

# SIMULACION --------------------------------------------------------------

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M, 
                               InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run)

load(file = "Parametros_Calibrados.Rdata") # Parametros obtenidos en la calibracion
Param <- Param# Ingreso de Parametros calibrados
OutputsModel <- RunModel_GR2M(InputsModel = InputsModel, 
                              RunOptions = RunOptions, 
                              Param = Param)
# Visualizando el modelamiento GR2M aplicando Cambio Climatico
plot(OutputsModel)

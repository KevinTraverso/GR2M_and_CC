'DOWNSCALLING POR LA METODOLOGIA DE QUANTILE MAPPING

Revisar el articulo: Technical Note: Downscaling RCM precipitation
to the station scale using statistical transformations - a comparison of methods
L. Gudmundsson et al. (2012)
La correccion se realiza cada tres meses, el codio esta realizado para la 
informacion a nivel mensual

Adaptado por:
@autor: Kevin Arnold Traverso
mail: tkevinarnold@gmail.com
Verificado: Waldo Lavado
mail: waldo.lavado@gmail.com'

# CARGAR/INSTALAR LOS PAQUETES NECESARIOS ---------------------------------

if(!require(qmap)) install.packages("qmap")
if(!require(zoo)) install.packages("zoo")
if(!require(latticeExtra)) install.packages("latticeExtra")
if(!require(raster)) install.packages("Raster")

# wd <- "D:/2019/2_SENAMHI/Application_GR2M_and_CC_Kevin/"; setwd(wd); getwd() # Definir Directorio
rm(list = ls()) # Liberar Enviroment

# CARGAR LA INFORMACION NECESARIA -----------------------------------------

# Informacion Historica PRECIPITACIONES
load(file = "BasinCC.Rdata")
head(BasinCC)
date <- as.Date(BasinCC[,1])
OBS_hist <- zoo(x = BasinCC[,2],
                order.by = seq(date[1],
                               date[length(date)], 
                               by = "month")) # Precip/Temp
# str(OBS_hist)
# plot(OBS_hist)

# Informacion Historica CGM
load(file = "MIROC5_historical_pr.Rdata") # Precipitacion Historica
GCM_hist <- serie2
# str(GCM_hist)
# plot(GCM_hist)

# Informacion Futura del modelo RCP
load(file = "MIROC5_rcp85_pr.Rdata")
GCM_fut <- serie2
# str(GCM_fut)
# plot(GCM_fut)

GCM_model <- rbind(GCM_hist,GCM_fut)
data_at <- cbind(OBS_hist, GCM_model)
data_wt <- cbind(OBS_hist, GCM_hist)

# ..............................................................................
rm(list = ls()) # Liberar Enviroment
# Informacion Historica TEMPERATURA/ACTIVAR PARA SU USO
load(file = "BasinTEMP.Rdata")
head(BasinTEMP)
date <- as.Date(BasinTEMP[,1])
OBS_hist <- zoo(x = BasinTEMP[,4],
                order.by = seq(date[1],
                               date[length(date)],
                               by = "month")) # Precip/Temp
str(OBS_hist)
# plot(OBS_hist)

# Informacion Historica CGM
load(file = "MIROC5_historical_tas.Rdata") # Precipitacion Historica
GCM_hist <- serie2/30
# str(GCM_hist)
# plot(GCM_hist)

# Informacion Futura del modelo RCP
load(file = "MIROC5_rcp85_tas.Rdata")
GCM_fut <- serie2/30
# str(GCM_fut)
# plot(GCM_fut)

GCM_model <- rbind(GCM_hist,GCM_fut)
data_at <- cbind(OBS_hist, GCM_model)
data_wt <- cbind(OBS_hist, GCM_hist)

# ..............................................................................
# GRAFICANDO LAS SERIES DE TIEMPO
plot(data_at, plot.type = "single", col = c(1, 2), lwd = 1, 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)",
     xlab = "years", );legend("topleft", 
                              legend=c("OBS", "GCM"), col=c(1, 2), lty=1, cex=0.8)

plot(data_wt, plot.type = "single", col = c(1, 2),  lwd = 1, 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)",
     xlab = "years");legend("topleft",
                            legend=c("OBS", "GCM"), col=c(1, 2), lty=1, cex=0.8)

plot(data_wt[months(time(data_wt)) %in% c("December","January","February")], 
     plot.type = "single", col = c(1, 2),  lwd = 1, type='p', 
     main = c("OBS vs GCM"), ylab = "pp (mm/mes)",
     xlab = "years");legend("topleft", 
                            legend=c("OBS", "GCM"), col=c(1, 2), pch=1, cex=0.8)

# GRAFICANDO SCATTERPLOT
plot(GCM_hist~OBS_hist, coredata(data_wt), 
     col = c(1,2));legend("topleft", legend=c("OBS", "GCM"), 
                          col=c(1, 2), pch=1, cex=0.8)

plot(GCM_hist~OBS_hist,
     coredata(data_wt[months(time(data_wt)) %in% 
                        c("December", "January", "February")]),
     col = c(1,2)) ;legend("topleft", legend=c("OBS", "GCM"),
                           col=c(1, 2), pch=1, cex=0.8)

# GRAFICANDO ECDF
ecdfplot(~ OBS_hist +  GCM_hist, data = data.frame(data_wt), 
         lwd = 2, col = c(1, 2))

ecdfplot(~ OBS_hist +  GCM_hist,
         data = data.frame(data_wt[months(time(data_wt)) %in% 
                                     c("December", "January", "February")]), 
         lwd = 2, col = c(1, 2))

# ..............................................................................
# SEASONAL QUANTILE EMPIRICAL MAPPING

data_wt$gcm_dowscaled <- data_wt$GCM_hist
seasons_by_year <- list(c("December","January","February"), 
                        c("March","April","May"), 
                        c("June","July","August"),
                        c("September","October","November"))
seasonal_qm_fit_model <- list()

for(i in 1:4) {
  obs_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$OBS_hist
  mod_sl <- data_wt[months(time(data_wt)) %in% seasons_by_year[[i]]]$GCM_hist
  #MODEL, read!: L. Gudmundsson et al. (2012)
  qm_fit <- fitQmapQUANT(obs = coredata(obs_sl),
                         coredata(mod_sl),
                         qstep = 0.01,
                         nboot = 1, 
                         wet.day = TRUE, # Adaptado para temperatuas negativas
                         type = "linear")
  mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), qm_fit, type = "linear")
  data_wt$gcm_dowscaled[ months(time(data_wt)) %in%
                           seasons_by_year[[i]]] <- mod_sl_qmapped
  seasonal_qm_fit_model[[i]] <- qm_fit 
}

# ..............................................................................
# GRAFICANDO
plot(data_wt, plot.type = "single", col = c(1, 2, 4),  lwd = 2, 
     main = c("OBS vs GCM vs GCM downscaled"),
     xlab = "Years", ylab = "pp (mm/mes)")
plot(data_wt[months(time(data_wt)) %in% c("December",
                                          "January",
                                          "February")], 
     plot.type = "single", col = c(1, 2, 4),  lwd = 2, type = "p",
     main = c("OBS vs GCM vs GCM downscaled"), 
     xlab = "Years", ylab = "pp (mm/mes)")
plot(gcm_dowscaled~OBS_hist, coredata(data_wt), col = c(4,1))
plot(gcm_dowscaled~OBS_hist, coredata(data_wt[months(time(data_wt)) %in%
                                                c("December",
                                                  "January",
                                                  "February")]), 
     col = c(4,1))
ecdfplot(~ OBS_hist +  GCM_hist + gcm_dowscaled, data = data.frame(data_wt), 
         lwd = 3, col = c(1, 2, 4))
ecdfplot(~ OBS_hist +  GCM_hist + gcm_dowscaled,
         data = data.frame(data_wt[months(time(data_wt)) %in%
                                     c("December","January","February")]), 
         lwd = 3, col = c(1, 2, 4))
ecdfplot(~ OBS_hist + gcm_dowscaled, data = data.frame(data_wt), 
         lwd = 3, col = c(1, 2, 4))

# ..............................................................................
# INTERPOLANDO LA INFORMACION

data_at$GCM_downscaled <- data_at$GCM_model
for(i in 1:4) {
  mod_sl <- data_at[months(time(data_at)) %in% seasons_by_year[[i]]]$GCM_model
  mod_sl_qmapped <- doQmapQUANT(coredata(mod_sl), seasonal_qm_fit_model[[i]], type = "linear")
  data_at$GCM_downscaled[ months(time(data_at)) %in% seasons_by_year[[i]]] <- mod_sl_qmapped
}
# verificar la nueva tabla y el grafico
View(data_at)
plot(data_at, plot.type = "single", col = c(1, 2, 4), lwd = 1, 
     main = c("OBS vs GCM VS GCM DOWNSCALED"),
     ylab = "pp (mm/mes)", xlab = "years")

# ..............................................................................
# GUARDANDO LOS DATOS AJUSTADOS PARA GCM
gcmdw <- data_at[,3]
# save(gcmdw, file = "GCM_DOWNSCALED_PR.Rdata") # Activara para Precipitacion
save(gcmdw, file = "GCM_DOWNSCALED_TS.Rdata") # Activar para Temperatura
# write.zoo(data_at[,3], file = "GCM_DOWNSCALED_1981_2100.csv", sep = ",")
print("El proceso se ha completado satisfactoriamente!")

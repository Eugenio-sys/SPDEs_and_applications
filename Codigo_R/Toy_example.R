# =============================================================================
# TOY MODEL SPDE — Red Sísmica PRSN
# Intensidad espacial de sismos agregados por celda de 0.1°
# INLA inla.spde2.pcmatern (alpha=2, nu=1 fijo)
# =============================================================================
options(timeout = 100000) 
install.packages("INLA",
                 repos = c(INLA = "https://inla.r-inla-download.org/R/testing"),
                 dependencies = TRUE)
library(tidyverse)
library(ggplot2)
library(viridis)
library(fmesher)
library(inlabru)
library(jsonlite)
library(INLA)
library(rSPDE) #INLA::SPDEhowto	in HELP
library(splancs)
library(brinla)
library(httr)

hola

# 1. LECTURA DE DATOS JSON 
DATA_DIR <- if (dir.exists("datos")) "datos" else "."
ARCHIVOS <- list.files(DATA_DIR, pattern = "\\.json$", full.names = TRUE)
if (!length(ARCHIVOS))
  stop("No encuentro archivos .json en '", normalizePath(DATA_DIR), "'.")

COLS <- c("ID","BULLETIN_NUM","DATE","TIME","LAT","LON","DEPTH","MAGNITUDE",
          "MAG_TYPE","STATIONS","PHASES","FE_REGION","STATUS","NETWORK",
          "INTENSITY","TIMEZONE","ERZ","ERH","RMS","GAP","QUALITY","SOURCE",
          "PRSN_ID")

leer_json <- function(f) {
  x <- jsonlite::fromJSON(f, simplifyDataFrame = TRUE)
  if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  nm <- toupper(trimws(names(x)))
  nm <- gsub("[^A-Z0-9]+", "_", nm)
  nm <- gsub("^_+|_+$",   "", nm)
  names(x) <- nm
  x[] <- lapply(x, as.character)
  falt <- setdiff(COLS, names(x))
  for (cc in falt) x[[cc]] <- NA_character_
  x <- x[, COLS, drop = FALSE]
  x$archivo <- basename(f)
  attr(x, "campos_ausentes") <- falt
  x
}

crudos <- lapply(ARCHIVOS, leer_json)

ESQUEMA <- do.call(rbind, lapply(crudos, function(d) {
  data.frame(
    archivo  = d$archivo[1],
    n        = nrow(d),
    ausentes = paste(attr(d, "campos_ausentes"), collapse = ", "),
    stringsAsFactors = FALSE
  )
}))
ESQUEMA$ausentes[ESQUEMA$ausentes == ""] <- "(ninguno)"

num <- function(x) suppressWarnings(as.numeric(x))

datos <- dplyr::bind_rows(crudos) |>
  dplyr::mutate(
    dt        = as.POSIXct(paste(DATE, TIME), tz = "UTC",
                           format = "%Y-%m-%d %H:%M:%OS"),
    fecha     = as.Date(dt),
    anio      = as.integer(format(dt, "%Y")),
    lat       = num(LAT),
    lon       = num(LON),
    prof      = num(DEPTH),
    mag       = num(MAGNITUDE),
    est       = num(STATIONS),
    fases     = num(PHASES),
    erz       = num(ERZ),
    erh       = num(ERH),
    rms       = num(RMS),
    gap       = num(GAP),
    tipo_mag  = toupper(trimws(MAG_TYPE)),
    calidad   = toupper(trimws(QUALITY)),
    estado    = toupper(trimws(STATUS)),
    fuente    = toupper(trimws(SOURCE)),
    id_tiempo = dplyr::coalesce(PRSN_ID, ID)
  ) |>
  dplyr::filter(!is.na(dt)) |>
  dplyr::arrange(dt)

summary(datos)
str(datos)

# 2. FILTRO Y AGREGACIÓN POR CELDA
# Filtro: bounding box PR + umbral de completitud PRSN (mag >= 2.5)
# Agregación: celdas de 0.1° (~11 km) para tener conteos con varianza real
# Con y_i = 1 por sismo individual el campo SPDE colapsa (sigma ~ 0)

sismos_raw <- datos |>
  dplyr::filter(
    !is.na(lat), !is.na(lon), !is.na(mag),
    lat  >= 17.5,  lat  <= 19.0,
    lon  >= -68.0, lon  <= -64.5,
    mag  >= 2.5,
    prof >= 0
  )

cat("\nSismos filtrados (mag >= 2.5):", nrow(sismos_raw), "\n")

# Agregar en celdas de 0.1 grados
sismos <- sismos_raw |>
  dplyr::mutate(
    lon_bin = round(lon / 0.1) * 0.1,
    lat_bin = round(lat / 0.1) * 0.1
  ) |>
  dplyr::count(lon_bin, lat_bin, name = "n_sismos")

cat("Celdas con al menos 1 sismo:", nrow(sismos),    "\n")
cat("Conteo mínimo:", min(sismos$n_sismos),
    "| máximo:",     max(sismos$n_sismos),
    "| media:",      round(mean(sismos$n_sismos), 1), "\n")

# Coordenadas de los centroides de celda
coords <- as.matrix(sismos[, c("lon_bin", "lat_bin")])

# 3. MESH FEM 
# max.edge interior 0.30° ~ 33 km (PR mide ~170 x 60 km)
# cutoff 0.08° ~ 9 km evita triángulos degenerados en zonas densas (SW)
# offset extiende el dominio para alejar la condición de frontera

mesh <- fmesher::fm_mesh_2d_inla(
  loc      = coords,
  max.edge = c(0.30, 1.0),
  cutoff   = 0.08,
  offset   = c(0.4, 1.5)
)

cat("\nNodos del mesh:  ", mesh$n, "\n")
cat("Triángulos:      ", nrow(mesh$graph$tv), "\n")

# Diagnóstico visual
plot(mesh, main = "Mesh FEM — Puerto Rico (PRSN)")
points(coords, pch = 20, cex = 0.5, col = "red")

# 4. MODELO SPDE 
# SPDE Matérn con PC priors (Simpson et al. 2017)
# alpha = 2 → nu = alpha - d/2 = 2 - 1 = 1 (nu fijo)
# SPDE: (kappa² - Delta) x(s) = W(s)
# PC prior rango: P(rho < 0.5°) = 0.5  →  prior penaliza rangos cortos
# PC prior sigma: P(sigma > 1)  = 0.5  →  prior penaliza varianzas grandes

spde_model <- INLA::inla.spde2.pcmatern(
  mesh        = mesh,
  alpha       = 2,
  prior.range = c(0.5, 0.5),
  prior.sigma = c(1.0, 0.5)
)

# 5. MATRIZ A E ÍNDICES 
# A: n_celdas × n_nodos — pesos baricéntricos
# x(s_i) = A_i · w  (w = campo en nodos del mesh)

A   <- INLA::inla.spde.make.A(mesh = mesh, loc = coords)
idx <- INLA::inla.spde.make.index(name = "field", n.spde = spde_model$n.spde)

cat("Dimensiones de A:", dim(A), "\n")

# 6. INLA STACK 
# Respuesta: n_sismos por celda (conteo Poisson)
# Predictor: log(lambda_i) = mu + x(s_i)
# El offset log(area) = log(0.01 grados^{2}) corrige por el área de cada celda

stk <- INLA::inla.stack(
  data    = list(y = sismos$n_sismos,
                 offset = rep(log(0.01), nrow(sismos))),
  A       = list(A, 1),
  effects = list(
    idx,
    list(intercepto = rep(1, nrow(sismos)))
  ),
  tag = "est"
)

# 7. AJUSTE DEL MODELO
# log(E[y_i]) = offset_i + mu + x(s_i)
# Familia Poisson: sismos son conteos en celdas del espacio
# El campo x(s) captura la heterogeneidad espacial de la intensidad sísmica

formula <- y ~ -1 + intercepto + f(field, model = spde_model)

result <- INLA::inla(
  formula,
  family            = "poisson",
  data              = INLA::inla.stack.data(stk),
  offset            = INLA::inla.stack.data(stk)$offset,
  control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE),
  control.compute   = list(dic  = TRUE, waic = TRUE,
                           mlik = TRUE, config = TRUE),
  verbose           = FALSE
)

# 8. DIAGNÓSTICO DEL AJUSTE 
cat("\n Ajuste\n")
cat("DIC:  ", round(result$dic$dic,   2), "\n")
cat("WAIC: ", round(result$waic$waic, 2), "\n")
cat("Log-verosimilitud marginal:", round(result$mlik[1], 2), "\n")

cat("\n Efectos fijos \n")
print(round(result$summary.fixed[, c("mean","sd","0.025quant","0.975quant")], 4))

cat("\n Hiperparámetros (escala interna INLA) \n")
print(round(result$summary.hyperpar[, c("mean","sd","0.025quant","0.975quant")], 4))

#### Extraer predictor lineal (eta) del resultado
eta_fitted <- result$summary.linear.predictor
head(eta_fitted)



# 9. PARÁMETROS DEL CAMPO SPDE (escala original) 
# inla.spde2.result() extrae rho y sigma en escala original
# do.transf = TRUE aplica la transformación inversa desde la escala interna

spde_result <- INLA::inla.spde2.result(
  inla      = result,
  name      = "field",
  spde      = spde_model,
  do.transf = TRUE
)

# Rango espacial rho (en grados, ~111 km/grado)
rho_mean <- INLA::inla.emarginal(function(x) x,
                                 spde_result$marginals.range.nominal[[1]])
rho_q    <- INLA::inla.qmarginal(c(0.025, 0.975),
                                 spde_result$marginals.range.nominal[[1]])

# Varianza marginal sigma² del campo
sig2_mean <- INLA::inla.emarginal(function(x) x,
                                  spde_result$marginals.variance.nominal[[1]])
sig2_q    <- INLA::inla.qmarginal(c(0.025, 0.975),
                                  spde_result$marginals.variance.nominal[[1]])

cat("\n Parámetros del campo SPDE \n")
cat(sprintf("rho   (rango):  media = %.3f°   IC95%% [%.3f, %.3f]\n",
            rho_mean, rho_q[1], rho_q[2]))
cat(sprintf("sigma (desv.):  media = %.3f    IC95%% [%.3f, %.3f]\n",
            sqrt(sig2_mean), sqrt(sig2_q[1]), sqrt(sig2_q[2])))

# Marginales posteriores de rho y sigma^{2}
par(mfrow = c(1, 2))
plot(spde_result$marginals.range.nominal[[1]], type = "l",
     main = "Rango espacial rho (grados)",
     xlab = "rho", ylab = "densidad posterior")
abline(v = rho_mean, col = "red", lty = 2)
abline(v = rho_q,    col = "red", lty = 3)

plot(spde_result$marginals.variance.nominal[[1]], type = "l",
     main = "Varianza marginal sigma²",
     xlab = "sigma²", ylab = "densidad posterior")
abline(v = sig2_mean, col = "red", lty = 2)
abline(v = sig2_q,    col = "red", lty = 3)
par(mfrow = c(1, 1))

#  10. CAMPO ESTIMADO EN GRILLA REGULAR 
# Proyección de la media posterior del campo x(s) a una grilla de 0.04°

lon_grid <- seq(-68.0, -64.5, by = 0.04)
lat_grid <- seq(17.5,  19.0,  by = 0.04)
grilla   <- expand.grid(lon = lon_grid, lat = lat_grid)

A_pred <- INLA::inla.spde.make.A(
  mesh = mesh,
  loc  = as.matrix(grilla[, c("lon", "lat")])
)

field_mean <- result$summary.random$field$mean
field_sd   <- result$summary.random$field$sd

grilla$x_mean <- as.numeric(A_pred %*% field_mean)
grilla$x_sd   <- as.numeric(A_pred %*% field_sd)

# Intensidad esperada de sismos por celda de 0.01 grados²
mu_hat        <- result$summary.fixed["intercepto", "mean"]
grilla$lambda <- exp(log(0.01) + mu_hat + grilla$x_mean)

# 11. MAPAS 

# 11a. Log-intensidad media posterior 
p_campo <- ggplot() +
  geom_raster(data = grilla,
              aes(lon, lat, fill = x_mean),
              interpolate = TRUE) +
  scale_fill_viridis_c(option = "plasma",
                       name   = "log-intensidad\n(media post.)") +
  geom_point(data = sismos,
             aes(lon_bin, lat_bin, size = n_sismos),
             alpha = 0.25, colour = "white", shape = 1) +
  scale_size_continuous(range = c(0.5, 4), name = "Sismos\npor celda") +
  coord_fixed() +
  xlim(-68.0, -64.5) + ylim(17.5, 19.0) +
  labs(
    title    = "Campo SPDE-Matérn: Log-Intensidad Sísmica — PRSN",
    subtitle = sprintf(
      "rho = %.2f°  |  sigma = %.2f  |  %d celdas  |  %d sismos",
      rho_mean, sqrt(sig2_mean), nrow(sismos), nrow(sismos_raw)
    ),
    x       = "Longitud",
    y       = "Latitud",
    caption = "INLA inla.spde2.pcmatern (alpha=2) | PC priors | Celdas 0.1°"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p_campo)

#  11b. Incertidumbre (SD posterior del campo) 
p_sd <- ggplot() +
  geom_raster(data = grilla,
              aes(lon, lat, fill = x_sd),
              interpolate = TRUE) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name   = "SD posterior") +
  coord_fixed() +
  xlim(-68.0, -64.5) + ylim(17.5, 19.0) +
  labs(
    title    = "Incertidumbre del Campo SPDE (SD posterior)",
    subtitle = "Mayor SD donde hay menos sismos observados",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p_sd)

#  11c. Intensidad esperada de sismos por celda 
p_lambda <- ggplot() +
  geom_raster(data = grilla,
              aes(lon, lat, fill = lambda),
              interpolate = TRUE) +
  scale_fill_viridis_c(option = "inferno",
                       name   = "lambda\n(sismos/celda)") +
  coord_fixed() +
  xlim(-68.0, -64.5) + ylim(17.5, 19.0) +
  labs(
    title    = "Intensidad Esperada de Sismos por Celda",
    subtitle = "exp(offset + mu + x(s))",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p_lambda)

# 12. EXPORTAR RESULTADOS A JSON
res_json <- list(
  modelo = list(
    tipo       = "SPDE Matérn — INLA inla.spde2.pcmatern (alpha=2, nu=1)",
    familia    = "Poisson (log link) con offset log(area)",
    n_sismos   = nrow(sismos_raw),
    n_celdas   = nrow(sismos),
    resolucion = "0.1 grados (~11 km)",
    n_nodos    = mesh$n,
    dic        = round(result$dic$dic,   2),
    waic       = round(result$waic$waic, 2),
    mlik       = round(result$mlik[1],   2)
  ),
  efectos_fijos = lapply(
    rownames(result$summary.fixed), function(nm) list(
      nombre = nm,
      mean   = round(result$summary.fixed[nm, "mean"],       4),
      sd     = round(result$summary.fixed[nm, "sd"],         4),
      q025   = round(result$summary.fixed[nm, "0.025quant"], 4),
      q975   = round(result$summary.fixed[nm, "0.975quant"], 4)
    )
  ),
  campo_spde = list(
    rho = list(
      descripcion = "Rango espacial en grados (~111 km/grado)",
      mean        = round(rho_mean,    4),
      q025        = round(rho_q[1],    4),
      q975        = round(rho_q[2],    4)
    ),
    sigma = list(
      descripcion = "Desviacion estandar marginal del campo",
      mean        = round(sqrt(sig2_mean), 4),
      q025        = round(sqrt(sig2_q[1]), 4),
      q975        = round(sqrt(sig2_q[2]), 4)
    )
  )
)

jsonlite::write_json(res_json, "resultados_spde_prsn.json",
                     pretty = TRUE, auto_unbox = TRUE)
cat("\nResultados exportados → resultados_spde_prsn.json\n")


# ============================================================
#  CATALOGO SISMICO RSPR - PANEL DE EXPLORACION Y AUDITORIA
#  Red Sismica de Puerto Rico (UPR-Mayaguez)
#
#  USO:  coloque este archivo en una carpeta junto a los .json
#        (o ponga los .json en una subcarpeta "datos/") y ejecute
#        shiny::runApp()
#HOLAMUN
# ============================================================

# ---- 0. PAQUETES ------------------------------------------
pkgs <- c("shiny","bslib","dplyr","tidyr","ggplot2","plotly","DT",
          "scales","htmltools","jsonlite","leaflet")
faltan <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltan)) install.packages(faltan)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. CARGA Y NORMALIZACION -----------------------------

# Busca los .json en "datos/" y, si no existe, en la carpeta de la app.
DATA_DIR <- if (dir.exists("datos")) "datos" else "."
ARCHIVOS <- list.files(DATA_DIR, pattern = "\\.json$", full.names = TRUE)
if (!length(ARCHIVOS))
  stop("No encuentro archivos .json en '", normalizePath(DATA_DIR), "'.")

# Union de todos los campos vistos en los catalogos RSPR.
# Los archivos 2025+ no traen ERZ/ERH/RMS/GAP/QUALITY y si traen PRSN_ID.
COLS <- c("ID","BULLETIN_NUM","DATE","TIME","LAT","LON","DEPTH","MAGNITUDE",
          "MAG_TYPE","STATIONS","PHASES","FE_REGION","STATUS","NETWORK",
          "INTENSITY","TIMEZONE","ERZ","ERH","RMS","GAP","QUALITY","SOURCE",
          "PRSN_ID")

leer_json <- function(f) {
  x <- jsonlite::fromJSON(f, simplifyDataFrame = TRUE)
  if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  # normaliza nombres: "MAG TYPE" -> "MAG_TYPE", "PRSN ID" -> "PRSN_ID"
  nm <- toupper(trimws(names(x)))
  nm <- gsub("[^A-Z0-9]+", "_", nm)
  nm <- gsub("^_+|_+$", "", nm)
  names(x) <- nm
  x[] <- lapply(x, as.character)
  # completa los campos ausentes para poder apilar esquemas distintos
  falt <- setdiff(COLS, names(x))
  for (cc in falt) x[[cc]] <- NA_character_
  x <- x[, COLS, drop = FALSE]
  x$archivo <- basename(f)
  # registro de que campos venian realmente en el archivo
  attr(x, "campos_ausentes") <- falt
  x
}

crudos <- lapply(ARCHIVOS, leer_json)

# inventario de esquema por archivo (para la pestana de auditoria)
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
    # STATUS cambia de vocabulario y de caja entre archivos
    estado    = toupper(trimws(STATUS)),
    fuente    = toupper(trimws(SOURCE)),
    # identificador estable: el ID viejo es la marca de tiempo,
    # en 2025+ eso se movio a PRSN_ID
    id_tiempo = dplyr::coalesce(PRSN_ID, ID)
  ) |>
  dplyr::filter(!is.na(dt)) |>
  dplyr::arrange(dt)

# duplicados por marca de tiempo (solapamiento entre archivos)
N_DUP <- sum(duplicated(datos$id_tiempo))

RANGO_FECHA <- range(datos$fecha, na.rm = TRUE)
RANGO_MAG   <- range(datos$mag,   na.rm = TRUE)
RANGO_PROF  <- range(datos$prof,  na.rm = TRUE)
TIPOS_MAG   <- sort(unique(na.omit(datos$tipo_mag)))
ESTADOS     <- sort(unique(na.omit(datos$estado)))

# ---- vista inicial del mapa (fija, calculada una sola vez) ----
# Se usan cuantiles y no el rango completo: el catalogo RSPR puede incluir
# eventos telesismicos, y un ajuste automatico a los extremos alejaria la
# vista hasta perder Puerto Rico de pantalla.
qq <- function(v, p) as.numeric(stats::quantile(v, p, na.rm = TRUE))
BB <- list(
  lng1 = qq(datos$lon, .005), lng2 = qq(datos$lon, .995),
  lat1 = qq(datos$lat, .005), lat2 = qq(datos$lat, .995)
)
# margen minimo por si el filtro inicial es muy estrecho
if (BB$lng2 - BB$lng1 < 0.5) { BB$lng1 <- BB$lng1 - 0.25; BB$lng2 <- BB$lng2 + 0.25 }
if (BB$lat2 - BB$lat1 < 0.5) { BB$lat1 <- BB$lat1 - 0.25; BB$lat2 <- BB$lat2 + 0.25 }

# Escalas FIJAS al catalogo completo: al filtrar cambia que se dibuja,
# no la correspondencia color-magnitud ni tamano-magnitud.
MAG_MIN  <- min(datos$mag, na.rm = TRUE)
PAL_MAG  <- leaflet::colorNumeric("YlOrRd",  domain = RANGO_MAG,  na.color = "#bbbbbb")
PAL_PROF <- leaflet::colorNumeric("viridis", domain = RANGO_PROF, na.color = "#bbbbbb",
                                  reverse = TRUE)

# ---- 2. CONSTANTES Y TEMA ---------------------------------
APP_TITULO   <- "Cat\u00e1logo S\u00edsmico RSPR"
APP_SUBTITLE <- "Red S\u00edsmica de Puerto Rico | exploraci\u00f3n y auditor\u00eda de catalogo"

UPR_ROJO   <- "#E4002B"
UPR_ROJO_2 <- "#B30021"
UPR_GRIS   <- "#6D6E71"

tema <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position  = "bottom")

gp <- function(p) plotly::ggplotly(p) |>
  plotly::config(displayModeBar = TRUE, displaylogo = FALSE)

# Variante WebGL: para dispersiones con decenas de miles de puntos.
# Permite dibujarlos todos sin muestrear.
gpw <- function(p) plotly::ggplotly(p) |>
  plotly::toWebGL() |>
  plotly::config(displayModeBar = TRUE, displaylogo = FALSE)

# validate(need(...)) NO funciona dentro de plotly::renderPlotly: la condicion
# de validacion de Shiny se propaga mal y produce "is.character(txt) is not TRUE".
# En su lugar se devuelve un grafico vacio con el mensaje como anotacion.
p_msg <- function(txt) {
  plotly::plot_ly() |>
    plotly::layout(
      xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
      annotations = list(list(text = txt, showarrow = FALSE,
                              xref = "paper", yref = "paper", x = .5, y = .5,
                              font = list(size = 13, color = "#777777")))
    ) |>
    plotly::config(displayModeBar = FALSE)
}

# ---- 3. FUNCIONES DE SISMOLOGIA ESTADISTICA ---------------

# Distribucion frecuencia-magnitud (no acumulada y acumulada)
fmd <- function(m, dm = 0.1) {
  m <- m[is.finite(m)]
  if (length(m) < 10) return(NULL)
  mr   <- round(m / dm) * dm
  grid <- seq(min(mr), max(mr), by = dm)
  n    <- as.integer(table(factor(round(mr, 3), levels = round(grid, 3))))
  data.frame(mag = grid, n = n, cum = rev(cumsum(rev(n))))
}

# Mc por maxima curvatura (Wiemer & Wyss 2000) con correccion habitual
mc_maxc <- function(f, corr = 0.2) {
  if (is.null(f)) return(NA_real_)
  f$mag[which.max(f$n)] + corr
}

# b de Aki (1965) con error de Shi & Bolt (1982)
b_aki <- function(m, mc, dm = 0.1) {
  x <- m[is.finite(m) & m >= mc - dm/2]
  n <- length(x)
  if (n < 30) return(list(b = NA_real_, se = NA_real_, a = NA_real_, n = n))
  mbar <- mean(x)
  den  <- mbar - (mc - dm/2)
  if (!is.finite(den) || den <= 0)
    return(list(b = NA_real_, se = NA_real_, a = NA_real_, n = n))
  b  <- log10(exp(1)) / den
  se <- 2.30 * b^2 * sqrt(sum((x - mbar)^2) / (n * (n - 1)))
  list(b = b, se = se, a = log10(n) + b * mc, n = n)
}

# ---- 4. ESTILO --------------------------------------------
CSS <- sprintf("
.navbar { background-color: %s !important; border: none; }
.navbar .navbar-brand, .navbar .nav-link { color: #ffffff !important; }
.navbar .nav-link.active { color: #ffffff !important; font-weight: 700;
  border-bottom: 3px solid #ffffff; }
.navbar .nav-link:hover { color: #ffe3e3 !important; }
.app-title { font-weight: 800; font-size: 1.15rem; letter-spacing: .3px; }
.app-subheader { background-color: %s; color: #ffffff; padding: 6px 18px;
  font-size: .85rem; letter-spacing: .3px; font-weight: 700; }
.app-badge { display: inline-block; border: 1.5px solid #ffffff; border-radius: 20px;
  padding: 3px 12px; color: #ffffff; font-size: .78rem; font-weight: 600; }
.card-header { font-weight: 600; }
.nota { color:#555; font-size:.8rem; margin-top:8px; }
.aviso { background:#fff3cd; border-left:4px solid #E4002B; padding:10px 14px;
  font-size:.85rem; margin-bottom:10px; }
", UPR_GRIS, UPR_ROJO_2)

JS <- "$(document).on('shown.bs.tab', function(){ window.dispatchEvent(new Event('resize')); });"

# ---- 5. UI ------------------------------------------------
ui <- bslib::page_navbar(
  title = htmltools::span(APP_TITULO, class = "app-title"),
  theme = bslib::bs_theme(version = 5, primary = UPR_ROJO),
  fillable = TRUE,
  header = htmltools::tagList(
    htmltools::tags$head(htmltools::tags$style(htmltools::HTML(CSS)),
                         htmltools::tags$script(htmltools::HTML(JS))),
    htmltools::div(class = "app-subheader", APP_SUBTITLE)
  ),
  
  # ---- filtros globales, compartidos por todas las pestanas ----
  sidebar = bslib::sidebar(
    width = 300,
    dateRangeInput("f_fecha", "Rango de fechas",
                   start = RANGO_FECHA[1], end = RANGO_FECHA[2],
                   min = RANGO_FECHA[1], max = RANGO_FECHA[2],
                   language = "es", separator = " a "),
    sliderInput("f_mag", "Magnitud",
                min = floor(RANGO_MAG[1] * 10) / 10,
                max = ceiling(RANGO_MAG[2] * 10) / 10,
                value = c(floor(RANGO_MAG[1] * 10) / 10,
                          ceiling(RANGO_MAG[2] * 10) / 10),
                step = 0.1),
    sliderInput("f_prof", "Profundidad (km)",
                min = 0, max = ceiling(max(RANGO_PROF, na.rm = TRUE)),
                value = c(0, ceiling(max(RANGO_PROF, na.rm = TRUE)))),
    selectInput("f_tipo", "Tipo de magnitud", choices = c("Todos", TIPOS_MAG)),
    selectInput("f_estado", "STATUS", choices = c("Todos", ESTADOS)),
    hr(),
    htmltools::strong("Filtros de calidad", style = "font-size:.85rem;"),
    numericInput("f_gap", "GAP azimutal m\u00e1ximo (\u00b0)", value = 360,
                 min = 0, max = 360, step = 10),
    numericInput("f_erh", "ERH m\u00e1ximo (km)", value = NA, min = 0, step = 1),
    numericInput("f_est", "M\u00ednimo de estaciones", value = 0, min = 0, step = 1),
    checkboxInput("f_conserva_na",
                  "Conservar eventos sin campos de calidad", value = TRUE),
    htmltools::div(class = "nota",
                   "Los archivos 2025+ no traen ERZ/ERH/RMS/GAP/QUALITY. Si desmarca la casilla, ",
                   "esos eventos quedan excluidos al aplicar filtros de calidad.")
  ),
  
  # ===================== AUDITORIA =====================
  bslib::nav_panel(
    "Auditor\u00eda",
    bslib::layout_columns(
      col_widths = c(3,3,3,3),
      bslib::value_box("Eventos (filtrados)", textOutput("vb_n"), theme = "primary"),
      bslib::value_box("Rango temporal",      textOutput("vb_rango")),
      bslib::value_box("Magnitud m\u00e1xima", textOutput("vb_mmax")),
      bslib::value_box("Duplicados por marca de tiempo", textOutput("vb_dup"))
    ),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(
        bslib::card_header("Esquema por archivo"),
        DT::DTOutput("t_esquema"),
        htmltools::div(class = "nota",
                       "Campos ausentes por archivo. Un esquema no uniforme obliga a decidir ",
                       "explicitamente que hacer con los eventos que carecen de metricas de calidad.")
      ),
      bslib::card(
        bslib::card_header("Cobertura temporal y huecos"),
        plotly::plotlyOutput("p_cobertura", height = "300px"),
        htmltools::div(class = "nota",
                       "Conteo diario. Las franjas sin barras son dias sin eventos registrados; ",
                       "un hueco largo suele indicar truncamiento de descarga, no quiescencia real.")
      )
    ),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(
        bslib::card_header("Huecos m\u00e1s largos sin eventos"),
        DT::DTOutput("t_huecos")
      ),
      bslib::card(
        bslib::card_header("Composici\u00f3n por a\u00f1o"),
        DT::DTOutput("t_composicion"),
        htmltools::div(class = "nota",
                       "Cambios de STATUS, SOURCE o tipo de magnitud a lo largo del tiempo ",
                       "introducen heterogeneidad correlacionada con el tiempo.")
      )
    )
  ),
  
  # ===================== MAPA =====================
  bslib::nav_panel(
    "Mapa",
    bslib::card(
      bslib::card_header(
        htmltools::div(
          style = "display:flex; align-items:center; gap:14px; flex-wrap:wrap;",
          htmltools::span("Distribuci\u00f3n epicentral"),
          radioButtons("map_color", NULL, inline = TRUE,
                       choices = c("Color por magnitud"    = "mag",
                                   "Color por profundidad" = "prof"),
                       selected = "mag"),
          actionButton("map_reset", "Vista inicial",
                       icon = shiny::icon("crosshairs"),
                       class = "btn-sm btn-outline-secondary")
        )
      ),
      shiny::uiOutput("aviso_mapa"),
      leaflet::leafletOutput("mapa", height = "620px"),
      htmltools::div(class = "nota",
                     "Color y radio siguen la magnitud en escala FIJA al cat\u00e1logo completo: ",
                     "al filtrar cambia qu\u00e9 eventos se dibujan, no la correspondencia. ",
                     "La vista no se reencuadra sola; use el bot\u00f3n para volver al encuadre inicial. ",
                     "Se dibujan todos los eventos filtrados, sin muestreo.")
    )
  ),
  
  # ===================== SERIE TEMPORAL =====================
  bslib::nav_panel(
    "Serie temporal",
    bslib::layout_columns(
      col_widths = c(12),
      bslib::card(
        bslib::card_header("Tasa de ocurrencia"),
        radioButtons("ts_res", NULL, inline = TRUE,
                     choices = c("Diaria" = "dia", "Semanal" = "semana",
                                 "Mensual" = "mes")),
        plotly::plotlyOutput("p_tasa", height = "300px")
      )
    ),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(
        bslib::card_header("Conteo acumulado"),
        plotly::plotlyOutput("p_acum", height = "330px"),
        htmltools::div(class = "nota",
                       "Los saltos verticales marcan secuencias; los tramos de pendiente ",
                       "constante marcan sismicidad de fondo aproximadamente estacionaria.")
      ),
      bslib::card(
        bslib::card_header("Magnitud contra tiempo"),
        plotly::plotlyOutput("p_mag_t", height = "330px"),
        htmltools::div(class = "nota",
                       "Un borde inferior que sube tras un evento grande es incompletitud ",
                       "dependiente de la tasa.")
      )
    )
  ),
  
  # ===================== MAGNITUD-FRECUENCIA =====================
  bslib::nav_panel(
    "Mc y valor b",
    bslib::layout_columns(
      col_widths = c(4,8),
      bslib::card(
        bslib::card_header("Par\u00e1metros"),
        checkboxInput("gr_auto", "Mc autom\u00e1tico (m\u00e1xima curvatura)", TRUE),
        sliderInput("gr_mc", "Mc manual", min = 0, max = 5, value = 2.5, step = 0.1),
        numericInput("gr_dm", "Ancho de bin (\u0394M)", value = 0.1,
                     min = 0.05, max = 0.5, step = 0.05),
        htmltools::hr(),
        htmltools::strong("Ajuste"),
        verbatimTextOutput("gr_res"),
        htmltools::div(class = "nota",
                       "Mc: maxima curvatura + 0.2. b: estimador de maxima verosimilitud ",
                       "de Aki (1965); error de Shi & Bolt (1982).")
      ),
      bslib::card(
        bslib::card_header("Distribuci\u00f3n frecuencia-magnitud"),
        plotly::plotlyOutput("p_fmd", height = "440px")
      )
    ),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(
        bslib::card_header("Estabilidad de b frente a Mc"),
        plotly::plotlyOutput("p_bstab", height = "300px"),
        htmltools::div(class = "nota",
                       "Se elige el Mc donde b se estabiliza. Si b nunca se aplana, el ",
                       "catalogo no soporta un ajuste Gutenberg-Richter en esa ventana.")
      ),
      bslib::card(
        bslib::card_header("Mc y b por a\u00f1o"),
        DT::DTOutput("t_mc_anio"),
        htmltools::div(class = "nota",
                       "Un Mc que cambia entre a\u00f1os es el problema central para cualquier ",
                       "estimacion temporal: parte de la variacion de tasa es de la red, no del proceso.")
      )
    )
  ),
  
  # ===================== CALIDAD =====================
  bslib::nav_panel(
    "Calidad de localizaci\u00f3n",
    shiny::uiOutput("aviso_calidad"),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(bslib::card_header("Error horizontal (ERH, km)"),
                  plotly::plotlyOutput("p_erh", height = "280px")),
      bslib::card(bslib::card_header("GAP azimutal (\u00b0)"),
                  plotly::plotlyOutput("p_gap", height = "280px"))
    ),
    bslib::layout_columns(
      col_widths = c(6,6),
      bslib::card(
        bslib::card_header("ERH frente a n\u00famero de estaciones"),
        plotly::plotlyOutput("p_erh_est", height = "300px"),
        htmltools::div(class = "nota",
                       "Si el error tipico de localizacion es comparable al rango espacial que ",
                       "se quiere estimar, la malla del SPDE no puede ser mas fina que ese error.")
      ),
      bslib::card(
        bslib::card_header("ERH mediano por a\u00f1o"),
        plotly::plotlyOutput("p_erh_anio", height = "300px")
      )
    )
  ),
  
  # ===================== SECUENCIA =====================
  bslib::nav_panel(
    "Secuencia / decaimiento",
    bslib::layout_columns(
      col_widths = c(4,8),
      bslib::card(
        bslib::card_header("Origen de la secuencia"),
        uiOutput("sq_sel"),
        numericInput("sq_dias", "Ventana (d\u00edas)", value = 60, min = 1, step = 1),
        numericInput("sq_mmin", "Magnitud m\u00ednima", value = 2, step = 0.1),
        sliderInput("sq_p", "Pendiente de referencia p", min = 0.5, max = 2,
                    value = 1, step = 0.05),
        htmltools::div(class = "nota",
                       "La linea punteada es una referencia de Omori-Utsu con la p elegida, ",
                       "no un ajuste. El ajuste formal es parte del modelo, no del panel.")
      ),
      bslib::card(
        bslib::card_header("Tasa de eventos tras el origen (log-log)"),
        plotly::plotlyOutput("p_omori", height = "440px")
      )
    )
  ),
  
  # ===================== TABLA =====================
  bslib::nav_panel(
    "Tabla",
    bslib::card(
      bslib::card_header("Eventos filtrados"),
      bslib::card_body(fillable = FALSE, DT::DTOutput("t_datos"))
    )
  ),
  
  bslib::nav_spacer(),
  bslib::nav_item(htmltools::span("UPR R\u00edo Piedras", class = "app-badge"))
)

# ---- 6. SERVER --------------------------------------------
server <- function(input, output, session) {
  
  # ---------- filtro global ----------
  d_f <- reactive({
    req(input$f_fecha)
    d <- datos |>
      dplyr::filter(
        fecha >= input$f_fecha[1], fecha <= input$f_fecha[2],
        !is.na(mag), mag >= input$f_mag[1], mag <= input$f_mag[2]
      )
    
    # profundidad (conserva NA)
    d <- dplyr::filter(d, is.na(prof) |
                         (prof >= input$f_prof[1] & prof <= input$f_prof[2]))
    
    if (!identical(input$f_tipo,   "Todos")) d <- dplyr::filter(d, tipo_mag == input$f_tipo)
    if (!identical(input$f_estado, "Todos")) d <- dplyr::filter(d, estado   == input$f_estado)
    
    keep_na <- isTRUE(input$f_conserva_na)
    
    if (!is.na(input$f_gap) && input$f_gap < 360)
      d <- dplyr::filter(d, (is.na(gap) & keep_na) | gap <= input$f_gap)
    
    if (!is.na(input$f_erh))
      d <- dplyr::filter(d, (is.na(erh) & keep_na) | erh <= input$f_erh)
    
    if (!is.na(input$f_est) && input$f_est > 0)
      d <- dplyr::filter(d, (is.na(est) & keep_na) | est >= input$f_est)
    
    d
  })
  
  # ---------- valores resumen ----------
  output$vb_n     <- renderText(format(nrow(d_f()), big.mark = ","))
  output$vb_rango <- renderText({
    d <- d_f(); if (!nrow(d)) return("-")
    paste(format(min(d$fecha), "%Y-%m-%d"), "a", format(max(d$fecha), "%Y-%m-%d"))
  })
  output$vb_mmax  <- renderText({
    d <- d_f(); if (!nrow(d)) return("-")
    sprintf("%.2f", max(d$mag, na.rm = TRUE))
  })
  output$vb_dup   <- renderText(format(N_DUP, big.mark = ","))
  
  # ---------- auditoria ----------
  output$t_esquema <- DT::renderDT({
    DT::datatable(ESQUEMA, rownames = FALSE,
                  colnames = c("Archivo","Eventos","Campos ausentes"),
                  options = list(dom = "t", paging = FALSE, scrollX = TRUE))
  })
  
  output$p_cobertura <- plotly::renderPlotly({
    d <- d_f(); req(nrow(d) > 0)
    x <- dplyr::count(d, fecha, name = "n")
    p <- ggplot2::ggplot(x, ggplot2::aes(fecha, n)) +
      ggplot2::geom_col(width = 1, fill = UPR_ROJO) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(x = NULL, y = "eventos/d\u00eda") + tema
    gp(p)
  })
  
  output$t_huecos <- DT::renderDT({
    d <- d_f(); req(nrow(d) > 1)
    f <- sort(unique(d$fecha))
    dif <- as.integer(diff(f))
    if (!length(dif)) return(NULL)
    o <- order(dif, decreasing = TRUE)[seq_len(min(15, length(dif)))]
    h <- data.frame(
      desde = as.character(f[o]),
      hasta = as.character(f[o + 1]),
      dias  = dif[o],
      stringsAsFactors = FALSE
    )
    h <- h[h$dias > 1, , drop = FALSE]
    DT::datatable(h, rownames = FALSE,
                  colnames = c("\u00daltimo evento","Siguiente evento","D\u00edas sin eventos"),
                  options = list(dom = "tp", pageLength = 8))
  })
  
  output$t_composicion <- DT::renderDT({
    d <- d_f(); req(nrow(d) > 0)
    x <- d |>
      dplyr::group_by(anio) |>
      dplyr::summarise(
        n        = dplyr::n(),
        estados  = paste(sort(unique(na.omit(estado))),   collapse = "/"),
        tipos    = paste(sort(unique(na.omit(tipo_mag))), collapse = "/"),
        fuentes  = paste(sort(unique(na.omit(fuente))),   collapse = "/"),
        con_erh  = sum(!is.na(erh)),
        .groups  = "drop"
      )
    DT::datatable(x, rownames = FALSE,
                  colnames = c("A\u00f1o","Eventos","STATUS","Tipo mag.","SOURCE","Con ERH"),
                  options = list(dom = "tp", pageLength = 10, scrollX = TRUE))
  })
  
  # ---------- mapa ----------
  # Se dibujan TODOS los eventos filtrados, sin muestreo.
  d_mapa <- reactive({
    d_f() |> dplyr::filter(is.finite(lat), is.finite(lon))
  })
  
  output$aviso_mapa <- renderUI({
    n   <- nrow(d_mapa())
    sin <- nrow(d_f()) - n
    htmltools::div(
      class = "nota",
      sprintf("Dibujando %s eventos%s.%s",
              format(n, big.mark = ","),
              if (sin > 0) sprintf(" (%s sin coordenadas, excluidos)",
                                   format(sin, big.mark = ",")) else "",
              if (n > 40000) " Con este volumen el primer dibujado tarda unos segundos." else "")
    )
  })
  
  output$mapa <- leaflet::renderLeaflet({
    # Se dibuja UNA sola vez y sin dependencias reactivas: por eso la vista
    # no se mueve al cambiar los filtros. Los marcadores entran por proxy.
    # preferCanvas: renderizado por canvas en vez de SVG; sin esto el navegador
    # no aguanta decenas de miles de marcadores.
    leaflet::leaflet(options = leaflet::leafletOptions(preferCanvas = TRUE)) |>
      leaflet::addProviderTiles(
        "CartoDB.Positron",
        options = leaflet::providerTileOptions(updateWhenZooming = FALSE,
                                               updateWhenIdle = TRUE)) |>
      leaflet::fitBounds(BB$lng1, BB$lat1, BB$lng2, BB$lat2)
  })
  
  observe({
    d   <- d_mapa()
    por <- input$map_color
    req(!is.null(por))
    
    pal <- if (identical(por, "prof")) PAL_PROF else PAL_MAG
    val <- if (identical(por, "prof")) d$prof   else d$mag
    ttl <- if (identical(por, "prof")) "Prof. (km)" else "Magnitud"
    
    p <- leaflet::leafletProxy("mapa", data = d) |>
      leaflet::clearMarkers() |>
      leaflet::clearControls()
    
    if (!nrow(d)) return(invisible(NULL))
    
    withProgress(message = sprintf("Dibujando %s eventos", format(nrow(d), big.mark = ",")),
                 value = 1, {
                   p |>
                     leaflet::addCircleMarkers(
                       lng = ~lon, lat = ~lat,
                       # radio anclado al minimo GLOBAL, no al del filtro
                       radius = ~pmax(2.5, (mag - MAG_MIN + 0.3) * 2.2),
                       color = pal(val), stroke = FALSE, fillOpacity = 0.55,
                       popup = ~sprintf(
                         "<b>%s</b><br/>%s %s UTC<br/>M %.2f (%s)<br/>prof %s km<br/>GAP %s | ERH %s | est %s",
                         id_tiempo, DATE, TIME, mag, tipo_mag, DEPTH,
                         ifelse(is.na(gap), "n/d", as.character(gap)),
                         ifelse(is.na(erh), "n/d", as.character(erh)),
                         ifelse(is.na(est), "n/d", as.character(est)))
                     ) |>
                     leaflet::addLegend("bottomright", pal = pal,
                                        values = if (identical(por, "prof")) RANGO_PROF else RANGO_MAG,
                                        title = ttl, opacity = 0.8)
                 })
  })
  
  # boton de reencuadre
  observeEvent(input$map_reset, {
    leaflet::leafletProxy("mapa") |>
      leaflet::flyToBounds(BB$lng1, BB$lat1, BB$lng2, BB$lat2)
  })
  
  # ---------- serie temporal ----------
  output$p_tasa <- plotly::renderPlotly({
    d <- d_f(); req(nrow(d) > 0)
    d$per <- switch(input$ts_res,
                    dia    = d$fecha,
                    semana = as.Date(cut(d$fecha, "week")),
                    mes    = as.Date(cut(d$fecha, "month")))
    x <- dplyr::count(d, per, name = "n")
    p <- ggplot2::ggplot(x, ggplot2::aes(per, n)) +
      ggplot2::geom_col(width = 1, fill = UPR_ROJO) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(x = NULL, y = "eventos") + tema
    gp(p)
  })
  
  output$p_acum <- plotly::renderPlotly({
    d <- d_f(); req(nrow(d) > 0)
    d <- dplyr::arrange(d, dt); d$acum <- seq_len(nrow(d))
    p <- ggplot2::ggplot(d, ggplot2::aes(dt, acum)) +
      ggplot2::geom_step(linewidth = .6, colour = UPR_ROJO_2) +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(x = NULL, y = "eventos acumulados") + tema
    gp(p)
  })
  
  output$p_mag_t <- plotly::renderPlotly({
    d <- d_f(); req(nrow(d) > 0)
    p <- ggplot2::ggplot(d, ggplot2::aes(dt, mag)) +
      ggplot2::geom_point(alpha = .25, size = .7, colour = UPR_GRIS) +
      ggplot2::labs(x = NULL, y = "magnitud") + tema
    gpw(p)
  })
  
  # ---------- Mc y b ----------
  fmd_r <- reactive({
    d <- d_f(); req(nrow(d) > 10)
    fmd(d$mag, dm = input$gr_dm)
  })
  
  mc_r <- reactive({
    if (isTRUE(input$gr_auto)) mc_maxc(fmd_r()) else input$gr_mc
  })
  
  observeEvent(fmd_r(), {
    f <- fmd_r(); req(!is.null(f))
    updateSliderInput(session, "gr_mc",
                      min = min(f$mag), max = max(f$mag),
                      value = if (isTRUE(input$gr_auto)) mc_maxc(f) else input$gr_mc)
  })
  
  b_r <- reactive({
    d <- d_f(); req(nrow(d) > 10)
    b_aki(d$mag, mc_r(), dm = input$gr_dm)
  })
  
  output$gr_res <- renderText({
    b <- b_r(); mc <- mc_r()
    if (is.na(b$b))
      return(sprintf("Mc = %.2f\nEventos sobre Mc: %d\nInsuficientes para estimar b.",
                     mc, b$n))
    sprintf("Mc = %.2f\nb  = %.3f  (\u00b1 %.3f)\na  = %.2f\nN(\u2265Mc) = %s",
            mc, b$b, b$se, b$a, format(b$n, big.mark = ","))
  })
  
  output$p_fmd <- plotly::renderPlotly({
    f <- fmd_r(); req(!is.null(f))
    b <- b_r(); mc <- mc_r()
    f$cum[f$cum == 0] <- NA
    p <- ggplot2::ggplot(f, ggplot2::aes(mag)) +
      ggplot2::geom_col(ggplot2::aes(y = n), fill = "#cccccc", width = input$gr_dm * .9) +
      ggplot2::geom_point(ggplot2::aes(y = cum), colour = UPR_ROJO_2, size = 1.5) +
      ggplot2::geom_vline(xintercept = mc, linetype = 2, colour = UPR_ROJO) +
      ggplot2::scale_y_log10(labels = scales::comma) +
      ggplot2::labs(x = "magnitud", y = "frecuencia (log)") + tema
    if (!is.na(b$b)) {
      gr <- data.frame(mag = seq(mc, max(f$mag), by = input$gr_dm))
      gr$cum <- 10^(b$a - b$b * gr$mag)
      p <- p + ggplot2::geom_line(data = gr, ggplot2::aes(mag, cum),
                                  colour = "#000000", linewidth = .6)
    }
    gp(p)
  })
  
  output$p_bstab <- plotly::renderPlotly({
    d <- d_f(); f <- fmd_r(); req(!is.null(f), nrow(d) > 50)
    mcs <- seq(min(f$mag), quantile(d$mag, .98, na.rm = TRUE), by = input$gr_dm)
    res <- lapply(mcs, function(m) {
      bb <- b_aki(d$mag, m, dm = input$gr_dm)
      data.frame(mc = m, b = bb$b, se = bb$se)
    })
    x <- dplyr::bind_rows(res) |> dplyr::filter(!is.na(b))
    req(nrow(x) > 0)
    p <- ggplot2::ggplot(x, ggplot2::aes(mc, b)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = b - se, ymax = b + se),
                           fill = "#dddddd") +
      ggplot2::geom_line(colour = UPR_ROJO_2, linewidth = .7) +
      ggplot2::geom_vline(xintercept = mc_r(), linetype = 2, colour = UPR_ROJO) +
      ggplot2::labs(x = "Mc de corte", y = "b") + tema
    gp(p)
  })
  
  output$t_mc_anio <- DT::renderDT({
    d <- d_f(); req(nrow(d) > 0)
    x <- lapply(split(d, d$anio), function(g) {
      f  <- fmd(g$mag, dm = input$gr_dm)
      mc <- mc_maxc(f)
      bb <- if (is.na(mc)) list(b = NA, se = NA, n = 0) else b_aki(g$mag, mc, input$gr_dm)
      data.frame(anio = g$anio[1], n = nrow(g),
                 Mc = round(mc, 2), b = round(bb$b, 3), se = round(bb$se, 3),
                 n_sobre_Mc = bb$n)
    })
    DT::datatable(dplyr::bind_rows(x), rownames = FALSE,
                  colnames = c("A\u00f1o","Eventos","Mc","b","EE(b)","N(\u2265Mc)"),
                  options = list(dom = "tp", pageLength = 10))
  })
  
  # ---------- calidad ----------
  hay_calidad <- reactive(sum(!is.na(d_f()$erh)) > 0)
  
  output$aviso_calidad <- renderUI({
    d <- d_f()
    sin <- sum(is.na(d$erh)); tot <- nrow(d)
    if (tot == 0) return(NULL)
    if (sin == tot)
      htmltools::div(class = "aviso",
                     "Ninguno de los eventos filtrados trae campos de calidad. ",
                     "Los archivos 2025+ no los incluyen.")
    else if (sin > 0)
      htmltools::div(class = "aviso",
                     sprintf("%s de %s eventos filtrados (%.1f%%) no traen ERH/GAP/RMS.",
                             format(sin, big.mark = ","), format(tot, big.mark = ","),
                             100 * sin / tot))
  })
  
  hist_q <- function(v, lab, bins = 50) {
    force(v); force(lab)
    plotly::renderPlotly({
      d <- d_f(); x <- d[[v]]; x <- x[is.finite(x)]
      if (!length(x)) return(p_msg("Sin datos de calidad en la selecci\u00f3n."))
      p <- ggplot2::ggplot(data.frame(x = x), ggplot2::aes(x)) +
        ggplot2::geom_histogram(bins = bins, fill = UPR_ROJO, colour = NA) +
        ggplot2::labs(x = lab, y = "eventos") + tema
      gp(p)
    })
  }
  output$p_erh <- hist_q("erh", "ERH (km)")
  output$p_gap <- hist_q("gap", "GAP (\u00b0)")
  
  output$p_erh_est <- plotly::renderPlotly({
    d <- d_f() |> dplyr::filter(is.finite(erh), is.finite(est))
    if (!nrow(d)) return(p_msg("Sin datos de calidad en la selecci\u00f3n."))
    p <- ggplot2::ggplot(d, ggplot2::aes(est, erh)) +
      ggplot2::geom_point(alpha = .2, size = .8, colour = UPR_GRIS) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = "estaciones", y = "ERH (km, log)") + tema
    gpw(p)
  })
  
  output$p_erh_anio <- plotly::renderPlotly({
    d <- d_f() |> dplyr::filter(is.finite(erh))
    if (!nrow(d)) return(p_msg("Sin datos de calidad en la selecci\u00f3n."))
    x <- d |> dplyr::group_by(anio) |>
      dplyr::summarise(mediana = median(erh), q1 = quantile(erh, .25),
                       q3 = quantile(erh, .75), .groups = "drop")
    p <- ggplot2::ggplot(x, ggplot2::aes(anio, mediana)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q1, ymax = q3), fill = "#dddddd") +
      ggplot2::geom_line(colour = UPR_ROJO_2, linewidth = .8) +
      ggplot2::geom_point(colour = UPR_ROJO_2) +
      ggplot2::labs(x = NULL, y = "ERH mediano (km)") + tema
    gp(p)
  })
  
  # ---------- secuencia / decaimiento ----------
  output$sq_sel <- renderUI({
    d <- d_f(); req(nrow(d) > 0)
    top <- d |> dplyr::arrange(dplyr::desc(mag)) |> head(20)
    ch <- setNames(as.character(top$dt),
                   sprintf("M%.2f  %s  (%.3f, %.3f)",
                           top$mag, format(top$dt, "%Y-%m-%d %H:%M"),
                           top$lat, top$lon))
    selectInput("sq_t0", "Evento origen (20 mayores del filtro)", choices = ch)
  })
  
  output$p_omori <- plotly::renderPlotly({
    req(input$sq_t0)
    t0 <- as.POSIXct(input$sq_t0, tz = "UTC")
    d  <- d_f() |>
      dplyr::filter(dt > t0,
                    dt <= t0 + input$sq_dias * 86400,
                    mag >= input$sq_mmin)
    if (nrow(d) <= 20)
      return(p_msg("Muy pocos eventos tras el origen para ver decaimiento."))
    
    tt <- as.numeric(difftime(d$dt, t0, units = "days"))
    tt <- tt[tt > 0]
    if (length(tt) < 5 || diff(range(tt)) <= 0)
      return(p_msg("Ventana temporal insuficiente."))
    brk <- 10^seq(log10(min(tt)), log10(max(tt)), length.out = 25)
    ct  <- cut(tt, brk, include.lowest = TRUE)
    n   <- as.integer(table(ct))
    anc <- diff(brk)
    ctr <- sqrt(brk[-1] * brk[-length(brk)])
    x   <- data.frame(t = ctr, tasa = n / anc)
    x   <- x[is.finite(x$tasa) & x$tasa > 0, ]
    if (nrow(x) <= 3) return(p_msg("Bins insuficientes."))
    
    ref <- data.frame(t = x$t, tasa = max(x$tasa) * (x$t / min(x$t))^(-input$sq_p))
    p <- ggplot2::ggplot(x, ggplot2::aes(t, tasa)) +
      ggplot2::geom_point(colour = UPR_ROJO_2, size = 1.8) +
      ggplot2::geom_line(data = ref, linetype = 2, colour = "#000000") +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::labs(x = "d\u00edas desde el origen (log)",
                    y = "eventos/d\u00eda (log)") + tema
    gp(p)
  })
  
  # ---------- tabla ----------
  output$t_datos <- DT::renderDT({
    d <- d_f() |>
      dplyr::select(id_tiempo, DATE, TIME, lat, lon, prof, mag, tipo_mag,
                    est, fases, erh, erz, rms, gap, calidad, estado,
                    fuente, archivo)
    DT::datatable(
      d, rownames = FALSE, filter = "top", extensions = "Buttons",
      colnames = c("ID tiempo","Fecha","Hora UTC","Lat","Lon","Prof","Mag",
                   "Tipo","Est.","Fases","ERH","ERZ","RMS","GAP","Cal.",
                   "STATUS","SOURCE","Archivo"),
      options = list(pageLength = 25, scrollX = TRUE, dom = "Blfrtip",
                     buttons = c("copy","csv","excel"))
    )
  })
}

shinyApp(ui, server)
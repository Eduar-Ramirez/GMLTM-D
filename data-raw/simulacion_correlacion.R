# =====================================================================
# Estudio de recuperacion de parametros para GMLTM_corr()
# =====================================================================
#
# QUE HACE ESTE SCRIPT
# --------------------
# Simula datos bajo el GMLTM-D con estructura latente CORRELACIONADA
# (theta_i ~ N_M(0, Sigma), Sigma con correlacion rho_true entre los dos
# componentes cognitivos), ajusta GMLTM_corr() a cada conjunto simulado,
# y evalua que tan bien el modelo recupera los parametros verdaderos.
#
# Es analogo a la Fase 1 del Capitulo 3 (recuperacion de parametros del
# GMLTM-D base), pero el foco aqui es especificamente Sigma[1,2]: la
# correlacion entre los componentes "global" y "local", que es la
# innovacion (iii) del modelo (Sigma general en vez de Sigma = sigma^2*I).
#
# Se usan la MISMA matriz Q (27 items x 5 reglas) y los MISMOS componentes
# (global = reglas 1-3, local = reglas 4-5) del tutorial del Capitulo 6 /
# la vinieta del paquete, para que los resultados sean directamente
# comparables con el resto de la tesis.
#
# COMO LEER LAS METRICAS DE RECUPERACION (para preparar la defensa)
# -------------------------------------------------------------------
# Para cada parametro (theta, eta, alpha, c, Sigma[1,2]) se calculan tres
# metricas estandar, comparando el valor ESTIMADO (media posterior, EAP)
# contra el valor VERDADERO usado para simular los datos:
#
#   - bias (sesgo medio) = mean(estimado - verdadero)
#       Un sesgo cercano a 0 indica que el modelo no sobre- ni sub-estima
#       sistematicamente el parametro. Un sesgo positivo/negativo grande
#       indica una tendencia sistematica del modelo a sobre/subestimar.
#
#   - RMSE (raiz del error cuadratico medio) = sqrt(mean((estimado-verdadero)^2))
#       Mide la magnitud TIPICA del error de estimacion (sesgo + varianza).
#       A diferencia del sesgo, el RMSE no se cancela cuando los errores
#       son de distinto signo, por lo que es la metrica principal de
#       "que tan cerca esta la estimacion del valor real, en promedio".
#
#   - cobertura del IC 95% = proporcion de casos en que el valor VERDADERO
#       cae dentro del intervalo de credibilidad del 95% estimado
#       (cuantiles 0.025 y 0.975 del posterior).
#       Si el modelo esta bien calibrado, la cobertura deberia ser
#       cercana a 0.95 (no significativamente menor). Cobertura muy por
#       debajo de 0.95 sugiere que los intervalos son demasiado angostos
#       (exceso de confianza); cobertura muy por encima sugiere
#       intervalos demasiado anchos (conservadores).
#
# IMPORTANTE sobre Sigma[1,2] y `n_reps`: theta, eta, alpha y c tienen
# MUCHAS instancias dentro de una sola corrida (un theta por sujeto y
# componente, un alpha por regla identificada, etc.), asi que su
# cobertura ya es una proporcion informativa incluso con una sola
# replica por condicion (n_reps = 1). Sigma[1,2], en cambio, es UN SOLO
# numero por corrida: con n_reps = 1 la "cobertura" para Sigma es solo un
# indicador binario (0 o 1: el IC 95% de esa corrida particular contuvo o
# no el rho_true usado), no una proporcion real. Para obtener una
# cobertura de Sigma[1,2] interpretable como porcentaje (p. ej. "8 de 10
# replicas contuvieron el valor verdadero") hay que subir `n_reps` (p.
# ej. 10-20), lo cual multiplica el tiempo de computo por ese factor.
# Por defecto se deja n_reps = 1 por costo computacional; subirlo si el
# tiempo disponible lo permite.
#
# COMO AJUSTAR EL COSTO COMPUTACIONAL
# ------------------------------------
# Los parametros de la seccion "CONFIGURACION" controlan directamente el
# tiempo de computo. Los valores por defecto (N = 300, cadenas completas
# 2000+1000, 2 cadenas) siguen la especificacion del estudio de
# recuperacion, pero ajustar iters/iter_warmup/chains hacia abajo (o
# n_reps hacia arriba) segun el tiempo disponible es razonable para una
# corrida exploratoria antes de la corrida final que se reporte.
#
# =====================================================================

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
})

# ============================== CONFIGURACION =========================

N_subj      <- 300                 # sujetos por condicion simulada
rho_true    <- c(0, 0.3, 0.6)       # correlaciones verdaderas a evaluar
n_reps      <- 1                    # replicas por condicion (ver nota arriba)
iters       <- 2000                 # iteraciones post-warmup por cadena
iter_warmup <- 1000                 # iteraciones de warmup por cadena
chains      <- 2                    # numero de cadenas
cores       <- min(chains, parallel::detectCores() - 1)
credible    <- 0.95                 # ancho del IC usado para cobertura
seed_base   <- 20240819             # semilla base (se le suma un offset por condicion/replica)

output_dir  <- "data-raw"
csv_path    <- file.path(output_dir, "simulacion_correlacion_resultados.csv")
plot_path   <- file.path(output_dir, "simulacion_correlacion_resultados.png")

# Misma Q (27 items x 5 reglas) y componentes del tutorial del Capitulo 6
Q <- structure(
  c(0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,
    1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
    1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,
    0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,
    1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1),
  dim = c(27L, 5L),
  dimnames = list(NULL, c("rot_fig", "rot_trap", "reflection",
                          "subt_seg", "mov_point")))
components <- list(global = c(1, 2, 3), local = c(4, 5))

# ============================== FUNCIONES AUXILIARES ==================

# Simula un conjunto de datos binarios bajo el GMLTM-D con Sigma general.
# Los parametros verdaderos se muestrean de las MISMAS priors que usa el
# modelo (eta ~ N(0,1), alpha ~ N+(0,1) es decir half-normal, c ~ Beta(3,20)),
# y theta se genera con la correlacion rho_true especificada entre
# "global" y "local".
simulate_gmltm_corr_data <- function(Q, components, N_subj, rho_true, seed) {
  set.seed(seed)

  binary_eta    <- get_eta(components)
  C             <- get_C(Q, components)
  indexes_eta   <- locate_eta(binary_eta)
  indexes_alpha <- locate_alpha(Q, components)
  M       <- length(components)
  K       <- ncol(Q)
  N_item  <- nrow(Q)
  n_eta   <- nrow(indexes_eta)
  # ATENCION: nrow(indexes_alpha) es el numero de pares (item, componente)
  # con alpha no nulo, NO el numero de parametros alpha DISTINTOS del
  # modelo (varios items pueden compartir el mismo alpha, ver locate_alpha()
  # en R/GMLTM.R). El numero real de parametros alpha que el modelo estima
  # es max(indexes_alpha[, 3]) -- ese es el que hay que muestrear aqui.
  n_alpha     <- nrow(indexes_alpha)
  max_n_alpha <- max(indexes_alpha[, 3])

  # --- parametros verdaderos, muestreados de las priors del modelo ---
  eta_true   <- stats::rnorm(n_eta, mean = 0, sd = 1)
  alpha_true <- abs(stats::rnorm(max_n_alpha, mean = 0, sd = 1))   # half-normal(0,1)
  c_true     <- stats::rbeta(N_item, shape1 = 3, shape2 = 20)

  # --- theta con correlacion rho_true entre componentes (M = 2 aqui) ---
  Sigma_true <- diag(M)
  Sigma_true[1, 2] <- Sigma_true[2, 1] <- rho_true
  L <- t(chol(Sigma_true))                 # L %*% t(L) = Sigma_true
  theta_raw <- matrix(stats::rnorm(N_subj * M), N_subj, M)
  theta_true <- theta_raw %*% t(L)         # theta_true[i, ] = L %*% theta_raw[i, ]

  # --- reconstruir eta_matrix, beta, alpha_matrix a partir de los indices ---
  eta_matrix <- matrix(0, K, M)
  for (i in seq_len(n_eta)) eta_matrix[indexes_eta[i, 1], indexes_eta[i, 2]] <- eta_true[i]
  beta_true <- Q %*% eta_matrix

  alpha_matrix <- matrix(0, N_item, M)
  for (i in seq_len(n_alpha)) alpha_matrix[indexes_alpha[i, 1], indexes_alpha[i, 2]] <- alpha_true[indexes_alpha[i, 3]]

  # --- generar respuestas binarias bajo el modelo GMLTM-D ---
  ID   <- rep(seq_len(N_subj), times = N_item)
  item <- rep(seq_len(N_item), each = N_subj)
  mu <- plogis(alpha_matrix[item, ] * (theta_true[ID, ] - beta_true[item, ])) ^ C[item, ]
  p  <- c_true[item] + (1 - c_true[item]) * apply(mu, 1, prod)
  y  <- stats::rbinom(length(p), size = 1, prob = p)

  data_sim <- as.data.frame(matrix(y, nrow = N_subj, ncol = N_item))
  colnames(data_sim) <- paste0("item", seq_len(N_item))
  rownames(data_sim) <- paste0("subj", seq_len(N_subj))

  list(
    data = data_sim,
    truth = list(
      theta = theta_true, eta = eta_true, alpha = alpha_true, c = c_true,
      Sigma_12 = rho_true
    ),
    indexes = list(indexes_eta = indexes_eta, indexes_alpha = indexes_alpha,
                   n_eta = n_eta, n_alpha = n_alpha)
  )
}

# Sesgo medio, RMSE y cobertura del IC de `credible`, dado un vector de
# valores verdaderos y los vectores EAP/limite-inferior/limite-superior
# correspondientes (mismo orden y longitud).
recovery_metrics <- function(true, est, lower, upper) {
  data.frame(
    bias     = mean(est - true),
    rmse     = sqrt(mean((est - true)^2)),
    coverage = mean(true >= lower & true <= upper),
    n        = length(true)
  )
}

# ============================== SIMULACION + AJUSTE ====================

results <- list()
row_i <- 0

for (cond_i in seq_along(rho_true)) {
  rho <- rho_true[cond_i]

  # acumuladores para agrupar (pool) las n_reps replicas de esta condicion
  pool <- list(
    theta = list(true = c(), est = c(), lower = c(), upper = c()),
    eta   = list(true = c(), est = c(), lower = c(), upper = c()),
    alpha = list(true = c(), est = c(), lower = c(), upper = c()),
    c     = list(true = c(), est = c(), lower = c(), upper = c()),
    Sigma_12 = list(true = c(), est = c(), lower = c(), upper = c())
  )

  for (rep_i in seq_len(n_reps)) {
    seed <- seed_base + 1000 * cond_i + rep_i
    message(sprintf(
      "== Condicion rho_true = %.2f | replica %d/%d | seed = %d ==",
      rho, rep_i, n_reps, seed
    ))

    sim <- simulate_gmltm_corr_data(Q, components, N_subj = N_subj,
                                     rho_true = rho, seed = seed)

    fit <- GMLTM_corr(sim$data, Q, components,
                       iters = iters, iter_warmup = iter_warmup,
                       chains = chains, cores = cores, quantiles = c(0.025, 0.5, 0.975))

    # --- theta: N_subj x M valores ---
    pool$theta$true  <- c(pool$theta$true,  as.vector(sim$truth$theta))
    pool$theta$est   <- c(pool$theta$est,   as.vector(fit$EAP$theta))
    pool$theta$lower <- c(pool$theta$lower, as.vector(fit$quantiles$theta[["0.025"]]))
    pool$theta$upper <- c(pool$theta$upper, as.vector(fit$quantiles$theta[["0.975"]]))

    # --- eta: n_eta valores (un valor por combinacion regla-componente) ---
    idx_eta <- sim$indexes$indexes_eta
    eta_est   <- fit$EAP$eta[idx_eta[, c(1, 2), drop = FALSE]]
    eta_lower <- fit$quantiles$eta[["0.025"]][idx_eta[, c(1, 2), drop = FALSE]]
    eta_upper <- fit$quantiles$eta[["0.975"]][idx_eta[, c(1, 2), drop = FALSE]]
    pool$eta$true  <- c(pool$eta$true,  sim$truth$eta)
    pool$eta$est   <- c(pool$eta$est,   eta_est)
    pool$eta$lower <- c(pool$eta$lower, eta_lower)
    pool$eta$upper <- c(pool$eta$upper, eta_upper)

    # --- alpha: n_alpha valores (un valor por grupo de items con el mismo alpha) ---
    idx_alpha <- sim$indexes$indexes_alpha
    first_row_per_alpha <- !duplicated(idx_alpha[, 3])
    idx_unique <- idx_alpha[first_row_per_alpha, c(1, 2), drop = FALSE]
    alpha_est   <- fit$EAP$alpha[idx_unique]
    alpha_lower <- fit$quantiles$alpha[["0.025"]][idx_unique]
    alpha_upper <- fit$quantiles$alpha[["0.975"]][idx_unique]
    pool$alpha$true  <- c(pool$alpha$true,  sim$truth$alpha)
    pool$alpha$est   <- c(pool$alpha$est,   alpha_est)
    pool$alpha$lower <- c(pool$alpha$lower, alpha_lower)
    pool$alpha$upper <- c(pool$alpha$upper, alpha_upper)

    # --- c: N_item valores ---
    pool$c$true  <- c(pool$c$true,  sim$truth$c)
    pool$c$est   <- c(pool$c$est,   fit$EAP$guessing)
    pool$c$lower <- c(pool$c$lower, fit$quantiles$guessing[["0.025"]])
    pool$c$upper <- c(pool$c$upper, fit$quantiles$guessing[["0.975"]])

    # --- Sigma[1,2]: un valor por replica (foco principal del estudio) ---
    corr_out <- extract_correlation(fit, credible = credible, plot = FALSE)
    pool$Sigma_12$true  <- c(pool$Sigma_12$true,  rho)
    pool$Sigma_12$est   <- c(pool$Sigma_12$est,   corr_out$pairs$estimate[1])
    pool$Sigma_12$lower <- c(pool$Sigma_12$lower, corr_out$pairs$lower[1])
    pool$Sigma_12$upper <- c(pool$Sigma_12$upper, corr_out$pairs$upper[1])
  }

  for (param_name in names(pool)) {
    row_i <- row_i + 1
    metrics <- recovery_metrics(pool[[param_name]]$true, pool[[param_name]]$est,
                                 pool[[param_name]]$lower, pool[[param_name]]$upper)
    results[[row_i]] <- cbind(
      data.frame(rho_true = rho, parameter = param_name, n_reps = n_reps),
      metrics
    )
  }
}

recovery_table <- do.call(rbind, results)
rownames(recovery_table) <- NULL

print(recovery_table)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(recovery_table, csv_path, row.names = FALSE)
message("Tabla de recuperacion guardada en: ", csv_path)

# ============================== GRAFICO =================================
# Sesgo y cobertura por condicion de rho_true, un panel por parametro.
# Estilo consistente con el resto del paquete (theme_minimal, paleta
# azul/rojo usada en R/plot_ICC.R).

param_labels <- c(
  theta = "theta", eta = "eta (dificultad de reglas)",
  alpha = "alpha (discriminacion)", c = "c (adivinanza)",
  Sigma_12 = "Sigma[1,2] (correlacion global-local)"
)
recovery_table$parameter_label <- factor(
  param_labels[recovery_table$parameter],
  levels = param_labels
)

p_bias <- ggplot2::ggplot(recovery_table, ggplot2::aes(x = factor(rho_true), y = bias)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(size = 3, color = "#1f77b4") +
  ggplot2::facet_wrap(~parameter_label, scales = "free_y") +
  ggplot2::labs(title = "Sesgo medio de recuperacion por condicion de correlacion verdadera",
                x = expression(rho[true]), y = "Sesgo (estimado - verdadero)") +
  ggplot2::theme_minimal(base_size = 13)

p_coverage <- ggplot2::ggplot(recovery_table, ggplot2::aes(x = factor(rho_true), y = coverage)) +
  ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
  ggplot2::geom_point(size = 3, color = "#d62728") +
  ggplot2::facet_wrap(~parameter_label) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(title = sprintf("Cobertura del IC %.0f%% por condicion de correlacion verdadera", 100 * credible),
                x = expression(rho[true]), y = "Cobertura observada") +
  ggplot2::theme_minimal(base_size = 13)

grDevices::png(plot_path, width = 1100, height = 900, res = 120)
gridExtra::grid.arrange(p_bias, p_coverage, nrow = 2)
grDevices::dev.off()
message("Grafico de sesgo/cobertura guardado en: ", plot_path)

# recovery_table queda disponible en el entorno para inspeccion interactiva.

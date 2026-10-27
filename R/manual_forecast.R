# Original educational functions by Sarra El Faiez and Naoufal Berquecha.
# See AUTHORS.md for contributions.
Arima_sim <- function(n, p=2, q=1, d=1, ar, ma, c=0, sigma=1,burnin=300, seed=NULL, x0=0)
{
  stopifnot(d == 1)  # Dans notre projet on fixe d=1 (stationnaire en différence)
  if (!is.null(seed)) set.seed(seed)
  
  m <- n + burnin                    # longueur totale simulée avant tronquage
  eps <- rnorm(m, mean = 0, sd = sigma)
  
  # Simule y_t = diff(X_t) via ARMA(p,q) avec constante c
  y <- numeric(m)
  e <- eps                           # innovations (on garde le nom 'e' pour clarté)
  for (t in 1:m){
    # partie AR
    ar_part <- 0
    if (p > 0){
      for (i in 1:p){
        if (t - i > 0) ar_part <- ar_part + ar[i] * y[t - i]
      }
    }
    # partie MA : e_t + sum theta_j e_{t-j}
    ma_part <- e[t]
    if (q > 0){
      for (j in 1:q){
        if (t - j > 0) ma_part <- ma_part + ma[j] * e[t - j]
      }
    }
    y[t] <- c + ar_part + ma_part
  }
  
  # Intégration pour obtenir X_t (niveau)
  x_full <- x0 + cumsum(y)
  
  # Retrait du burn-in et sortie
  x <- x_full[(burnin+1):(burnin+n)]
  y_keep <- y[(burnin+1):(burnin+n)]
  
  out <- list(
    x = ts(x),
    y = ts(y_keep),         # la série différenciée simulée
    ar = ar, ma = ma, c = c, sigma = sigma,
    order = c(p, 1, q),
    d = 1, x0 = x0, burnin = burnin
  )
  class(out) <- c("ArimaSim", class(out))
  return(out)
}
# Calcule les poids psi_k (ARMA) puis psi*_k (ARIMA d=1) jusqu'à k_max
.psi_weights <- function(ar, ma, k_max) {
  p <- length(ar); q <- length(ma)
  psi <- numeric(k_max + 1)
  psi[1] <- 1  # psi_0
  # pour k >= 1
  for (k in 1:k_max) {
    theta_k <- if (k <= q) ma[k] else 0
    acc <- theta_k
    # somme phi_i * psi_{k-i}
    if (p > 0) {
      for (i in 1:min(p, k)) acc <- acc + ar[i] * psi[k - i + 1]
    }
    psi[k + 1] <- acc
  }
  # psi* = cumul des psi (car (1-B)^-1)
  psi_star <- cumsum(psi)
  list(psi = psi, psi_star = psi_star)
}
# Renvoie les écarts-types des erreurs de prédiction pour h=1..h :
#  - se_y[h]  = sd(Y_{T+h} - Yhat_{T+h|T})
#  - se_x[h]  = sd(X_{T+h} - Xhat_{T+h|T})
# sigma : sd des innovations (si NULL on l'estime par sd des eps_hat de Forecast_Per)
arima_d1_pred_se <- function(ar, ma, sigma, h) {
  w <- .psi_weights(ar, ma, k_max = h - 1)
  psi      <- w$psi      # longueur h
  psi_star <- w$psi_star # longueur h
  var_y <- sigma^2 * cumsum(psi^2)        # somme_{j=0}^{h-1} psi_j^2
  var_x <- sigma^2 * cumsum(psi_star^2)   # somme_{j=0}^{h-1} (psi*_j)^2
  list(
    se_y = sqrt(var_y),  # h valeurs
    se_x = sqrt(var_x)   # h valeurs
  )
}

####### --------- 2) Prévision: Forecast_Per (formules ARIMA(p,1,q)) ---------######
# x  : série au niveau (X_1, ..., X_T)
# ar : vecteur (phi_1,...,phi_p)
# ma : vecteur (theta_1,...,theta_q)
# c  : constante dans l'équation ARMA sur Y_t = diff(x)
# h  : horizon de prévision
Forecast_Per <- function(x, ar, ma, c = 0, h = 2, level = 0.95, sigma = NULL)
{
  x  <- as.numeric(x)
  Tn <- length(x)
  stopifnot(Tn >= 1)
  p  <- length(ar)
  q <- length(ma)
  
  # 1) Différencier : Y_t = X_t - X_{t-1} (ARMA(p,q))
  y  <- diff(x)                # longueur Ty = Tn-1
  Ty <- length(y)
  
  # 2) Reconstruire les innovations eps_hat (t=1..Ty)
  # Modèle ARMA sur Y:  y_t = c + sum_{i=1}^p phi_i y_{t-i} + eps_t + sum_{j=1}^q theta_j eps_{t-j}
  # <=> eps_t = y_t - c - sum phi_i y_{t-i} - sum theta_j eps_{t-j}
  eps_hat <- numeric(Ty)
  for (t in 1:Ty) {
    ar_part <- 0
    if (p > 0) {
      for (i in 1:p) {
        if (t - i > 0) ar_part <- ar_part + ar[i] * y[t - i]
      }
    }
    ma_part_past <- 0
    if (q > 0) {
      for (j in 1:q) {
        if (t - j > 0) ma_part_past <- ma_part_past + ma[j] * eps_hat[t - j]
      }
    }
    eps_hat[t] <- y[t] - c - ar_part - ma_part_past
  }
  
  # 3) Prévisions récursives de Y_{T+k} pour k=1..h
  # E[Y_{T+k} | I_T] = c + sum phi_i * Y_{T+k-i} + sum theta_j * E[eps_{T+k-j} | I_T]
  # où E[eps_{T+k-j}|I_T] = 0 si k>j  (innovation future), = eps_hat_{T - (j-k)} si k<=j (innovation passée)
  y_buf <- c(y, rep(NA_real_, h))      # tampon pour utiliser observés + prévus
  e_buf <- c(eps_hat, rep(0, h))       # futures innovations = 0 (optimum MSE)
  yhat  <- numeric(h)
  
  for (k in 1:h) {
    t <- Ty + k
    
    ar_part <- 0
    if (p > 0) {
      for (i in 1:p) {
        if (t - i > 0) ar_part <- ar_part + ar[i] * y_buf[t - i]
      }
    }
    
    ma_part <- 0
    if (q > 0) {
      for (j in 1:q) {
        if (t - j > 0) ma_part <- ma_part + ma[j] * e_buf[t - j]
      }
    }
    
    yhat[k]  <- c + ar_part + ma_part         # eps_{T+k} a pour espérance 0
    y_buf[t] <- yhat[k]                        # sauvegarder la prévision pour réutilisation
    # e_buf[t] reste 0 (innovation future à espérance nulle)
  }
  
  # 4) Retour au niveau : X_{T+h} = X_T + sum_{m=1}^h Y_{T+m}
  xhat <- x[Tn] + cumsum(yhat)
  
  # 5) Estimation de sigma si besoin (écart-type des innovations)
  if (is.null(sigma)) 
  {
    # ddl approx : Ty - (p + q + (c != 0))
    k_par <- p + q + (c != 0)
    sigma <- sqrt(sum((eps_hat)^2, na.rm = TRUE) / max(1, Ty - k_par))
  }
  
  # 6) Écarts-types & intervalles de prédiction au niveau
  ses  <- arima_d1_pred_se(ar, ma, sigma, h = h)
  z    <- qnorm(0.5 + level/2)
  pi_lower <- xhat - z * ses$se_x
  pi_upper <- xhat + z * ses$se_x
  
  list(
    yhat = yhat,
    xhat = xhat,
    eps_hat = eps_hat,
    sigma = sigma,
    se_y = ses$se_y,
    se_x = ses$se_x,
    level = level,
    pi_lower = pi_lower,
    pi_upper = pi_upper
  )
}


####### --------- 3) Démo rapide (h=2) ----------------------######
demo_partie1 <- function(n=400, p=2, q=2,ar=c(0.4,0.45),ma=c(1,0.25), seed=123, c=0.2, sigma=1, h=2, level = 0.95, plot=TRUE) {
  # 1) Simulation
  sim <- Arima_sim(n=n, p=p, q=q,ar=ar,ma=ma, d=1, c=c, sigma=sigma, seed=seed)
  x   <- sim$x
  
  # 2) Prévision manuelle (oracle: on utilise les vrais paramètres de la simulation)
  frc_man <- Forecast_Per(x, ar=sim$ar, ma=sim$ma, c=sim$c, h=h ,level = level)
  
  # 3) Affichages numériques
  df <- data.frame(
    h    = seq_len(h),
    point = as.numeric(frc_man$xhat),
    lo    = as.numeric(frc_man$pi_lower),
    hi    = as.numeric(frc_man$pi_upper),
    se    = as.numeric(frc_man$se_x)
  )
  message(sprintf("Prévisions manuelles (oracle) à h=1..%d (niveau %.0f%%)", h, 100*level))
  print(df)
  message(sprintf("Estimation sigma des innovations (si non fournie à la main) : %.4f", frc_man$sigma))
  
  
  # 4) Graphiques
  if (plot) {
    g1 <- autoplot(x) + ggtitle("X_t simulée — Demo 1 (manuel + intervalles)")
    print(g1)
    
    g2 <- ggAcf(x) + ggtitle("ACF de la série simulée X_t")
    g3 <- ggPacf(x) + ggtitle("PACF de la série simulée X_t")
    print(g2); print(g3)
    
    # Visualisation locale avec bande d'IP
    Tn    <- length(x)
    kshow <- max(50, h + 10)  # fenêtre d'affichage fin de série
    idx0  <- max(1, Tn - kshow + 1)
    x_last <- ts(x[idx0:Tn], start = time(x)[idx0], frequency = frequency(x))
    
    df_future <- data.frame(
      t    = (Tn + 1):(Tn + h),
      point = as.numeric(frc_man$xhat),
      lo    = as.numeric(frc_man$pi_lower),
      hi    = as.numeric(frc_man$pi_upper)
    )
    
    p2 <- autoplot(x_last) +
      # bande d'intervalle de prédiction
      geom_ribbon(
        data = df_future,
        aes(x = t, ymin = lo, ymax = hi),
        alpha = 0.2, inherit.aes = FALSE
      ) +
      # point forecast
      geom_line(
        data = df_future,
        aes(x = t, y = point),
        linetype = "dashed", linewidth = 0.7, inherit.aes = FALSE
      ) +
      geom_point(
        data = df_future,
        aes(x = t, y = point),
        size = 2, inherit.aes = FALSE
      ) +
      ggtitle(sprintf("Prévisions manuelles (h = 1..%d) avec IP %.0f%% — Demo 1", h, 100*level)) +
      xlab("t") + ylab("X_t")
    print(p2)
  }
  
  invisible(list(sim = sim, frc_man = frc_man, table = df))
}


####### --------- 4) Évaluation MSE (rolling-origin, horizon h=2) ---------######
mse_h2_manual <- function(x, ar, ma, c, window=250) {
  Tn <- length(x)
  stopifnot(Tn > window + 2)
  
  preds_man <- c()
  truths    <- c()
  
  for (t in seq(window, Tn - 2)) {
    x_train <- x[1:t]
    
    # Prévision manuelle à h=2 avec paramètres connus (oracle)
    fr_man <- Forecast_Per(x_train, ar=ar, ma=ma, c=c, h=2)$xhat[2]
    
    preds_man <- c(preds_man, as.numeric(fr_man))
    truths    <- c(truths,    as.numeric(x[t + 2]))
  }
  
  data.frame(
    model = "Manuel (oracle, Forecast_Per)",
    MSE   = mean((preds_man - truths)^2, na.rm = TRUE),
    MAE   = mean(abs(preds_man - truths), na.rm = TRUE),
    RMSE  = sqrt(mean((preds_man - truths)^2, na.rm = TRUE))
  )
}

mse_h2_compare <- function(x, p, q, window = 250, h = 2, include_drift = TRUE) {
  stopifnot(is.numeric(x) || is.ts(x))
  Tn <- length(x)
  stopifnot(Tn > window + h, h >= 1)
  preds  <- numeric(0)
  truths <- numeric(0)
  # --- 1) Boucle rolling-origin ---
  for (t in seq(window, Tn - h)) 
  {
    x_train <- x[1:t]
    # ARIMA(p,1,q) sur X (drift estimé automatiquement)
    fit <- forecast::Arima(x_train, order = c(p, 1, q), include.drift = include_drift)
    # Prévision à horizon h
    fr  <- forecast::forecast(fit, h = h)$mean[h]
    preds  <- c(preds, as.numeric(fr))
    truths <- c(truths, as.numeric(x[t + h]))
  }
  # --- 2) Visualisation locale (facultative) ---
  # On refait une seule estimation finale sur toute la série pour tracer les prévisions h étapes
  fit_final <- forecast::Arima(x, order = c(p, 1, q), include.drift = include_drift)
  fr_final  <- forecast::forecast(fit_final, h = h, level = 95)
  kshow <- max(50, h + 10)  # fenêtre d'affichage fin de série, identique à la démo
  plot <- autoplot(fr_final, include = kshow) +ggtitle(sprintf("Prévisions R::Arima (h = 1..%d) — fenêtre fin de série", h)) +xlab("t") + ylab("X_t")
  print(plot)
  # --- 3) Calcul des métriques ---
  mse  <- mean((preds - truths)^2, na.rm = TRUE)
  mae  <- mean(abs(preds - truths), na.rm = TRUE)
  rmse <- sqrt(mse)
  # --- 4) Résumé des performances ---
  out <- data.frame(
    model = sprintf("R::Arima sur X_t (p=%d,d=1,q=%d)%s", p, q, if (include_drift) " + drift" else ""),
    MSE   = mse,
    MAE   = mae,
    RMSE  = rmse,
    n_eval = length(truths)
  )
  rownames(out) <- NULL
  return(out)
}



poly_shifted <- function(coeffs, lag = 1, sign = -1L) {
  # construit (1 + sign * sum_i coeffs[i] * B^{lag*i}) si length(coeffs)=m
  # pour AR: sign=-1 ⇒ (1 - Σ φ_i B^{i*lag})
  # pour MA: sign=+1 ⇒ (1 + Σ θ_i B^{i*lag})
  out <- c(1)
  m <- length(coeffs)
  if (m == 0) return(out)
  # on ajoute en une fois: 1 + sign * c1 B^lag + sign * c2 B^{2*lag} + ...
  # ce n'est pas une convolution multiple, c'est juste un poly à trous
  poly2 <- rep(0, 1 + lag * m)
  poly2[1] <- 1
  for (i in seq_len(m)) poly2[lag*i + 1] <- sign * coeffs[i]
  out <- poly2
  out
}

poly_mult <- function(a, b) {
  # convolution (a * b)
  na <- length(a); nb <- length(b)
  out <- rep(0, na + nb - 1)
  for (i in seq_len(na)) {
    for (j in seq_len(nb)) {
      out[i + j - 1] <- out[i + j - 1] + a[i] * b[j]
    }
  }
  out
}

# Construit les polynômes AR et MA combinés (non-saisonnier × saisonnier)
# ar  = (φ1,...,φp),   ma  = (θ1,...,θq)
# sar = (Φ1,...,ΦP),   sma = (Θ1,...,ΘQ), période s
build_combined_polys <- function(ar = numeric(0), ma = numeric(0),
                                 sar = numeric(0), sma = numeric(0), s = 4) {
  # AR: A(B) = (1 - Σ φ_i B^i) * (1 - Σ Φ_k B^{k s})
  A_non <- poly_shifted(ar,  lag = 1, sign = -1L)
  A_sea <- poly_shifted(sar, lag = s, sign = -1L)
  A <- poly_mult(A_non, A_sea)  # 1 - Σ a_ℓ B^ℓ
  # MA: M(B) = (1 + Σ θ_j B^j) * (1 + Σ Θ_k B^{k s})
  M_non <- poly_shifted(ma,  lag = 1, sign = +1L)
  M_sea <- poly_shifted(sma, lag = s, sign = +1L)
  M <- poly_mult(M_non, M_sea)  # 1 + Σ m_ℓ B^ℓ
  
  # Pour la récurrence y_t = Σ a_ℓ y_{t-ℓ} + ε_t + Σ m_ℓ ε_{t-ℓ},
  # on extrait les coefficients de lag (hors terme constant 1).
  # Attention aux longueurs: si A = (1, -a1, -a2, ...), alors ar_comb = (a1, a2, ...)
  # Ici, par construction A = 1 - Σ a_ℓ B^ℓ ⇒ ar_comb = a_ℓ (positifs selon les signes issus de l’expansion).
  ar_comb <- if (length(A) > 1) (+1) * (-A[-1]) else numeric(0)  # car A = 1 + Σ c_ℓ B^ℓ avec c_ℓ = -a_ℓ
  ma_comb <- if (length(M) > 1) (+1) * ( M[-1]) else numeric(0)  # M = 1 + Σ m_ℓ B^ℓ ⇒ ma_comb = m_ℓ
  list(ar_comb = ar_comb, ma_comb = ma_comb,
       A = A, M = M)
}

#=============================
# SARIMA forecast (wrapper)
#=============================
# x    : série au niveau (X_t)
# ar,ma: non-saisonniers (φ, θ)
# sar,sma: saisonniers (Φ, Θ); s = période; d=1, D=1 (imposés ici)
# c    : constante dans l’équation sur w_t = Δ(Δ_s X_t) (souvent 0 après double diff)
# h    : horizon ; level : niveau des IP
# sigma: sd des innovations (optionnel, sinon estimé via Forecast_Per sur w_t)
sarima_forecast_manual <- function(x, ar = numeric(0), ma = numeric(0),
                                   sar = numeric(0), sma = numeric(0),
                                   s = 4, d = 1, D = 1,
                                   c = 0, h = 8, level = 0.95, sigma = NULL)
{
  stopifnot(d == 1, D == 1)     # cette version gère d=1, D=1 comme dans ton projet
  x <- as.numeric(x)
  Tn <- length(x)
  stopifnot(Tn >= s + 2)
  
  # 1) Différence saisonnière: x_s = ∇_s X
  x_s <- diff(x, lag = s)              # longueur Tn - s
  
  # 2) Polynômes combinés AR/MA (multiplicatifs)
  polys <- build_combined_polys(ar = ar, ma = ma, sar = sar, sma = sma, s = s)
  ar_c <- polys$ar_comb
  ma_c <- polys$ma_comb
  # Important: Forecast_Per applique Δ sur l'argument 'x' (ici x_s).
  # Donc il va modéliser w_t = Δ x_s = Δ(Δ_s X) avec ARMA(ar_c, ma_c).
  
  # 3) Prévisions sur x_s via Forecast_Per 
  fr_s <- Forecast_Per(x = x_s, ar = ar_c, ma = ma_c, c = c,
                       h = h, level = level, sigma = sigma)
  
  # À ce stade:
  #  - fr_s$xhat = prévisions de S_{T+k} := ∇_s X_{T+k}
  #  - fr_s$pi_lower / pi_upper = IP approx pour S_{T+k}
  
  # 4) Re-cumul saisonnier pour revenir à X:
  #    X_{T+k} = S_{T+k} + X_{T+k-s}   (récursif)
  X_hist <- x                             # on part de l'historique
  S_hat  <- as.numeric(fr_s$xhat)         # prévisions de ∇_s X
  X_hat  <- numeric(h)
  
  for (k in 1:h) {
    idx_back <- Tn + k - s
    base <- if (idx_back <= Tn) X_hist[idx_back] else X_hat[idx_back - Tn]
    X_hat[k] <- S_hat[k] + base
  }
  
  # Intervalles de prédiction (approx) sur X :
  # On propage les bornes S vers X en supposant 'base' déterministe (approx standard).
  lo_S <- as.numeric(fr_s$pi_lower)
  hi_S <- as.numeric(fr_s$pi_upper)
  X_lo <- numeric(h); X_hi <- numeric(h)
  for (k in 1:h) {
    idx_back <- Tn + k - s
    base <- if (idx_back <= Tn) X_hist[idx_back] else X_hat[idx_back - Tn]
    X_lo[k] <- lo_S[k] + base
    X_hi[k] <- hi_S[k] + base
  }
  
  list(
    xhat      = X_hat,            # prévisions niveau X_{T+1..T+h}
    pi_lower  = X_lo,             # IP (approx) pour X
    pi_upper  = X_hi,
    Shat      = S_hat,            # prévisions saisonnières ∇_s X
    Shat_lo   = lo_S,
    Shat_hi   = hi_S,
    ar_comb   = ar_c,             # pour inspection
    ma_comb   = ma_c,
    inner_forecast = fr_s         # objet Forecast_Per sur x_s
  )
}


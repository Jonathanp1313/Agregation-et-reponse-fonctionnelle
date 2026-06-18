library(shiny)
library(ggplot2)

# ─── Modèles de Holling ────────────────────────────────────────────────────
holling2 <- function(N, a, b) (a * N)   / (1 + b * N)
holling3 <- function(N, a, b) (a * N^2) / (1 + b * N^2)
holling4 <- function(N, a, b, c) (a * N) / (1 + b * N + c * N^2)

# ─── Métriques ──────────────────────────────────────────────────────────────
rmse <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) == 0) return(Inf)
  sqrt(mean((obs[ok] - pred[ok])^2))
}

r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  if (sum(ok) < 2) return(NA)
  1 - sum((obs[ok] - pred[ok])^2) / sum((obs[ok] - mean(obs[ok]))^2)
}

aic <- function(obs, pred, k) {
  ok <- is.finite(obs) & is.finite(pred)
  n <- sum(ok)
  rss <- sum((obs[ok] - pred[ok])^2)
  if (n <= k + 1 || rss <= 0) return(Inf)
  n * log(rss / n) + 2 * k
}

# ─── Ajustement générique par minimisation du RMSE ─────────────────────────
# cost_fn(params) -> rmse ; on tire de nombreux points de départ aléatoires
# dans [lower, upper] (en log pour rester positifs), on optimise chacun avec
# L-BFGS-B, puis on polit le meilleur résultat avec Nelder-Mead (sans bornes,
# souvent plus précis pour affiner un optimum déjà localisé).
fit_model <- function(cost_fn, lower, upper, n_starts = 150) {
  best <- list(value = Inf, par = NULL)
  set.seed(1)
  
  starts <- matrix(runif(n_starts * length(lower), lower, upper),
                   ncol = length(lower), byrow = TRUE)
  starts <- rbind((lower + upper) / 2, starts)  # un point de départ central
  
  for (i in seq_len(nrow(starts))) {
    res <- tryCatch(
      optim(starts[i, ], cost_fn, method = "L-BFGS-B",
            lower = lower, upper = upper,
            control = list(maxit = 2000)),
      error = function(e) NULL
    )
    if (!is.null(res) && is.finite(res$value) && res$value < best$value) {
      best <- res
    }
  }
  
  # Polissage final : Nelder-Mead à partir du meilleur point trouvé
  if (!is.null(best$par)) {
    polish <- tryCatch(
      optim(best$par, cost_fn, method = "Nelder-Mead",
            control = list(maxit = 5000, reltol = 1e-12)),
      error = function(e) NULL
    )
    if (!is.null(polish) && is.finite(polish$value) && polish$value < best$value) {
      # On garde le résultat poli seulement s'il respecte les bornes
      if (all(polish$par >= lower) && all(polish$par <= upper)) {
        best <- polish
      }
    }
  }
  
  best
}

couleurs_modele <- c(H2 = "#E07B39", H3 = "#3A86B4", H4 = "#4CAF7D")
noms_modele      <- c(H2 = "Holling II", H3 = "Holling III", H4 = "Holling IV")

# ─── UI ─────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  titlePanel("Réponse fonctionnelle — Ajustement Holling II / III / IV"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("file", "Fichier de simulation (.txt, colonnes : N capture)",
                accept = ".txt"),
      
      checkboxGroupInput("modeles", "Modèles à afficher",
                         choices  = noms_modele,
                         selected = noms_modele),
      
      actionButton("fit", "🎯 Ajuster (minimiser RMSE)",
                   class = "btn-primary", width = "100%"),
      
      br(), br(),
      uiOutput("statut")
    ),
    
    mainPanel(
      width = 9,
      uiOutput("meilleur_fit_bandeau"),
      plotOutput("plot", height = "500px"),
      hr(),
      tableOutput("table_metrics")
    )
  )
)

# ─── SERVER ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  data_sim <- reactive({
    req(input$file)
    df <- read.table(input$file$datapath, comment.char = "#",
                     col.names = c("N", "capture"))
    df[is.finite(df$N) & is.finite(df$capture), ]
  })
  
  # Stocke les meilleurs paramètres trouvés pour chaque modèle
  params <- reactiveValues(
    H2 = c(a = 0.1, b = 0.01),
    H3 = c(a = 0.1, b = 0.01),
    H4 = c(a = 0.1, b = 0.01, c = 0.01)
  )
  
  metrics <- reactiveVal(NULL)
  meilleur_modele <- reactiveVal(NULL)
  
  observeEvent(input$file, {
    metrics(NULL)
    meilleur_modele(NULL)
  })
  
  observeEvent(input$fit, {
    df <- tryCatch(data_sim(), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      showNotification("Charge d'abord un fichier de données.", type = "warning")
      return()
    }
    
    N <- df$N
    y <- df$capture
    
    withProgress(message = "Ajustement en cours…", value = 0, {
      
      # Holling II : a,b > 0 -> on optimise log(a), log(b)
      incProgress(0.2, detail = "Holling II")
      fit2 <- fit_model(
        cost_fn = function(lp) rmse(y, holling2(N, exp(lp[1]), exp(lp[2]))),
        lower = c(-15, -15), upper = c(8, 8)
      )
      params$H2 <- c(a = exp(fit2$par[1]), b = exp(fit2$par[2]))
      
      # Holling III
      incProgress(0.3, detail = "Holling III")
      fit3 <- fit_model(
        cost_fn = function(lp) rmse(y, holling3(N, exp(lp[1]), exp(lp[2]))),
        lower = c(-20, -15), upper = c(8, 8)
      )
      params$H3 <- c(a = exp(fit3$par[1]), b = exp(fit3$par[2]))
      
      # Holling IV : a,c > 0 (forme en cloche), b libre
      incProgress(0.3, detail = "Holling IV")
      fit4 <- fit_model(
        cost_fn = function(lp) {
          a <- exp(lp[1]); b <- lp[2]; c <- exp(lp[3])
          pred <- holling4(N, a, b, c)
          if (any(!is.finite(pred) | pred < 0)) return(1e6)
          rmse(y, pred)
        },
        lower = c(-15, -10, -15), upper = c(8, 10, 8)
      )
      params$H4 <- c(a = exp(fit4$par[1]), b = fit4$par[2], c = exp(fit4$par[3]))
      
      incProgress(0.2, detail = "Terminé")
    })
    
    # Calcul des métriques finales pour les 3 modèles
    preds <- list(
      H2 = holling2(N, params$H2["a"], params$H2["b"]),
      H3 = holling3(N, params$H3["a"], params$H3["b"]),
      H4 = holling4(N, params$H4["a"], params$H4["b"], params$H4["c"])
    )
    k <- c(H2 = 2, H3 = 2, H4 = 3)
    
    met <- do.call(rbind, lapply(names(preds), function(nm) {
      data.frame(
        Modèle = noms_modele[nm],
        a      = formatC(params[[nm]]["a"], format = "e", digits = 3),
        b      = formatC(params[[nm]]["b"], format = "e", digits = 3),
        c      = if (nm == "H4") formatC(params[[nm]]["c"], format = "e", digits = 3) else "—",
        R2     = formatC(r2(y, preds[[nm]]), format = "f", digits = 12),
        RMSE   = formatC(rmse(y, preds[[nm]]), format = "f", digits = 12),
        AIC    = round(aic(y, preds[[nm]], k[nm]), 2)
      )
    }))
    metrics(met)
    
    # Meilleur modèle = RMSE minimal (sans pénalité de complexité)
    rmse_vals <- sapply(names(preds), function(nm) rmse(y, preds[[nm]]))
    code_best <- names(which.min(rmse_vals))
    meilleur_modele(list(
      code = code_best,
      nom  = noms_modele[code_best],
      r2   = r2(y, preds[[code_best]]),
      rmse = rmse_vals[code_best]
    ))
    
    showNotification("Ajustement terminé.", type = "message")
  })
  
  output$statut <- renderUI({
    df <- tryCatch(data_sim(), error = function(e) NULL)
    if (is.null(df)) {
      tags$p(style = "color:#999; font-style:italic;", "Aucun fichier chargé.")
    } else {
      tags$p(style = "color:#4CAF7D; font-weight:600;",
             paste("✅", nrow(df), "points chargés"))
    }
  })
  
  output$meilleur_fit_bandeau <- renderUI({
    mb <- meilleur_modele()
    if (is.null(mb)) return(NULL)
    
    col <- couleurs_modele[mb$code]
    
    tags$div(
      style = sprintf(
        "background: %s12; border: 2px solid %s; border-radius: 10px;
         padding: 16px 24px; margin-bottom: 16px; text-align: center;",
        col, col
      ),
      tags$div(
        style = "font-size:13px; color:#888; text-transform:uppercase; letter-spacing:0.06em; font-weight:600;",
        "Meilleur ajustement (RMSE minimal)"
      ),
      tags$div(
        style = sprintf("font-size:28px; font-weight:800; color:%s; margin-top:4px;", col),
        mb$nom
      ),
      tags$div(
        style = "font-size:20px; font-weight:700; color:#333; margin-top:4px;",
        sprintf("R² = %.4f   ·   RMSE = %.5f", mb$r2, mb$rmse)
      )
    )
  })
  
  output$plot <- renderPlot({
    df <- tryCatch(data_sim(), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    N_seq <- seq(min(df$N), max(df$N), length.out = 300)
    
    sel <- input$modeles
    if (is.null(sel)) sel <- character(0)
    sel_codes <- names(noms_modele)[noms_modele %in% sel]
    
    courbes <- do.call(rbind, lapply(sel_codes, function(nm) {
      fitted <- switch(nm,
                       H2 = holling2(N_seq, params$H2["a"], params$H2["b"]),
                       H3 = holling3(N_seq, params$H3["a"], params$H3["b"]),
                       H4 = holling4(N_seq, params$H4["a"], params$H4["b"], params$H4["c"])
      )
      data.frame(N = N_seq, fitted = fitted, modele = nm)
    }))
    
    p <- ggplot() +
      geom_point(data = df, aes(N, capture), color = "#333", alpha = 0.7, size = 2)
    
    if (!is.null(courbes) && nrow(courbes) > 0) {
      p <- p + geom_line(data = courbes, aes(N, fitted, color = modele),
                         linewidth = 1.3, na.rm = TRUE)
    }
    
    p +
      scale_color_manual(values = couleurs_modele, labels = noms_modele) +
      theme_minimal(base_size = 14) +
      labs(x = "N (nombre de proies)", y = "Taux de capture", color = "Modèle")
  })
  
  output$table_metrics <- renderTable({
    metrics()
  }, striped = TRUE, bordered = TRUE)
}

shinyApp(ui, server)
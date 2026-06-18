
library(plotly)
library(magrittr)
library(circular)
library(dplyr)



########################### Dimension du Tore ###########################

torex <- 600  
torey <- 600
torez <- 600



################ On décide ici si les bancs se forment ##################

banc_proie <- TRUE
banc_pred  <- FALSE



######################### Paramètres de base ############################

delta_t  <- 0.1
N <- 100  # Nb proies
P <- 20  # Nb de prédateurs
T <- 60 # temps total de simulation
Tstart <- 50  # Début des interactions proies-prédateurs




############################### Vitesses ################################

vinertie <- 20
vproie <- 15
vpred <-  40




######################### Paramètres Von mises ##########################

kappa_proie_ref = 200.0  # Valeur que l'on aurait pour delta_t = 1
kappa_pred_ref  = 200.0

kappa_proie = kappa_proie_ref / delta_t
kappa_pred  = kappa_pred_ref  / delta_t




##################### Paramètres de banc des proies #####################

Rproie_att <- 100 ; Cproie_att <- 20
Rproie_al  <- 30 ; Cproie_al  <- 5
Rproie_rep <- 5 ; Cproie_rep <- 10




################### Paramètres de banc des prédateurs ####################

Rpred_att <- 200 ; Cpred_att <- 20
Rpred_al  <- 80 ; Cpred_al  <- 5
Rpred_rep <- 30 ; Cpred_rep <- 10




####### Décroissance exponentielle de l'attraction pour les bancs ########

alpha_att <- 0.05  # Décroissance exponentielle de l'attraction

epsilon = 0.01 #




##################### Interactions proies-prédateurs #####################

Rfuite     <- 200  ;   Cfuite <- 300  ;   beta_fuite <- 0.025
Rchasse <- 200     ;   Cchasse <- 250  ;   gamma_chasse  <- 0.055
Rcapt      <- 5   #distance de capture




###################### Temps de manip et estomac #########################

estomacmax  <- 2  
energie_proie <- 1 
decroissance_estomac = 0.95
Tmanip      <- 1   




############################# Fonctions annexes ################################


###################### Conditions aux bords ##############################


# Renvoie un individu sur le tore si il en sort

veriftore <- function(xt, yt, zt) {
  if (xt >  torex) xt <- xt -2 * torex
  if (xt < -torex) xt <- xt +2 * torex
  if (yt >  torey) yt <- yt -2 * torey
  if (yt < -torey) yt <- yt +2 * torey
  if (zt > torez)  zt <- zt -2 * torez
  if (zt < -torez) zt <- zt +2 * torez
  return(c(xt, yt,zt))
}

veriftore_vec <- function(position) veriftore(position[1], position[2], position[3])



# Calcule delta_x et delta_y en prenant en compte le fait qu'on est dans un tore

dxdy_tore <- function(x1, y1, z1, x2, y2, z2) {
  dx <- x2-x1
  dy <- y2-y1
  dz <- z2-z1
  
  #prise en compte du tore
  if (dx >  torex) dx <- dx - 2 * torex
  if (dx < -torex) dx <- dx + 2 * torex
  if (dy >  torey) dy <- dy - 2 * torey
  if (dy < -torey) dy <- dy + 2 * torey
  if (dz >  torez) dz <- dz - 2 * torez
  if (dz < -torez) dz <- dz + 2 * torez
  return(c(dx, dy, dz))
}

#################  Construction loi de Von Mises-Fisher #################

vmf_sample_W <- function(kappa) {
  b  <- (-2*kappa + sqrt(4*kappa^2 + 4)) / 2
  x0 <- (1 - b) / (1 + b)
  c  <- kappa * x0 + 2 * log(1 - x0^2)
  repeat {
    z <- cos(pi * runif(1))
    u <- runif(1)
    W <- (1 - (1+b)*z) / (1 - (1-b)*z)
    if (kappa * W + 2 * log(1 - x0*W) - c >= log(u)) return(W)
  }
}

von_mises_fisher <- function(theta_mu, phi_mu, kappa) {
  
  W <- vmf_sample_W(kappa)
  
  psi <- runif(1, 0, 2*pi)
  V1  <- cos(psi)
  V2  <- sin(psi)
  

  mx <- sin(theta_mu) * cos(phi_mu)
  my <- sin(theta_mu) * sin(phi_mu)
  mz <- cos(theta_mu)
  

  if (abs(mx) <= abs(my) && abs(mx) <= abs(mz)) {
    e1 <- c(0, -mz, my)
  } else if (abs(my) <= abs(mz)) {
    e1 <- c(mz, 0, -mx)
  } else {
    e1 <- c(-my, mx, 0)
  }
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- c(my*e1[3] - mz*e1[2],
          mz*e1[1] - mx*e1[3],
          mx*e1[2] - my*e1[1])
  

  s  <- sqrt(max(0, 1 - W^2))
  r  <- s * (V1*e1 + V2*e2) + W * c(mx, my, mz)
  
  theta_new <- acos(max(-1, min(1, r[3])))
  rho       <- sqrt(r[1]^2 + r[2]^2)
  phi_new   <- if (rho > 1e-10) atan2(r[2], r[1]) else phi_mu
  
  return(c(theta_new, phi_new))
}

###################### Formation des bancs ##############################


interaction_vec <- function(xmtn, ymtn, zmtn, theta_mtn, phi_mtn, indiv, Ratt, Ral, Rrep,Catt, Cal, Crep) {
  
  taille_pop <- length(xmtn)
  
  ############### Initialisation ##################
  
  
  vx_att <- 0
  vy_att <- 0
  vz_att <- 0
  
  vx_al <- 0
  vy_al <- 0
  vz_al <- 0
  
  vx_rep <- 0
  vy_rep <- 0
  vz_rep <- 0
  
  inter <- FALSE
  
  for (conj in seq_len(taille_pop)[-indiv]) {
    if (is.na(xmtn[conj]) || is.na(ymtn[conj]) || is.na(zmtn[conj])) next
    
    dist <- dxdy_tore(xmtn[indiv], ymtn[indiv], zmtn[indiv], xmtn[conj], ymtn[conj], zmtn[conj])
    
    dx <- dist[1]
    dy <- dist[2]
    dz <- dist[3]
    
    distcongenere2 <- dx*dx + dy*dy + dz*dz
    
    ################ Attraction ###################
    
    if (distcongenere2 < Ratt^2 && distcongenere2 > Ral^2) {
      poids <- Catt * exp(-alpha_att*sqrt(distcongenere2))
      
      vx_att <- vx_att + poids * dx
      vy_att <- vy_att + poids * dy
      vz_att <- vz_att + poids * dz
      
      inter <- TRUE
    }
    
    ################ Alignement ###################
    
    if (!is.na(theta_mtn[conj]) &&
        distcongenere2 < Ral^2 && distcongenere2 > Rrep^2) {
      poids <- Cal
      
      vx_al <- vx_al + poids * sin(theta_mtn[conj]) * cos(phi_mtn[conj])
      vy_al <- vy_al + poids * sin(theta_mtn[conj])* sin(phi_mtn[conj])
      vz_al <- vz_al + poids * cos(theta_mtn[conj])
      
      
      inter <- TRUE
    }
    
    ################# Répulsion ###################
    
    if (distcongenere2 < Rrep^2) {
      poids <- Crep /(sqrt(distcongenere2) + epsilon)
      
      vx_rep <- vx_rep - poids * dx
      vy_rep <- vy_rep - poids * dy
      vz_rep <- vz_rep - poids * dz
      
      inter <- TRUE
    }
  } 
  return(list(inter = inter, vx_att = vx_att, vy_att = vy_att, vz_att = vz_att, vx_al = vx_al, vy_al = vy_al, vz_al = vz_al, vx_rep = vx_rep, vy_rep = vy_rep, vz_rep = vz_rep))
}


#################### Déplacement des proies #############################

deplacement_proie <- function(i, predateurs_actifs, xproie_mtn, yproie_mtn, zproie_mtn,
                              theta_mtn, phi_mtn,
                              xpred_mtn, ypred_mtn, zpred_mtn) {
  ################ 1) Inertie ###################
  
  teta_prev <- theta_mtn[i]
  phi_prev  <- phi_mtn[i]
  vx_inertie <- vinertie * sin(teta_prev) * cos (phi_prev)
  vy_inertie <- vinertie * sin(teta_prev) * sin (phi_prev)
  vz_inertie = vinertie * cos (teta_prev)
  
  ################ 2) Fuite ###################
  
  vx_fuite <- 0; vy_fuite <- 0; vz_fuite <- 0
  presence_pred <- FALSE
  
  if (predateurs_actifs) {
    for (k in seq_len(P)) {
      dang    <- dxdy_tore(xproie_mtn[i], yproie_mtn[i], zproie_mtn[i], xpred_mtn[k], ypred_mtn[k], zpred_mtn[k])
      danger2 <- dang[1]^2 + dang[2]^2 + dang[3]^2
      
      if (danger2 < Rfuite^2) {
        poids     <- Cfuite * exp(-beta_fuite * sqrt(danger2))
        
        vx_fuite  <- vx_fuite - poids * dang[1]
        vy_fuite  <- vy_fuite - poids * dang[2]
        vz_fuite <- vz_fuite - poids * dang[3]
        
        presence_pred <- TRUE
      }
    }
  }
  
  ############ 3) Interactions de Banc ###################
  
  vx_banc <- 0; vy_banc <- 0; vz_banc <- 0
  
  if (banc_proie) {
    res <- interaction_vec(xproie_mtn, yproie_mtn, zproie_mtn, theta_mtn, phi_mtn, indiv = i,
                           Ratt = Rproie_att, Ral = Rproie_al, Rrep = Rproie_rep, Catt = Cproie_att,
                           Cal = Cproie_al, Crep = Cproie_rep)
    
    vx_banc <- res$vx_att + res$vx_al + res$vx_rep 
    vy_banc <- res$vy_att + res$vy_al + res$vy_rep 
    vz_banc <- res$vz_att + res$vz_al + res$vz_rep
    
  } 
  
  ################ 4) Somme vectorielle ###################
  ################## Extraction angle #####################
  ######### Application d'une vitesse constante ###########
  
  vx_tot <-  vx_inertie + vx_banc + vx_fuite 
  vy_tot <-  vy_inertie + vy_banc + vy_fuite
  vz_tot <- vz_inertie + vz_banc + vz_fuite
  
  norme = sqrt(vx_tot^2 + vy_tot^2 + vz_tot^2)
  
  teta_tot <- acos(max(-1, min(1,vz_tot/norme)))
  
  phi_tot <- atan2 (vy_tot, vx_tot)
  
  
  bruit    <- von_mises_fisher(teta_tot, phi_tot, kappa_proie)
  teta_tot <- bruit[1]
  phi_tot  <- bruit[2]
  
  x_prochain <- xproie_mtn[i] + vproie * sin(teta_tot) * cos(phi_tot) * delta_t
  y_prochain <- yproie_mtn[i] + vproie * sin(teta_tot) * sin(phi_tot) * delta_t
  z_prochain <- zproie_mtn[i] + vproie * cos(teta_tot) * delta_t
  
  pos <- veriftore(x_prochain, y_prochain, z_prochain)
  return(c(pos[1], pos[2], pos[3], teta_tot, phi_tot))
  
}






#################### Déplacement des prédateurs #########################

deplacement_pred <- function(k, predateurs_actifs, xproie_mtn, yproie_mtn, zproie_mtn,
                             xpred_mtn, ypred_mtn, zpred_mtn,
                             theta_mtn, phi_mtn,
                             estomac, tempsmange) {
  
  ################ 1) Inertie ###################
  
  teta_prev <- theta_mtn[k]
  phi_prev <- phi_mtn[k]
  
  vx_inertie <- vinertie * sin (teta_prev) * cos (phi_prev)
  vy_inertie <- vinertie * sin(teta_prev) * sin(phi_prev)
  vz_inertie <- vinertie * cos(teta_prev)
  
  ################# 2) Chasse ###################
  
  vx_chasse <- 0; vy_chasse <- 0; vz_chasse <-0
  
  if (predateurs_actifs && estomac[k] < estomacmax && tempsmange[k] >= Tmanip) {
    for (i in seq_len(N)) {
      if (is.na(xproie_mtn[i]) || is.na(yproie_mtn[i])) next
      
      dist  <- dxdy_tore(xpred_mtn[k], ypred_mtn[k], zpred_mtn[k], xproie_mtn[i], yproie_mtn[i], zproie_mtn[i])
      dx    <- dist[1]; dy <- dist[2]; dz <- dist[3]
      dist2 <- dx*dx + dy*dy + dz*dz
      
      if (dist2 < Rchasse^2) {
        poids        <- Cchasse * exp(-gamma_chasse * sqrt(dist2))
        
        vx_chasse <- vx_chasse + poids * dx
        vy_chasse <- vy_chasse + poids * dy
        vz_chasse <- vz_chasse + poids * dz
        
      }
    }
  }
  
  ########## 3) Interactions de Banc ############
  
  vx_banc <- 0; vy_banc <- 0; vz_banc <- 0
  
  if (banc_pred) {
    res <- interaction_vec(xpred_mtn, ypred_mtn, zpred_mtn, theta_mtn, phi_mtn, indiv = k,
                           Ratt = Rpred_att, Ral = Rpred_al, Rrep = Rpred_rep, Catt = Cpred_att,
                           Cal = Cpred_al, Crep = Cpred_rep)
    
    vx_banc <- res$vx_att + res$vx_al + res$vx_rep 
    vy_banc <- res$vy_att + res$vy_al + res$vy_rep 
    vz_banc <- res$vz_att + res$vz_al + res$vz_rep
  }
  
  ################ 4) Somme vectorielle ###################
  ################## Extraction angle #####################
  ######### Application d'une vitesse constante ###########
  
  vx_tot <-   vx_inertie + vx_banc +  vx_chasse 
  vy_tot <-   vy_inertie + vy_banc +  vy_chasse 
  vz_tot <-   vz_inertie + vz_banc +  vz_chasse
  
  norme = sqrt(vx_tot^2 + vy_tot^2 + vz_tot^2)
  teta_tot = acos(max(-1, min(1, vz_tot/norme)))
  
  phi_tot = atan2(vy_tot,vx_tot)
  
  
  bruit    <- von_mises_fisher(teta_tot, phi_tot, kappa_pred)
  teta_tot <- bruit[1]
  phi_tot  <- bruit[2]
  
  
  x_prochain <- xpred_mtn[k] + vpred * sin(teta_tot) * cos(phi_tot) * delta_t
  y_prochain <- ypred_mtn[k] + vpred * sin(teta_tot) * sin(phi_tot) * delta_t
  z_prochain <- zpred_mtn[k] + vpred * cos(teta_tot) * delta_t
  
  pos <- veriftore(x_prochain, y_prochain, z_prochain)
  return(c(pos[1], pos[2], pos[3], teta_tot, phi_tot))
}




#################### Fonction principale de simulation #########################

simul <- function(N, P, T) {
  Nb_sous_pas <- round(1/delta_t)  # nombre de sous-pas par pas de temps entier
  
  capture <- 0
  capture_totale <- c(capture)
  
  Estomac_1er <- c(0) # Juste pour vérification dans le plot
  estomac <- rep(0, P)
  tempsmange <- rep(Tmanip, P)
  
  ################ Initialisation des positions et angle ###################
  
  xproie <- matrix(0, nrow = T, ncol = N)
  yproie <- matrix(0, nrow = T, ncol = N)
  zproie <- matrix(0, nrow = T, ncol = N)
  
  xproie[1, ] <- runif(N, min = -torex, max = torex)
  yproie[1, ] <- runif(N, min = -torey, max = torey)
  zproie[1, ] <- runif(N, min = -torez, max = torez)
  
  xpred <- matrix(0, nrow = T, ncol = P)
  ypred <- matrix(0, nrow = T, ncol = P)
  zpred <- matrix(0, nrow = T, ncol = P)
  
  xpred[1, ] <- runif(P, min = -torex, max = torex)
  ypred[1, ] <- runif(P, min = -torey, max = torey)
  zpred[1, ] <- runif(P, min = -torez, max = torez)
  
  theta_proie <- runif(N, 0, pi) 
  theta_pred  <- runif(P, 0, pi)
  
  phi_proie <- runif(N, 0, 2*pi)
  phi_pred  <- runif(P, 0, 2*pi)
  
  for (t in 2:T) {
    predateurs_actifs <- (t >= Tstart)
    
    # Matrices temporaires pour les sous-pas
    
    xproie_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = N)
    yproie_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = N)
    zproie_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = N)
    
    xproie_temporaire[1,] <- xproie[t-1,]
    yproie_temporaire[1,] <- yproie[t-1,]
    zproie_temporaire[1,] <- zproie[t-1,]
    
    xpred_temporaire <- matrix(0, nrow = Nb_sous_pas + 1, ncol = P)
    ypred_temporaire <- matrix(0, nrow = Nb_sous_pas + 1, ncol = P)
    zpred_temporaire <- matrix(0, nrow = Nb_sous_pas + 1, ncol = P)
    
    xpred_temporaire[1, ] <- xpred[t-1, ]
    ypred_temporaire[1, ] <- ypred[t-1, ]
    zpred_temporaire[1, ] <- zpred[t-1, ]
    
    theta_proie_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = N)
    phi_proie_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = N)
    
    theta_pred_temporaire  <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = P)
    phi_pred_temporaire <- matrix(NA, nrow = Nb_sous_pas + 1, ncol = P)
    
    theta_proie_temporaire[1, ] <- theta_proie
    phi_proie_temporaire[1,] <- phi_proie
    
    theta_pred_temporaire[1, ]  <- theta_pred
    phi_pred_temporaire[1, ] <- phi_pred
    
    
    for (s in seq_len(Nb_sous_pas)) {
      
      ##################### Déplacement des proies ########################
      
      for (i in seq_len(N)) {
        
        if (!is.na(xproie_temporaire[s, i])) {
          
          position <- deplacement_proie(i, predateurs_actifs,
                                        xproie_temporaire[s,], yproie_temporaire[s,], zproie_temporaire[s,],
                                        theta_proie_temporaire[s,], phi_proie_temporaire[s,],
                                        xpred_temporaire[s,], ypred_temporaire[s,], zpred_temporaire[s,])
          
          xproie_temporaire[s+1, i]   <- position[1]
          yproie_temporaire[s+1, i]   <- position[2]
          zproie_temporaire[s+1, i]   <- position[3]
          
          theta_proie_temporaire[s+1, i] <- position[4]
          phi_proie_temporaire[s+1, i] <- position[5]
        }
      }
      
      ################### Déplacement des prédateurs ######################
      
      for (k in seq_len(P)) {
        position <- deplacement_pred(k, predateurs_actifs, xproie_temporaire[s,], yproie_temporaire[s,], zproie_temporaire[s,],
                                     xpred_temporaire[s, ], ypred_temporaire[s, ], zpred_temporaire[s,],
                                     theta_pred_temporaire[s,], phi_pred_temporaire[s,],
                                     estomac, tempsmange)
        
        xpred_temporaire[s+1, k]  <- position[1]
        ypred_temporaire[s+1, k]  <- position[2]
        zpred_temporaire[s+1, k] <- position[3]
        
        theta_pred_temporaire[s+1, k] <- position[4]
        phi_pred_temporaire[s+1, k] <- position[5]
      }
      
      if (predateurs_actifs) {
        ############### Recensement des proies capturées ##################
        
        for (k in seq_len(P)) {
          if (tempsmange[k] >= Tmanip && estomac[k] < estomacmax){
            for (i in seq_len(N)) {
              
              if (!is.na(xproie[t-1,i]) && !is.na(xproie_temporaire[s,i])){
                dist <- dxdy_tore(xpred_temporaire[s+1, k], ypred_temporaire[s+1, k], zpred_temporaire[s+1,k],
                                  xproie_temporaire[s+1, i], yproie_temporaire[s+1,i], zproie_temporaire[s+1,i])
                dx <- dist[1]
                dy <- dist[2]
                dz <- dist[3]
                distcapture2 <- dx*dx + dy*dy + dz*dz
                
                if (distcapture2 < Rcapt^2) {
                  
                  for (zsp in (s+1):(Nb_sous_pas + 1)){ 
                    xproie_temporaire[zsp, i] <- NA
                    yproie_temporaire[zsp, i] <- NA
                    zproie_temporaire[zsp, i] <- NA
                  }
                  
                  xproie_temporaire[s, i] <- NA
                  yproie_temporaire[s, i] <- NA
                  zproie_temporaire[s, i] <- NA
                  
                  xproie[t-1, i] <- NA
                  yproie[t-1, i] <- NA
                  zproie[t-1, i] <- NA
                  
                  capture <- capture + 1
                  estomac[k] <- estomac[k] + energie_proie
                  tempsmange[k] <- 0
                  break
                }
              }
            }
          }
        }
        tempsmange <- tempsmange + delta_t
        estomac    <- decroissance_estomac^(delta_t) * estomac
      }
      
    }
    
    ################ Positions finales du pas de temps t ################
    ############# = dernière ligne des matrices temporaires #############
    
    xproie[t, ] <- xproie_temporaire[Nb_sous_pas+1,]
    yproie[t, ] <- yproie_temporaire[Nb_sous_pas+1,]
    zproie[t, ] <- zproie_temporaire[Nb_sous_pas+1,]
    
    xpred[t, ]  <- xpred_temporaire[Nb_sous_pas+1,]
    ypred[t, ]  <- ypred_temporaire[Nb_sous_pas+1,]
    zpred[t, ]  <- zpred_temporaire[Nb_sous_pas+1,]
    
    theta_proie <- theta_proie_temporaire[Nb_sous_pas + 1, ]
    theta_pred  <- theta_pred_temporaire[Nb_sous_pas + 1, ]
    
    phi_proie <- phi_proie_temporaire[Nb_sous_pas + 1, ]
    phi_pred  <- phi_pred_temporaire[Nb_sous_pas + 1, ]
    
    capture_totale       <- c(capture_totale, capture)
    Estomac_1er  <- c(Estomac_1er, estomac[1]) #Juste pour vérification dans le plot
    print(t)
  }
  return(list(xproie = xproie, yproie = yproie, zproie = zproie, xpred  = xpred,  ypred  = ypred, zpred = zpred, capture_totale = capture_totale, Estomac_1er = Estomac_1er))
}




res    <- simul(N, P, T)

xproie <- res$xproie
yproie <- res$yproie
zproie <- res$zproie

xpred  <- res$xpred
ypred  <- res$ypred
zpred  <- res$zpred


capture_totale   <- res$capture_totale
Estomac_1er   <- res$Estomac_1er




library(rgl)
library(htmlwidgets)

options(rgl.useNULL = TRUE)
open3d(windowRect = c(0, 0, 1200, 900))

view3d(theta = 45, phi = 15, fov = 150, zoom = 0.8)
aspect3d(1, 1, 1)

# Éclairage pour donner du relief
light3d(theta = 45, phi = 30, viewpoint.rel = TRUE)
light3d(theta = 135, phi = 20, viewpoint.rel = TRUE)


dim_tore <- c(torex, torey, torez)
ref_tore <- min(dim_tore)

rayon_proie <- 0.020 * ref_tore
rayon_pred  <- 0.050 * ref_tore


bg3d(color = "#1a5fa8")

static_ids <- c(
  axes3d(edges = "bbox", col = "white", labels = FALSE, tick = FALSE, box = TRUE),
  grid3d(c("x", "y", "z"), col = adjustcolor("white", alpha.f = 0.25), lwd = 0.4)
)


labels <- paste0("t = ", seq_len(T), " | Captures totales = ", capture_totale)


frame_ids <- vector("list", T)

for (t in seq_len(T)) {
  ids_t <- integer(0)
  
  # Proies vivantes
  idx_vivantes <- which(!is.na(xproie[t, ]))
  if (length(idx_vivantes) > 0) {
    id_proies <- spheres3d(
      x      = xproie[t, idx_vivantes],
      y      = yproie[t, idx_vivantes],
      z      = zproie[t, idx_vivantes],
      radius = rayon_proie,
      color  = "#FFE033",
      alpha  = 0.85
    )
    ids_t <- c(ids_t, id_proies)
  }
  
  # Prédateurs
  if (nrow(xpred) >= t && ncol(xpred) > 0) {
    id_pred <- spheres3d(
      x      = xpred[t, ],
      y      = ypred[t, ],
      z      = zpred[t, ],
      radius = rayon_pred,
      color  = "#FF3C3C",
      alpha  = 0.95
    )
    ids_t <- c(ids_t, id_pred)
  }
  
  # Texte affiché dans la scène
  id_txt <- texts3d(
    x = -0.95 * torex,
    y =  0.95 * torey,
    z =  0.95 * torez,
    texts = labels[t],
    color = "white",
    cex = 1.2,
    font = 2
  )
  
  frame_ids[[t]] <- c(static_ids, ids_t, id_txt)
}

# ─── Construction du widget animé ────────────────────────────────────────
w <- rglwidget(width = 1200, height = 900) |>
  playwidget(
    subsetControl(
      value   = 1,
      subsets = frame_ids,
      
    ),
    start      = 1,
    stop       = T,
    step       = 1,
    rate       = 6,
    loop       = TRUE,
    components = c("Reverse", "Play", "Slower", "Faster", "Reset", "Slider", "Label")
  )

# ─── Export HTML autonome ─────────────────────────────────────────────────
htmlwidgets::saveWidget(w, "simulation_banc.html", selfcontained = TRUE)
message("Exporté : simulation_banc.html")
browseURL("simulation_banc.html")
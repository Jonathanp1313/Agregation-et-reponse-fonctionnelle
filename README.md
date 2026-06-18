# Agregation-et-reponse-fonctionnelle
Codes permettant de tester 4 configurations proie-prédateur d'agrégation et de visualiser la réponse fonctionnelle émergente.

- opti_quadruple_3D.c   est le code principal.
  Il permet de générer 4 fichiers txt qui recencent les taux de capture pour chaque N (nombre de proie) allant de 1 à Nmax, et ce pour les 4 configurations
  Les noms des fichiers générés sont :
  solitaires.txt
  proie_banc.txt
  pred_banc.txt
  banc_banc.txt

  
- opti_quadruple_3D.R
  permet de générer une simulation animée de ce qui se passe pour un N spécifique. Cela permet notamment de vérifier si les valeurs des paramètres choisies ne donnent pas lieu à des comportements absurdes des individus.


- quadruple_rep_fonctionnelle.R
  permet de tracer simplement les nuages de points des données des 4 fichiers générées par opti_quadruple_3D.c


- quadruple_rep_fonctionnelle_new_regression.R
  permet de tracer les nuages de points des données des 4 fichier générées par opti_quadruple_3D.c + faire une régression paramétrique.
  La régression concerne : Holling II (   g(N) = aN / (1+bN)​   ) ; Holling III (   g(N)= aN² / (1+bN²)   ) ; Holling IV (   g(N)=aN / (1+bN+cN²)  )
  L'optimisation de la régression et basée sur la minimisation du RMSE.
  Il faut upload le fichier puis appuyer sur le bouton "Ajuster (minimiser RMSE)"

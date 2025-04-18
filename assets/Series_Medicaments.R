# Charger les packages nécessaires
install.packages("forecast")
install.packages("broom")
library(forecast)
library(ggplot2)
library(dplyr)
library(broom)

# Charger les données depuis le fichier CSV
df <- read.csv("data-_4_.csv")

# Convertir les colonnes en format Date et vérifier les noms de colonnes
# Supposons que MONTH est au format "YYYY-MM"
# Convertir MONTH en une date en ajoutant le jour "01" pour chaque mois
df$DATE <- as.Date(paste0(df$MONTH, "-01"))

# Extraire les années et les mois dans des colonnes distinctes si besoin
df$ANNEE <- format(df$DATE, "%Y")
df$MOIS <- format(df$DATE, "%m")


# Définir les dates de début et de fin
start_date <- as.Date("2015-01-01")
end_date <- as.Date("2022-03-31")
df <- df %>% filter(DATE >= start_date & DATE <= end_date)

# Créer la série temporelle mensuelle
ts_data <- ts(df$INDICATEUR, start = c(2015, 1), end = c(2022, 3), frequency = 12)

# Définir la date d'interruption (par exemple, mars 2020)
interruption_date <- as.Date("2020-03-01")
interruption_index <- which(df$DATE == interruption_date)

# Créer la variable d'intervention
intervention <- rep(0, length(ts_data))
intervention[interruption_index:length(ts_data)] <- 1

# Ajuster un modèle ARIMA avec intervention
model_arima_with_intervention <- Arima(ts_data, xreg = intervention, order = c(1, 0, 1))

# Prévisions avec intervention
forecast_with_intervention <- forecast(model_arima_with_intervention, xreg = intervention)

# Ajouter les prévisions au dataframe pour faciliter l'affichage
df$forecast <- c(fitted(model_arima_with_intervention), rep(NA, length(forecast_with_intervention$mean) - length(fitted(model_arima_with_intervention))))

# Créer le graphique avec les dates d'origine et date d'interruption
graph4 <- ggplot(df, aes(x = DATE, y = INDICATEUR)) +
  geom_line(color = "blue") +  # Ligne continue pour la série temporelle
  geom_line(aes(y = forecast), color = "darkgreen", linetype = "dotted") +  # Prévisions
  geom_vline(xintercept = as.numeric(interruption_date), linetype = "dashed", color = "red") +  # Ligne d'interruption
  labs(
    title = "Modèle ARIMA avec Intervention sur Plaquenil",
    subtitle = paste("Interruption détectée à la date :", interruption_date),
    x = "Date",
    y = "Indicateur"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +  # Paramètres pour afficher les dates
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotation des labels de date

# Afficher le graphique
print(graph3)

# Séparer les données en périodes avant et après l'interruption
df_before <- df %>% filter(DATE < interruption_date)
df_after <- df %>% filter(DATE >= interruption_date)

# Calculer les moyennes avant et après l'interruption
mean_before <- mean(df_before$INDICATEUR, na.rm = TRUE)
mean_after <- mean(df_after$INDICATEUR, na.rm = TRUE)

# Test statistique de comparaison de moyennes
mean_test <- t.test(df_before$INDICATEUR, df_after$INDICATEUR, var.equal = FALSE)
print(mean_test)

# Calculer les pentes avant et après l'interruption en ajustant des modèles linéaires simples
model_before <- lm(INDICATEUR ~ as.numeric(DATE), data = df_before)
slope_before <- coef(model_before)[2]

model_after <- lm(INDICATEUR ~ as.numeric(DATE), data = df_after)
slope_after <- coef(model_after)[2]

# Afficher les pentes avant et après
cat("Pente avant l'interruption:", slope_before, "\n")
cat("Pente après l'interruption:", slope_after, "\n")

# Comparaison statistique des pentes avec interaction
df <- df %>%
  mutate(
    period = ifelse(DATE < interruption_date, "Avant", "Après"),
    time_numeric = as.numeric(DATE)
  )

model_full <- lm(INDICATEUR ~ time_numeric * period, data = df)
summary_full <- summary(model_full)

# Résultat final
cat("Comparaison statistique des pentes :\n")
print(summary_full)


# Comparaison des moyennes avant et après l'interruption
mean_before2 <- mean(df_before$INDICATEUR, na.rm = TRUE)
mean_after2 <- mean(df_after$INDICATEUR, na.rm = TRUE)

# Test de Welch pour comparer les moyennes avant et après
mean_test <- t.test(df_before$INDICATEUR, df_after$INDICATEUR, var.equal = FALSE)
cat("Comparaison des moyennes avant et après l'interruption:\n")
print(mean_test)

# Calcul de la différence de pente avant et après l'interruption
model_before2 <- lm(INDICATEUR ~ as.numeric(DATE), data = df_before)
slope_before2 <- coef(model_before)[2]

model_after2 <- lm(INDICATEUR ~ as.numeric(DATE), data = df_after)
slope_after2 <- coef(model_after)[2]

# Comparaison des pentes avant et après l'interruption
cat("Pente avant l'interruption:", slope_before2, "\n")
cat("Pente après l'interruption:", slope_after2, "\n")

# Comparaison statistique des pentes avec un modèle de régression incluant l'interaction
df <- df %>%
  mutate(
    period = ifelse(DATE < interruption_date, "Avant", "Après"),
    time_numeric = as.numeric(DATE)
  )

model_full2 <- lm(INDICATEUR ~ time_numeric * period, data = df)
summary_full2 <- summary(model_full2)

cat("Analyse de la différence de pente (interaction) entre les périodes Avant et Après :\n")
print(summary_full2)


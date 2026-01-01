install.packages("BatchGetSymbols")
#---- Chargement des Packages ------
library(yfR)                 
library(xts)                 
library(zoo)                 
library(moments)             
library(FinTS)               
library(tseries)             
library(scales)              
library(urca)                
library(CADFtest)           
library(forecast)            
library(lmtest)              
library(TSA)                 
library(ghyp)               
library(PerformanceAnalytics) 
library(fGarch)              
library(rugarch)  
library(BatchGetSymbols)

#---- Initialisation RT/RTE/RTT  ---- 
my_ticker <- 'RIG'
first_date <- "2015-01-04"
last_date <- "2025-11-12"

df_yf <- yf_get(tickers = my_ticker, 
                first_date = first_date,
                last_date = last_date,
                freq_data = 'daily', 
                type_return = 'log')

pt <- df_yf$price_adjusted
dates <- df_yf$ref_date[-1]
rt <- df_yf$ret_adjusted_prices[-1]

N <- length(rt)
rte <- rt[1:1762]   
T_est <- length(rte)
rtt <- rt[1763:N]   

#---- Calcul des VaR -----
var_normale <- VaR(rte, p = 0.95, method = "gaussian")
var_cf      <- VaR(rte, p = 0.95, method = "modified")

#VaR Historique 

var_hist_static <- VaR(rte, p = 0.95, method = "historical")
print(paste("VaR Historique (statique sur rte) :", var_hist_static))



window_size <- 500  # Taille de l'historique glissant
var_hist_series <- numeric(length(rtt)) 

for (i in 1:length(rtt)) {
  idx_end <- T_est + i - 1
  idx_start <- idx_end - window_size + 1
  history_window <- rt[idx_start:idx_end]
  

  var_hist_series[i] <- VaR(history_window, p = 0.95, method = "historical")
}

dates_rtt <- dates[(T_est + 1):N]
var_hist_xts <- xts(var_hist_series, order.by = dates_rtt)

plot(var_hist_xts, main = "VaR Historique (Rolling Window)", col = "blue", lwd = 2)

#---- Prévision de la VaR Paramétrique (GARCH) - Fiche 4 ----

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)

spec_egarch = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(1,1)), 
                         distribution.model="ghyp")

roll = ugarchroll(spec_egarch, 
                  data = rt, 
                  n.ahead = 1, 
                  forecast.length = length(rtt), 
                  refit.every = 30, 
                  refit.window = "moving", 
                  solver = "hybrid", 
                  cluster = cl, 
                  fit.control = list(), 
                  calculate.VaR = TRUE, 
                  VaR.alpha = 0.05, 
                  keep.coef = TRUE)

stopCluster(cl)

dates_rtt_check <- dates[(T_est + 1):length(dates)]

valueatrisk_zoo <- zoo(roll@forecast$VaR[,1], order.by = dates_rtt_check)
reelles_zoo <- zoo(roll@forecast$VaR[,2], order.by = dates_rtt_check)

plot(reelles_zoo, type='b', col="grey", 
     main="Backtesting : Rendements vs VaR EGARCH (Rolling)", 
     xlab="Dates", ylab="Rendements", pch=20, cex=0.5)

lines(valueatrisk_zoo, type='l', col="red", lwd=2)

legend("bottomleft", legend=c("Rendements Réels", "VaR EGARCH 95%"), 
       col=c("grey", "red"), lty=c(1,1), lwd=c(1,2))

viols <- sum(reelles_zoo < valueatrisk_zoo)
total <- length(reelles_zoo)
print(paste("Nombre de violations :", viols, "sur", total))
print(paste("Pourcentage :", round(viols/total*100, 2), "%"))


#---- Choix de la Distribution ----
#On teste différentes distributions pour trouver celle qui correspond le mieux à rte. Lors du projet 1, l'étude de rte a révélé que le coefficient de skweness n’est pas significatif (p-value = 0.2679 > 0.05), donc on accepte l'hypothèse nulle du skew test : la
#skewness est statistiquement nulle, la distribution est symmétrique : il n’y a pas d’asymétrie statistiquement prouvée entre la probabilité des gains et celle des pertes.
#On estime donc uniquement les distributions symétriques. 
#Pour avoir une première idée de la distribution, on réalise un QQplot en supposant la normalité.

qqnorm(rte)
qqline(rte, col = 2)

#On voit bien sur le graphique à quel point les queues de distribution sont éloignées de la normale (en rouge).
# Les 2 queues de la distribution sont plus épaisses qu’une loi normale et semblent être aussi lourdes : la distance entre la courbe
#et la droite semble être autant importante pour les valeurs positives que négatives de rte. 


#1) Distribution Gausienne/Normale 

fitn<-fit.gaussuv(data=rte)
summary(fitn)

#2) Distribution student Symmétrique 

fitstu<-fit.tuv(rte,silent=T, symmetric = TRUE)
summary(fitstu)

#3) Distribution Hyperbolique symmétrique 
fithyp<-fit.hypuv(rte,silent=T, symmetric = TRUE)
summary(fithyp)

#4) Distribution hyperbolique généralisée symmétrique 
fitghypuv<-fit.ghypuv(rte,silent=T, symmetric = TRUE)
summary(fitghypuv)


#On sait que l'AIC est un critère qui a tendance à pénaliser les modèles les plus complexes (ie avec le plus de variables)
#Or, ici, le nombre de paramètres est égal à 3 pour chacun des tests de distribution. On peut donc utiliser l'AIC comme critère pour départager nos distributions. 
#Au sens de la minimisation de l'AIC, la distribution la plus adapéte à rte est la distribution de student symmétrique (AIC = -6101.57).  

#Vérification Graphique 

op <- par(mfrow = c(1,2))
plot(ghstfit, which = 1)
plot(ghstfit, which = 3)
par(op)

#L'analyse graphique de l'histogramme de rte suggère qu'une distribution symmétrique est effectivement la plus adaptée à rte. 

plot(density(rte))
lines(fitstu,col=2)#student
lines(fithyp,col=3)
lines(fitghypuv,col=4)
legend("topleft",legend =c("rte","student","hyp","ghyp"), col =1:4,lty=rep(1,4))

#L'analyse graphique de la distribution de rte nous montre que les plus proches de celle-ci sont les distributions student et ghyp ce que l'on avait déja pu comprendre en regardant l'AIC des tests.


# Modèles ARCH, GARCH symétriques 
#---- Avec distribution normale (Pédagogique) ----
spec1 = ugarchspec()
fit1 = ugarchfit(spec = spec1, data = rt,out.sample=length(rtt))
show(fit1)
#alpha1 et beta1 significatifs /Absence d'autocorrélation / Absence de clusters de volatilité / On ne rejette pas l(hypothèse nulle de stabilité /Effet taille mais pas effet signe / Tous les coeffs sont <0.05, normale pas adaptée logique elle est mesocurtique alors que rte leptokurtique)

#---- ARMA(1,1) + GARCH(1,1) 3157.359  ----
spec2 = ugarchspec(distribution.model="ghyp")
fit2 = ugarchfit(spec = spec1, data = rt,out.sample=length(rtt))
show(fit2)

#alpha1/beta1/ghlambda significatifs/Absence d'autocorrélation/ Absence de clusters de volatilité / On ne rejette pas l'hypothèse nuelle de stabilité/ Effet Taille des 2 chcos mais pas d'effet signe/ Tous les coeffs de Pearson sont >0.05 la distribution ghyp est adaptée à notre série rte. 

#Prise en compte des Effets Weekend/Janvier
jour=format(dates, format = "%A")
mois=format(dates, format = "%B")
moisrte=mois[1:1762]
juin=as.integer(moisrte=="juin")
jourrte=jour[1:1762]#comme rte
lundi=as.integer(jourrte=="lundi")
spec2bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
fit2bis = ugarchfit(spec = spec2bis, data = rt,out.sample=length(rtt))
show(fit2bis)

#omega/alpha1/beta1/ghlambda significatifs/Absence d'autocorrélation/ Absence de clusters de volatilité / On ne rejette pas l'hypothèse nulle de stabilité / Effet Taille maos pas effet signe / tous les coeffs sont >0.05

nigarch=newsimpact(z=NULL, fit1)
plot(nigarch$zx,nigarch$zy, xlab=nigarch$xexpr,ylab=nigarch$yexpr ,type="l", main = "Courbe des impacts")

#On a bien le même impact dans un GARCH des bonnes et des mauvaises nouvelles sur la volatilité 

prev = ugarchforecast(fit2, n.ahead=length(rtt),n.roll=length(rtt))
plot(prev,which="all")





#---- ARCH-M 3156.867 ----
spec3 = ugarchspec(mean.model=list(armaOrder=c(1,1),archm=TRUE),distribution.model="ghyp")
fit3 = ugarchfit(spec = spec3,data = rt,out.sample=length(rtt),solver="hybrid")
fit3
#spec3bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
#fit3bis = ugarchfit(spec = spec3bis, data = rt,out.sample=length(rtt))
#show(fit3bis)
#alpha1, beta1, ghlambda significatifs / Absence d'autocorrélation / Absence de clusters de volatilité / On ne rejette pas l'hypothèse nulle de stabilité / effet Taille mais pas d'effet signe / Toutes les p-val > 0.05 ghyp est adaptée.
#archm n'est pas significatif donc ce modèle n'est pas adéquat pour nos données. 

#---- IGARCH 3154.709 ----
spec4 = ugarchspec(variance.model=list(model="iGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(1,1)),distribution.model="ghyp")
fit4 = ugarchfit(spec = spec4, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit4)
#spec4bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
#fit4bis = ugarchfit(spec = spec4bis, data = rt,out.sample=length(rtt))
#show(fit4bis)

#omega, alpha1, ghlambda significatifs / absence d'autocorrélation / Absence de clusters de volatilité / On rejette l'hypothèse nulle de stabilité : au moins un coeff pas stable dans le temps mais on ne le voit pas ????? / Effet Taille mais pas signe / Tous les coeffs > 0.05 donc distribution adaptée


# Modèles ARCH, GARCH asymétriques
#---- Modèle EGARCH 3163.537 ----
spec5 = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                     mean.model=list(armaOrder=c(1,1)),distribution.model="ghyp")
fit5 = ugarchfit(spec = spec5, data = rt,out.sample=length(rtt),solver="hybrid")
show(fit5)
#spec5bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
#fit5bis = ugarchfit(spec = spec5bis, data = rt,out.sample=length(rtt))
#show(fit5bis)

#mu, skew, shape non significatifs / Absence d'autocorrélation / Absence de clusters de volatilité / On ne rejette pas H0 / Effet Taille mais pas signe / Toutes les p-val > 0.05 distribution adaptée

#spec5B = ugarchspec(variance.model=list(model="eGARCH", garchOrder=c(1,1)),
                       #mean.model=list(armaOrder=c(1,1),include.mean=F),distribution.model="ghyp")
#fit5B = ugarchfit(spec = spec5B, data = rt,out.sample=length(rtt),solver="hybrid")
#show(fit5B)

#---- Modèle GJR-GARCH 3161.784 ----

spec6 = ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(1,1)),distribution.model="ghyp")
fit6= ugarchfit(spec = spec6,data = rt,out.sample=length(rtt),solver="hybrid")
show(fit6)
#spec6bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
#fit6bis = ugarchfit(spec = spec6bis, data = rt,out.sample=length(rtt))
#show(fit6bis)

# alpha1, beta1, gamma1 significatifs / Absence d'autocorrélation / Absence de clusters de volatilité / On ne rejette pas l'hypothèse nulle de stabilité / Effet Taille mais pas d'effet signe / toutes les p-val > 0.05 donc ghyp adaptée 
# Gamma significatif donc ce modèle prend en compte l'effet de levier 
#---- Modèle APARCH 3164.102 ----
spec7 = ugarchspec(variance.model=list(model="apARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(1,1)),distribution.model="ghyp")
fit7= ugarchfit(spec = spec7,data = rt,out.sample=length(rtt),solver="hybrid")
show(fit7)
#spec7bis = ugarchspec(distribution.model ='ghyp',mean.model=list(external.regressors=as.matrix(cbind(lundi,juin))))
#fit7bis = ugarchfit(spec = spec7bis, data = rt,out.sample=length(rtt))
#show(fit7bis)

#alpha1, beta1, ghlambda, gamma1, delta, ar1, ma1 significatifs / Absence d'autocorrélation / Absence de clusters de volatilité / On ne rejette pas l'hypothèse nulle de stabilité / EFFEt taille mais pas d'effet signe / Toute sles p-val > 0.05 donc la distriution est la bonne. 


#On utilise le SBIC pour départager nos modèles et on choisit dont le modèle EGARCH(1,1), SBIC = -3.5796

#---- Backtesting des VaRs ----

actual_ret <- as.numeric(rtt)
n_test <- length(actual_ret)

var_norm_vec <- rep(as.numeric(var_normale), n_test)
var_cf_vec   <- rep(as.numeric(var_cf), n_test)
var_hist_vec <- as.numeric(var_hist_xts)
var_egarch_vec <- as.numeric(valueatrisk_zoo)

modeles <- list(
  "VaR Normale"        = var_norm_vec,
  "VaR Cornish-Fisher" = var_cf_vec,
  "VaR Historique"     = var_hist_vec,
  "VaR EGARCH"         = var_egarch_vec
)

df_backtest <- data.frame()

for (nom in names(modeles)) {
  
  test <- VaRTest(alpha = 0.05, 
                  actual = actual_ret, 
                  VaR = modeles[[nom]], 
                  conf.level = 0.95)
  
  ligne <- data.frame(
    Modele = nom,
    Kupiec_Pvalue = round(test$uc.LRp, 4),
    Decision_Kupiec = ifelse(test$uc.LRp > 0.05, "Validé", "Rejeté"),
    Christoffersen_Pvalue = round(test$cc.LRp, 4),
    Decision_Christoffersen = ifelse(test$cc.LRp > 0.05, "Validé", "Rejeté"),
    Violations_Est = test$actual.exceed,
    Violations_Theo = round(test$expected.exceed, 1),
    Taux_Est = paste0(round((test$actual.exceed / n_test) * 100, 2), "%"),
    Taux_Theo = "5%"
  )
  
  df_backtest <- rbind(df_backtest, ligne)
}

print(df_backtest)

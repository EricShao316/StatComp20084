#' @title bootstrap and jackknife
#' @name bootstrap.and.jackknife
#' @useDynLib StatComp20084
#' @examples 
#' \dontrun{
#' library(bootstrap)
#' law = law82[sample(1:nrow(law82), 15, replace = FALSE), c(2, 3)]
#' law = as.matrix(law)
#' # Jackknife estimation of bias and stand error
#' n = nrow(law)
#' jack_i = numeric(n)
#' theta_hat = cor(law[,1], law[,2])
#' for(i in 1:n){
#'   jack_i[i] = cor(law[-i,1], law[-i,2]) 
#' }
#' jack_mean = mean(jack_i)
#' jack_bias = (n - 1) * (jack_mean - theta_hat)
#' jack_sd = sqrt((n - 1) * mean((jack_i - jack_mean) ^ 2))
#' list(jack_bias = jack_bias, jack_sd = jack_sd)
#' }
#' \dontrun{
#' library(boot)
#' air = as.matrix(aircondit)
#' plot(density(air), ylab = " ", main = "Density of air-condition data")
#' ci.norm = ci.basic = ci.perc = ci.bca = matrix(NA, 1, 2)
#' me = boot(data = air, statistic = boot.mean, R = 1e3)
#' ci = boot.ci(me, type = c("norm", "basic", "perc", "bca"))
#' list(mean_hat = mean(air),
#'      norm = ci$norm[2:3],
#'      basic = ci$basic[4:5],
#’      perc = ci$percent[4:5],
#‘     BCa = ci$bca[4:5])
#' }
#' \dontrun{
#' library(bootstrap)
#' sc = as.matrix(scor)
#' n = nrow(sc)
#' cov_sc = cov(sc)
#' eicov = eigen(cov_sc)$values
#' theta_hat = max(eicov) / sum(eicov)
#' theta_jack = numeric(n)
#' for(i in 1:n){
#'   jack_cov = cov(sc[-i,])
#'   ei = eigen(jack_cov)$values
#'   theta_jack[i] = max(ei) / sum(ei)
#' }
#' bias_jack = (n - 1) * (mean(theta_jack) - theta_hat)
#' sd_jack = sqrt((n - 1) ^ 2 / n * var(theta_jack))
#' list(bias_jack = bias_jack, sd_jack = sd_jack)
#' }
#' \dontrun{
#' library(DAAG)
#' attach(ironslag)
#' n = length(magnetic)
#' # Leave-two-out cross validation
#' # Two options
#' fly = numeric(2)
#' for(i in 1:(n-1)){
#'   j = i + 1
#'   while(j <= n){
#'     fly = rbind(fly, c(i, j))
#'     j = j + 1
#'   }
#' }
#' fly = as.matrix(fly[-1,])
#' m = nrow(fly)
#' e1 = e2 = e3 = e4 = matrix(NA, m, 2)
#' for(i in 1:m){
#'   y = magnetic[-fly[i,]]
#'   x = chemical[-fly[i,]]
#'   J1 <- lm(y ~ x)
#'   yhat1 <- J1$coef[1] + J1$coef[2] * chemical[fly[i,]] 
#'   e1[i,] <- magnetic[fly[i,]] - yhat1 
#'   J2 <- lm(y ~ x + I(x^2))
#'   yhat2 <- J2$coef[1] + J2$coef[2] * chemical[fly[i,]] + 
#'     J2$coef[3] * chemical[fly[i,]]^2
#'   e2[i,] <- magnetic[fly[i,]] - yhat2 
#'   J3 <- lm(log(y) ~ x)
#'   logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[fly[i,]] 
#'   yhat3 <- exp(logyhat3)
#'   e3[i,] <- magnetic[fly[i,]] - yhat3 
#'   J4 <- lm(log(y) ~ log(x))
#'   logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[fly[i,]]) 
#'   yhat4 <- exp(logyhat4)
#'   e4[i,] <- magnetic[fly[i,]] - yhat4
#' }
#' list(mse1 = mean(e1^2), 
#'      mse2 = mean(e2^2),
#'      mse3 = mean(e3^2),
#'      mse4 = mean(e4^2))
#' L2 = lm(magnetic ~ chemical + I(chemical^2))
#' L2
#' }
NULL



#' @title boot mean
#' @description using to calculate bootstrap estimate
#' @param x,i vector and the chosen index
#' @export
boot.mean = function(x, i) mean(x[i])
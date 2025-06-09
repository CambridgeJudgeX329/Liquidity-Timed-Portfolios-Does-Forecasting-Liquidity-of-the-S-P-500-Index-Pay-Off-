###############################################################################
# MASTER SCRIPT — Liquidity-Timing Project
# ---------------------------------------
# A) Two-phase ARMA–(G)ARCH search for Δlog-ILLIQ
# B) Forecast-quality comparison:  LOCF  vs  ARIMA  vs  ARIMAX
# C) Liquidity-managed portfolios
#    C-1  uses *observed* ILLIQ  (upper-bound strategy)
#    C-2  uses *forecast* ILLIQ  (implementable strategy)
# D) Final figures & tables (1–4) plus regression output
#
###############################################################################

##################################  PACKAGES  #################################
# (Uncomment next line if any package is missing)
# install.packages(c("arrow","data.table","forecast","zoo","ggplot2","Metrics",
#   "xts","rugarch","tseries","FinTS","reshape2","quantmod",
#   "PerformanceAnalytics","patchwork","AER","TSA","dynlm","car",
#   "lmtest","sandwich","urca","gridExtra","grid"))

library(arrow);  library(data.table);  library(forecast);     library(zoo)
library(ggplot2);library(Metrics);     library(xts);          library(rugarch)
library(tseries);library(FinTS);       library(reshape2);     library(quantmod)
library(PerformanceAnalytics);         library(patchwork);    library(AER)
library(TSA);     library(dynlm);      library(car);          library(lmtest)
library(sandwich);library(urca);       library(gridExtra);    library(grid)

###############################################################################
# SECTION A —— ARMA–(G)ARCH MODEL SELECTION FOR Δlog-ILLIQ
###############################################################################
rm(list = ls()); graphics.off()

# A-1  Load OHLCV & construct Amihud ILLIQ ....................................
ohlcv <- as.data.table(read_parquet("sp500_ohlcv_trunc.parquet"))
setnames(ohlcv, names(ohlcv), sub("^\\('([^']+)',.*", "\\1", names(ohlcv)))
ohlcv[, ret        := log(Close) - log(shift(Close))]
ohlcv[, Volume     := as.numeric(Volume)]
ohlcv[, dollar_vol := Volume * Close]
ohlcv[, illiq      := 1e6 * abs(ret) / pmax(dollar_vol, 1)]
ohlcv              <- ohlcv[!is.na(ret) & illiq > 0]
ohlcv[, log_illiq  := log(illiq)]
setorder(ohlcv, Date)

# A-2  Δlog-ILLIQ & mean ARMA order ............................................
ohlcv[, dlog_illiq := log_illiq - shift(log_illiq)]
full_delta <- na.omit(ohlcv$dlog_illiq)
base_arima <- auto.arima(full_delta, stationary = TRUE,
                         seasonal = FALSE, max.p = 15, max.q = 15)
opt_p <- base_arima$arma[1];  opt_q <- base_arima$arma[2]

# A-3  Phase-1 grid: (p,q) for sGARCH-Normal ..................................
phase1 <- data.table(p = integer(), q = integer(), AIC = numeric())
for(p in 1:3) for(q in 1:3){
  sp <- ugarchspec(variance.model=list(model="sGARCH",
                                       garchOrder=c(p,q)),
                   mean.model=list(armaOrder=c(opt_p,opt_q),include.mean=TRUE),
                   distribution.model="norm")
  ft <- try(ugarchfit(spec=sp, data=full_delta, solver="hybrid"), silent=TRUE)
  phase1 <- rbind(phase1,
                  data.table(p=p,q=q,
                             AIC=if(inherits(ft,"try-error")) NA else
                               infocriteria(ft)[1]))
}
setorder(phase1,AIC); best_p<-phase1$p[1]; best_q<-phase1$q[1]

# A-4  Phase-2: hold (p,q) fixed, vary variance & distribution .................
var_mods <- c("sGARCH","eGARCH","gjrGARCH","apARCH")
dists    <- c("norm","ged","sstd","sged","nig","jsu")
phase2 <- data.table(Model=character(),Dist=character(),
                     AIC=numeric(), JB=numeric(), ARCH=numeric())
fits2  <- list()
for(vm in var_mods) for(ds in dists){
  sp <- ugarchspec(variance.model=list(model=vm,
                                       garchOrder=c(best_p,best_q)),
                   mean.model=list(armaOrder=c(opt_p,opt_q),include.mean=TRUE),
                   distribution.model=ds)
  nm <- paste0(vm,"(",best_p,",",best_q,")-",ds)
  ft <- try(ugarchfit(spec=sp, data=full_delta, solver="hybrid"), silent=TRUE)
  if(!inherits(ft,"try-error")){
    rz <- residuals(ft, standardize=TRUE)
    phase2 <- rbind(phase2,
                    data.table(Model=vm,Dist=ds,AIC=infocriteria(ft)[1],
                               JB   = jarque.bera.test(rz)$p.value,
                               ARCH = ArchTest(rz,lags=10)$p.value))
    fits2[[nm]] <- ft
  } else phase2 <- rbind(phase2,
                         data.table(Model=vm,Dist=ds,AIC=NA,JB=NA,ARCH=NA))
}
setorder(phase2,AIC)
best2 <- if(nrow(phase2[JB>0.05 & ARCH>0.05])>0)
  phase2[JB>0.05 & ARCH>0.05][1] else phase2[1]
best_vm<-best2$Model; best_dist<-best2$Dist

# A-5  Phase-3: optimise (p,q) for best variance/dist ..........................
phase3 <- data.table(p=integer(),q=integer(),AIC=numeric())
for(p in 1:3) for(q in 1:3){
  sp <- ugarchspec(variance.model=list(model=best_vm,
                                       garchOrder=c(p,q)),
                   mean.model=list(armaOrder=c(opt_p,opt_q),include.mean=TRUE),
                   distribution.model=best_dist)
  ft <- try(ugarchfit(spec=sp,data=full_delta,solver="hybrid"), silent=TRUE)
  phase3 <- rbind(phase3,
                  data.table(p=p,q=q,
                             AIC=if(inherits(ft,"try-error")) NA else
                               infocriteria(ft)[1]))
}
setorder(phase3,AIC); best3_p<-phase3$p[1]; best3_q<-phase3$q[1]

final_spec <- ugarchspec(
  variance.model=list(model=best_vm,garchOrder=c(best3_p,best3_q)),
  mean.model=list(armaOrder=c(opt_p,opt_q),include.mean=TRUE),
  distribution.model=best_dist)
final_fit <- ugarchfit(spec=final_spec, data=full_delta, solver="hybrid")

cat("\n===== Final GARCH Diagnostics =====\n",
    "Model: ",best_vm,"(",best3_p,",",best3_q,")-",best_dist,"\n",
    "AIC  : ", round(infocriteria(final_fit)[1],4),"\n",
    "JB p : ", signif(jarque.bera.test(residuals(final_fit,standardize=TRUE))$p.value,4),"\n",
    "ARCH p:", signif(ArchTest(residuals(final_fit,standardize=TRUE),lags=10)$p.value,4),"\n")

###############################################################################
# SECTION B —— Forecast-Error Comparison (LOCF  vs  ARIMA  vs  ARIMAX)
###############################################################################
###############################################################################

# 0 … packages already loaded above

## 1 ── Re-load data (keeps script stand-alone) --------------------------------
ohlcv <- as.data.table(read_parquet("sp500_ohlcv_trunc.parquet"))
setnames(ohlcv, names(ohlcv), sub("^\\('([^']+)',.*", "\\1", names(ohlcv)))
ohlcv[, ret        := log(Close) - log(shift(Close))]
ohlcv[, Volume     := as.numeric(Volume)]
ohlcv[, dollar_vol := Volume * Close]
ohlcv[, illiq      := 1e6 * abs(ret)/pmax(dollar_vol,1)]
ohlcv              <- ohlcv[!is.na(ret) & illiq>0]
ohlcv[, log_illiq  := log(illiq)]
ohlcv[, dlog_illiq := log_illiq - shift(log_illiq)]
ohlcv[, dlog_vol   := log(dollar_vol) - shift(log(dollar_vol))]
setorder(ohlcv, Date)

roll_win <- 252
ohlcv[, `:=`(
  fct_LOCF    = shift(illiq,1,"lag"),
  fct_ARIMA   = NA_real_,
  fct_ARIMAX  = NA_real_)]

full_delta <- na.omit(ohlcv$dlog_illiq)
base_arima <- auto.arima(full_delta, stationary=TRUE, seasonal=FALSE,
                         max.p=5,max.q=5)
p0<-base_arima$arma[1]; q0<-base_arima$arma[2]

for(i in (roll_win+1):nrow(ohlcv)){
  win <- (i-roll_win):(i-1)
  # ── ARIMA -------------------------------------------------------------------
  fit_a <- Arima(na.omit(ohlcv$dlog_illiq[win]),
                 order=c(p0,0,q0), include.constant=TRUE)
  d_hat <- forecast(fit_a,h=1)$mean
  ohlcv$fct_ARIMA[i] <- exp(ohlcv$log_illiq[i-1] + d_hat)
  # ── ARIMAX ------------------------------------------------------------------
  X <- na.omit(ohlcv$dlog_vol[win]);  Y <- na.omit(ohlcv$dlog_illiq[win])[seq_along(X)]
  fit_x <- Arima(Y, order=c(1,0,1), xreg=X, include.constant=TRUE)
  dx_hat <- forecast(fit_x, h=1, xreg = ohlcv$dlog_vol[i-1])$mean
  ohlcv$fct_ARIMAX[i] <- exp(ohlcv$log_illiq[i-1] + dx_hat)
}

metrics <- function(a,f){
  ok <- complete.cases(a,f)
  data.table(RMSE=rmse(a[ok],f[ok]),
             MAE = mae(a[ok],f[ok]),
             Bias=mean(f[ok]-a[ok]),
             Corr=cor(a[ok],f[ok]))}

res_table <- rbindlist(list(
  LOCF   = metrics(ohlcv$illiq, ohlcv$fct_LOCF),
  ARIMA  = metrics(ohlcv$illiq, ohlcv$fct_ARIMA),
  ARIMAX = metrics(ohlcv$illiq, ohlcv$fct_ARIMAX)
), idcol="Model")
print(res_table, digits=4)

# Plots (unchanged) -----------------------------------------------------------
log_actual <- log10(ohlcv$illiq)
plot_one <- function(col,title){
  ggplot(ohlcv,aes(Date))+
    geom_line(aes(y=log_actual),col="black",size=.3)+
    geom_line(aes_string(y=paste0("log10(",col,")")),
              col="red",alpha=.7,size=.3)+
    labs(title=title,y="log10(ILLIQ)",x="")+theme_minimal()}
p1<-plot_one("fct_LOCF","LOCF vs actual")
p2<-plot_one("fct_ARIMA","Rolling ARIMA(1,1) vs actual")
p3<-plot_one("fct_ARIMAX","Rolling ARIMAX(1,1) vs actual")
(p1/p2)|(p3)+plot_layout(guides="collect")&
  theme(legend.position="none")

###############################################################################
# SECTION C —— Liquidity-Weighted Portfolios
###############################################################################
# C-1  uses *observed* ILLIQ   (upper bound)
# C-2  uses *forecast* ILLIQ   (implementable)
###############################################################################

#### C-1 : OBSERVED ILLIQ  #####################################################

# 1 ── Re-load data again (for isolation) --------------------------------------
ohlcv <- as.data.table(read_parquet("sp500_ohlcv_trunc.parquet"))
setnames(ohlcv, names(ohlcv), sub("^\\('([^']+)',.*", "\\1", names(ohlcv)))
ohlcv[, ret        := log(Close) - log(shift(Close))]
ohlcv[, Volume     := as.numeric(Volume)]
ohlcv[, dollar_vol := Volume * Close]
ohlcv[, illiq      := 1e6 * abs(ret)/pmax(dollar_vol,1)]
ohlcv              <- ohlcv[!is.na(ret) & illiq>0]
ohlcv[, log_illiq  := log(illiq)]
setorder(ohlcv, Date)

# Rolling 5-year geometric mean of ILLIQ --------------------------------------
roll_n <- 252
ohlcv[, roll_m_liq     := zoo::rollapplyr(log_illiq, roll_n, mean, fill=NA)]
ohlcv[, roll_mean_illiq:= exp(roll_m_liq)]
ohlcv[, rel_illiq      := illiq / roll_mean_illiq]
ohlcv[, denom_alt      := rel_illiq^2 + rel_illiq]
c_target <- median(ohlcv$denom_alt,na.rm=TRUE)
ohlcv[, raw_alt   := c_target/denom_alt]
ohlcv[, weight_alt:= pmin(pmax(raw_alt,0),1)]

# Risk-free --------------------------------------------------------------------
getSymbols("DTB3", src="FRED")
rf_dt <- data.table(Date=index(DTB3), rf_pct=as.numeric(DTB3$DTB3))
rf_dt[, rf_daily := (rf_pct/100)/252]
ohlcv <- merge(ohlcv, rf_dt[,.(Date,rf_daily)], by="Date", all.x=TRUE)
ohlcv[, rf_daily := zoo::na.locf(rf_daily,na.rm=FALSE)]

ohlcv[, ret_alt  := weight_alt*ret + (1-weight_alt)*rf_daily]
mean_w <- mean(ohlcv$weight_alt,na.rm=TRUE)
ohlcv[, ret_const:= mean_w*ret + (1-mean_w)*rf_daily]

rets_xts <- xts(ohlcv[,.(ret,ret_alt,ret_const,rf_daily)], order.by=ohlcv$Date)
colnames(rets_xts) <- c("S&P500","LiqM","ConstW","RiskFree")

charts.PerformanceSummary(rets_xts[,c("S&P500","LiqM","ConstW","RiskFree")],
                          main="", wealth.index=TRUE)

table.AnnualizedReturns(rets_xts[,c("S&P500","LiqM","ConstW")],
                        Rf=rets_xts[,"RiskFree"], scale=252)

## Table 1 : HC3 robust regression output
reg_dt   <- na.omit(ohlcv[,.(ret,rel_illiq)])
quad_lm  <- lm(ret ~ rel_illiq + I(rel_illiq^2), data=reg_dt)
Table1 <- coeftest(quad_lm, vcov.=vcovHC(quad_lm,type="HC3"))
print(Table1)

## Figure 2:
ohlcv[, ret_next := shift(ret, type="lead")]

ggplot(ohlcv, aes(x = rel_illiq, y = ret)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = "loess", color = "red") +
  scale_x_log10() +
  labs(x="R_t", y="r_t",)

## Figure 3:
## Figure 3 : R_t vs R_{t+1}

## 1.  Create lead( rel_illiq )  ------------------------------------------------
ohlcv[, rel_illiq_next := shift(rel_illiq, type = "lead")]

## 2.  Keep pairs with no NA and (optionally) drop extreme outliers -------------
plot_dt <- na.omit(ohlcv[, .(Date, R_today = rel_illiq,
                             R_next  = rel_illiq_next)])

## 3.  Scatterplot on log–log scale ---------------------------------------------
ggplot(plot_dt, aes(x = R_today, y = R_next)) +
  geom_point(alpha = 0.25, size = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  geom_smooth(method = "loess", colour = "red", se = FALSE) +
  labs(x     = expression(R[t]~"(today)"),
       y     = expression(R[t+1]~"(next day)")) +
  theme_minimal()

#### C-2 : FORECAST ILLIQ  #####################################################

# Same weight kernel, but now R_t uses the one-step-ahead ARIMA forecast.
# Code identical to your second portfolio file; kept verbatim below.
# ------------------------------------------------------------------------------

# (All lines from “FULL SCRIPT: Liquidity-Weighted Portfolio Using
#  Normalised Expected ILLIQ²+ILLIQ Rule” are pasted here unmodified.)

# ─────────────────────────────────────────────────────────────────────────────
# 0  Packages already loaded above
# ---------------------------------------------------------------------------

# 1 …  (read data, build ILLIQ)  — identical to prior block
ohlcv <- as.data.table(read_parquet("sp500_ohlcv_trunc.parquet"))
setnames(ohlcv, names(ohlcv), sub("^\\('([^']+)',.*", "\\1", names(ohlcv)))
ohlcv[, ret        := log(Close) - log(shift(Close))]
ohlcv[, Volume     := as.numeric(Volume)]
ohlcv[, dollar_vol := Volume * Close]
ohlcv[, illiq      := 1e6 * abs(ret) / pmax(dollar_vol,1)]
ohlcv              <- ohlcv[!is.na(ret) & illiq>0]
ohlcv[, log_illiq  := log(illiq)]
setorder(ohlcv,Date)

# 2 …  Stationarity visuals (unchanged—kept for completeness) ----------------
if(!"log_vol"%in%names(ohlcv)) ohlcv[,log_vol:=log(dollar_vol)]
roll_n <- 252
ohlcv[,roll_m_liq:=zoo::rollapplyr(log_illiq,roll_n,mean,fill=NA)]

# 3 …  Merge daily risk-free ---------------------------------------------------
getSymbols("DTB3",src="FRED")
rf_dt <- data.table(Date=index(DTB3),rf_pct=as.numeric(DTB3$DTB3))
rf_dt[,rf_daily:=(rf_pct/100)/252]
ohlcv <- merge(ohlcv, rf_dt[,.(Date,rf_daily)], by="Date", all.x=TRUE)
ohlcv[,rf_daily:=zoo::na.locf(rf_daily,na.rm=FALSE)]

# 4 …  Rolling ARIMA(1,1) forecast of Δlog(ILLIQ) -----------------------------
ohlcv[, dlog_illiq := log_illiq - shift(log_illiq)]
full_series <- na.omit(ohlcv$dlog_illiq)
full_fit <- auto.arima(full_series, stationary=TRUE, seasonal=FALSE,
                       max.p=5, max.q=5)
opt_p <- full_fit$arma[1]; opt_q <- full_fit$arma[2]

roll_window <- 252
ohlcv[,`:=`(forecast_illiq=NA_real_)]
for(i in (roll_window+1):nrow(ohlcv)){
  win <- (i-roll_window):(i-1)
  fit <- Arima(na.omit(ohlcv$dlog_illiq[win]),
               order=c(opt_p,0,opt_q), include.constant=TRUE)
  d_hat <- forecast(fit,h=1)$mean
  ohlcv$forecast_illiq[i] <- exp(ohlcv$log_illiq[i-1] + d_hat)
}

ohlcv[,roll_mean_illiq:=exp(roll_m_liq)]
ohlcv[,rel_illiq      := forecast_illiq/roll_mean_illiq]
ohlcv[,denom_alt      := rel_illiq^2 + rel_illiq]
c_target <- median(ohlcv$denom_alt,na.rm=TRUE)
ohlcv[,raw_alt   := c_target/denom_alt]
ohlcv[,weight_alt:= pmin(pmax(raw_alt,0),1)]
ohlcv[,ret_alt   := weight_alt*ret + (1-weight_alt)*rf_daily]
mean_w <- mean(ohlcv$weight_alt,na.rm=TRUE)
ohlcv[,ret_const := mean_w*ret + (1-mean_w)*rf_daily]

ohlcv <- ohlcv[!is.na(ret_alt)]
rets_xts <- xts(ohlcv[,.(ret,ret_alt,ret_const,rf_daily)], order.by=ohlcv$Date)
colnames(rets_xts) <- c("S&P500","Alt","ConstW","RiskFree")

charts.PerformanceSummary(rets_xts[,c("S&P500","Alt","ConstW","RiskFree")],
                          main="", wealth.index=TRUE)
table.AnnualizedReturns(rets_xts[,c("S&P500","Alt","ConstW")],
                        Rf=rets_xts[,"RiskFree"], scale=252)

###############################################################################
# END OF MASTER SCRIPT  –  sourcing reproduces all results
###############################################################################

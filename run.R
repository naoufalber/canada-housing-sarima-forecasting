# Run from the repository root: Rscript run.R
suppressPackageStartupMessages(library(forecast))
suppressPackageStartupMessages(library(ggplot2))
source("R/manual_forecast.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)
data_path <- if (length(args)) args[1] else "data/QCAR628BIS.csv"
raw <- read.csv(data_path, stringsAsFactors = FALSE)
stopifnot(ncol(raw) == 2)
names(raw) <- c("date", "value")
raw$date <- as.Date(raw$date)
raw$value <- as.numeric(raw$value)
stopifnot(!anyNA(raw), all(raw$value > 0), !anyDuplicated(raw$date))
raw <- raw[order(raw$date), ]
quarter_index <- as.integer(format(raw$date, "%Y")) * 4L +
  (as.integer(format(raw$date, "%m")) - 1L) %/% 3L
stopifnot(all(diff(quarter_index) == 1L))
start_year <- as.integer(format(raw$date[1], "%Y"))
start_quarter <- (as.integer(format(raw$date[1], "%m")) - 1L) %/% 3L + 1L
x <- ts(log(raw$value), start = c(start_year, start_quarter), frequency = 4)
h <- 8L
stopifnot(length(x) > h + 40L)
train <- head(x, -h)
train <- ts(train, start = start(x), frequency = 4)
test <- as.numeric(tail(x, h))

# Model structures are specified before viewing the final eight quarters.
manual_fit <- forecast::Arima(train, order = c(1, 1, 0),
                             seasonal = list(order = c(0, 1, 1), period = 4))
auto_fit <- forecast::auto.arima(train, seasonal = TRUE,
                               stepwise = FALSE, approximation = FALSE)
predictions <- list(
  "SARIMA (1,1,0)(0,1,1)[4]" = forecast::forecast(manual_fit, h = h, level = 95),
  "Auto ARIMA" = forecast::forecast(auto_fit, h = h, level = 95),
  "Seasonal naive" = forecast::snaive(train, h = h, level = 95)
)
dir.create("reports/current", recursive = TRUE, showWarnings = FALSE)
scores <- do.call(rbind, lapply(names(predictions), function(name) {
  prediction <- as.numeric(predictions[[name]]$mean)
  data.frame(model = name, rmse_log = sqrt(mean((test - prediction)^2)),
             mae_log = mean(abs(test - prediction)), n_test = h)
}))
write.csv(scores, "reports/current/holdout_metrics.csv", row.names = FALSE)
print(scores)

# Refit the pre-specified SARIMA on all available observations for future forecasts.
fit <- forecast::Arima(x, order = c(1, 1, 0),
                       seasonal = list(order = c(0, 1, 1), period = 4))
future <- forecast::forecast(fit, h = h, level = 95)
manual <- sarima_forecast_manual(x, ar = unname(coef(fit)["ar1"]),
                                sma = unname(coef(fit)["sma1"]), h = h,
                                sigma = sqrt(fit$sigma2))
comparison <- data.frame(quarter_ahead = seq_len(h),
                         manual_log = manual$xhat,
                         forecast_package_log = as.numeric(future$mean),
                         lower_95_log = as.numeric(future$lower),
                         upper_95_log = as.numeric(future$upper))
write.csv(comparison, "reports/current/forecast_comparison.csv", row.names = FALSE)
# Use package intervals: the educational manual recursion approximates uncertainty
# beyond one seasonal cycle. Forecast values stay on the log-index scale.
plot <- autoplot(future) +
  ggtitle("Canada real housing prices | SARIMA forecast") +
  xlab("Year") + ylab("log(index, 2010 = 100)") + theme_minimal(base_size = 12)
ggsave("reports/current/forecast.png", plot, width = 11, height = 5, dpi = 160)
capture.output(sessionInfo(), file = "reports/current/sessionInfo.txt")
capture.output(list(data_path = data_path, observations = length(x),
                    first_date = raw$date[1], last_date = tail(raw$date, 1),
                    automatic_order = forecast::arimaorder(auto_fit)),
               file = "reports/current/run.txt")

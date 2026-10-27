# Methodology and execution status

## Academic work

The original study implements ARIMA simulation, innovation reconstruction, recursive forecasting, polynomial composition for seasonal models and Box–Jenkins diagnostics. The report discusses SARIMA(1,1,0)(0,1,1)[4] and comparison with another seasonal specification.

The historical script’s `fit_auto` variable is actually a fixed SARIMA specification in one section, although the report describes automatic selection. The new `run.R` explicitly calls `auto.arima()` for that comparison. The historical grid compares models with different differencing orders; their AIC values should not be treated as directly comparable across those transformations.

## New entry point

`run.R` preserves the first CSV observation (the old script dropped it even though `read.csv` had consumed the header), validates quarterly continuity and works on log prices. The final eight observations are held out chronologically. A pre-specified SARIMA, `auto.arima` fitted on training data only, and seasonal naive are evaluated on the same holdout using log-scale RMSE and MAE.

The selected structure from the academic work predates this added holdout. This is a retrospective benchmark, not a fully untouched prospective evaluation. Future work should use rolling-origin validation on fresh periods.

After evaluation, the fixed SARIMA is refitted on the complete series to forecast eight future quarters. Manual point predictions are compared with the package recursion; prediction intervals come from the `forecast` package. The manual seasonal function approximates uncertainty beyond a seasonal cycle by treating a forecast seasonal base as deterministic. Its bounds are educational and are not used as the main interval output. Log-scale point forecasts are not bias-adjusted arithmetic means on the original index scale.

The additional historical `IC_manual` function uses a constant residual standard deviation at every horizon and does not propagate forecast uncertainty. It is preserved only in `analysis/research_workflow.R`. The undefined `pred_manual` / `ic_manual` variables in that workflow were initialized during preparation.

## Validation status

The CSV was checked for schema, positive values, duplicate dates and consecutive quarters. Source functions and script dependencies were reviewed. The supplied report supports the original study, but **R is not installed in the preparation environment**, so the new R entry point has not been executed and no fresh forecast scores are claimed. `reports/historical_series.png` is a plot of the supplied observed data, not a model forecast. After running locally, retain `reports/current/sessionInfo.txt` with the resulting metrics.

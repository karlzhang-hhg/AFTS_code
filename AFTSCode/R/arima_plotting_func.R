library(fUnitRoots)
library(TSA)
library(forecast)
library(IRdisplay)
library(logger)

#' Install a package and direct output to a log file.
#'
#' Install package and direct output to a log file.
#'
#' @param lib A string of library name.
#' @param log_fpath A log file path to save installation output.
#' @return NULL
#' @export
#' @examples
#' install_redir_output("TSA", TSA_output.log)
install_redir_output <- function(lib, log_fpath) {
    # Open a file connection for writing
    sink(log_fpath, append = TRUE)

    # Install packages (output will be redirected to the file)
    install.packages(lib)

    # Close the file connection
    sink()
}

#' Plot time figure
#'
#' Plot time figure
#'
#' @param ts A ts object.
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @export
plot_time_fig <- function(ts, main = NULL, xlab = NULL, ylab = NULL) {
    par(bg = "white")
    plot(ts, type = "l", main = main, xlab = xlab, ylab = ylab)
    points(ts, pch = "*")
}

#' Calculate \phi_0 from AR model.
#'
#' @param ar_mod An AR model.
#' @param ord The order of AR model.
#' @param digits The number of digits to show float number.
#' @export
cal_phi_0 <- function(ar_mod, ord, digits = 6) {
    phi_0 <- as.numeric((1 - sum(ar_mod$coef[1:ord])) * ar_mod$coef["intercept"])
    print(format(phi_0, digits = digits))
    phi_0
}

#' Get \mu from ARIMA model.
#'
#' @param arima_mod An ARIMA model.
#' @param digits The number of digits to show float number.
#' @export
get_mu <- function(arima_mod, digits = 6) {
    mu <- as.numeric(arima_mod$coef["intercept"])
    print(format(mu, digits = digits))
    mu
}

#' Draw an ARIMA model and plot forecast results.
#' 
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param ylim ylim.
#' @param ts_fc_res An object already fit.
#' @import forecast
#' @export
draw_arima_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency, ts_fc_res,
    main=NULL, xlab=NULL, ylab=NULL, ylim = NULL
) {
    # Plot forecast
    par(mfrow = c(1, 1), bg = "white")
    if (is.null(npts)) {
        npts <- eotr
    }
    xmin <- time(da_ts)[eotr] - npts / frequency
    xmax <- time(da_ts)[eotr] + (max(h, length(da_ts) - eotr) + 1) / frequency
    cat(xmin, ";", xmax)
    tr_da_ts <- ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
    # # Label 1: Actual Observation Line
    plot(ts_fc_res, xlim = c(xmin, xmax), ylim = ylim, main = main, xlab = xlab, ylab = ylab)
    # Plot forecast mean (prepend the last observed data in the training dataset)
    dummy_1st_fmean_ts <- ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$mean)), frequency = frequency, start = end(tr_da_ts))
    # # Label -: NULL
    lines(dummy_1st_fmean_ts, col = "red")
    # # Label 2: Forecast Mean
    points(dummy_1st_fmean_ts, pch = 1)
    # Plot confidence interval (95%)
    dummy_1st_flower_ts <- ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$lower[, 2])), frequency = frequency, start = end(tr_da_ts))
    dummy_1st_fupper_ts <- ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$upper[, 2])), frequency = frequency, start = end(tr_da_ts))
    # # Label 3: Forecast 95% Lower Bound
    lines(dummy_1st_flower_ts, lty = 2)
    # # Label 4: Forecast 95% Upper Bound
    lines(dummy_1st_fupper_ts, lty = 2)
    # Plot original data
    orig_plot_ts <- ts(da_ts[(eotr - npts + 1):length(da_ts)], frequency = frequency, start = time(da_ts)[eotr] - (npts - 1) / frequency)
    # # Label -: NULL
    lines(orig_plot_ts)
    # # Label 5: Actual Observation Points
    points(orig_plot_ts, pch = 19)
    legend(
        "topleft", 
        legend = c(
            "Actual Obs Line", NULL, "Forecast Mean", "Forecast 95% Lower Bound",
            "Forecast 95% Upper Bound", NULL, "Actual Obs"
        ), 
        # col = c("black", "red", "blue"), 
        lty = c(1, NA, 2, 2, NA),
        pch = c(NA, 1, NA, NA, 19)
    )
}

#' Parameter 'seasonal' regeneration for wrapper functions.
#' 
#' @return Turn the `seasonal=NULL` to the default value for `arima(...)` and return a named list.
arima_seasonal_null_to_default <- function(...) {
    kwargs <- list(...)

    # Check if 'seasonal' is present
    if (('seasonal' %in% names(kwargs)) && is.null(kwargs$seasonal)) {
        # Modify the value of 'seasonal' if it exists and is NULL
        kwargs$seasonal <- list(order = c(0L, 0L, 0L), period = NA) # default value
    }
    
    return(kwargs)
}

#' Fit an ARIMA model and plot forecast results.
#'
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param xreg arima param (data.frame).
#' @param gof The lag for ts result diagnosis (36).
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param ylim ylim.
#' @param ts_fc_res An object already fit.
#' @import TSA forecast
#' @export
#' @examples
#' da = read.table("../AFTS_data/Ch02/q-jnj.txt", header = FALSE) #' needs to set `header = F`
#' jnj_er = da$V1
#' jnj_lg_er_ts = ts(jnj_lg_er, frequency = freq, start = c(1960, 1))
#' npts = 6
#' h = 8
#' eotr = length(jnj_lg_er_ts)-h
#' freq = 4
#' order = c(0,1,1) #' From the solution of 2-6
#' seasonal = list(order = order, period = freq)
#' fixed = NULL
#' include.mean = F
#' multi_seas_mod_jnj_ts = ts(jnj_lg_er_ts, frequency = freq, start = c(1960, 1))
#' multi_seas_mod_jnj_fc_res = plot_arima_forecast_fig(
#'     da_ts=multi_seas_mod_jnj_ts, eotr=eotr, h=h, npts=npts, frequency=freq, 
#'     main="Forecasts from ARIMA(0,1,1)(0,1,1)[4] with non-zero mean for\n(1-B)(1-B^4)ln(QtlyEarning) = (1-theta*B)(1-Theta*B^4)a_t of JnJ", 
#'     xlab="Year", ylab="QtlyEarning", ylim=c(2, 3.2),
#'     order=order, 
#'     seasonal=seasonal, 
#'     fixed=fixed, 
#'     method='ML', 
#'     include.mean=include.mean, 
#'     transform.pars=TRUE
#' )
#' summary(multi_seas_mod_jnj_fc_res)
#' 
#' multi_seas_mod_jnj_fc_tb = comb_forecast_res(
#'     multi_seas_mod_jnj_fc_res, 
#'     da_ts = multi_seas_mod_jnj_ts,
#'     eotr=eotr,
#'     freq=freq
#' )
#' multi_seas_mod_jnj_fc_df = as.data.frame(multi_seas_mod_jnj_fc_tb)
#' multi_seas_mod_jnj_fc_tb
plot_arima_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency, xreg=NULL, gof=36,
    main=NULL, xlab=NULL, ylab=NULL, ylim = NULL, ts_fc_res=NULL,
    ...
) {
    if (is.null(ts_fc_res)) {
        if (is.null(xreg)) {
            tr_xreg <- NULL
            fc_xreg <- NULL
        } else {
            stopifnot("xreg should be of type matrix and numeric"=(is.matrix(xreg) && is.numeric(xreg)))
            stopifnot("length(da_ts)!=dim(xreg)[1]"=(length(da_ts)==dim(xreg)[1]))
            tr_xreg <- xreg[1:eotr]
            fc_xreg <- xreg[(eotr+1):dim(xreg)[1]]
        }
        # arima model
        if (eotr > length(da_ts)) {
            err_msg = sprintf("The eotr (end-of-training) (%d) should not be larger than length(da_ts) (%d).", eotr, length(da_ts))
            log_error(err_msg)
            stop(err_msg)
        }
        tr_da_ts <- ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
        arima_kwargs <- arima_seasonal_null_to_default(...)
        arima_kwargs$y <- tr_da_ts
        if (!is.null(xreg)) {
            arima_kwargs$xreg <- tr_xreg
        }
        # Using TSA::arima won't be able to use the xreg in forecast.
        ts_fm <- do.call(forecast::Arima, arima_kwargs)
        print(ts_fm$nobs)
        # Forecast
        ts_fm$x <- tr_da_ts # https://stackoverflow.com/a/42464130/4307919
        if (is.null(xreg)) {
            ts_fc_res <- forecast(ts_fm, h = h)
        } else {
            ts_fc_res <- forecast(ts_fm, h = h, xreg = fc_xreg)
        }
        par(bg = 'white')
        tsdiag1(ts_fc_res$model, gof.lag = gof)
    }
    draw_arima_forecast_fig(da_ts, eotr, h, npts, frequency, ts_fc_res, main, xlab, ylab, ylim)
    ts_fc_res
}

#' Fit an ARIMA model using auto.arima and plot forecast results.
#'
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param xreg arima param (data.frame).
#' @param gof The lag for ts result diagnosis (36).
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param ylim ylim.
#' @param ts_fc_res An object already fit.
#' @param ... auto.arima params.
#' @import forecast logger
#' @export
#' @examples
#' da = read.table("../AFTS_sol/data/d-ibm3dxwkdays8008.txt", header = TRUE)
#' da[1:5,]
#' ew = da$ew * 100
#' ew_ts = ts(ew, frequency = 252, start = c(1980, 1, 2))
#' npts = 20
#' eotr = length(ew_ts)-npts
#' h = npts
#' freq = 252
#' xreg = as.matrix(da[, 8:11])
#' ew_fc_res = plot_auto_arima_forecast_fig(
#'     da_ts=ew_ts, eotr=eotr, h=h, npts=npts, frequency=freq,
#'     xreg=xreg,
#'     main="Forecasts from ARIMA(2,0,2)\nfor CRSP Equal-Weighted Index",
#'     xlab="Year", ylab="CRSP Equal-Weighted Index"
#'     d=0,
#'     D=0,
#'     max.p=2,
#'     max.q=2,
#'     max.P=1,
#'     max.Q=0,
#'     max.order=5,
#'     seasonal=TRUE,
#'     method="ML",
#'     allowmean=TRUE,
#'     stepwise=FALSE,
#'     parallel=TRUE,
#'     num.cores=12
#' )
#' summary(ew_fc_res)
#' 
#' ew_fc_tb = comb_forecast_res(ew_fc_res, ew_ts, eotr, freq)
#' ew_fc_tb
plot_auto_arima_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency, xreg=NULL, 
    gof=36, main=NULL, xlab=NULL, ylab=NULL, ylim = NULL, ts_fc_res = NULL,
    ...
) {
    if (is.null(ts_fc_res)) {
        if (is.null(xreg)) {
            tr_xreg <- NULL
            fc_xreg <- NULL
        } else {
            stopifnot("xreg should be of type matrix and numeric"=(is.matrix(xreg) && is.numeric(xreg)))
            stopifnot("length(da_ts)!=dim(xreg)[1]"=(length(da_ts)==dim(xreg)[1]))
            tr_xreg <- xreg[1:eotr]
            fc_xreg <- xreg[(eotr+1):dim(xreg)[1]]
        }
        # arima model
        if (eotr > length(da_ts)) {
            err_msg = sprintf("The eotr (end-of-training) (%d) should not be larger than length(da_ts) (%d).", eotr, length(da_ts))
            log_error(err_msg)
            stop(err_msg)
        }
        tr_da_ts <- ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
        if (is.null(xreg)) {
            # If I don't do this, I will get:
            # 
            # Error in eval(expr, p): object 'tr_xreg' not found
            # Traceback:

            # 1. plot_auto_arima_forecast_fig(da_ts = lg_gdpdef_ts, eotr = eotr, 
            # .     h = h, npts = npts, frequency = freq, main = NULL, xlab = "Year", 
            # .     ylab = "lg(GDP Deflator)", ylim = c(4.7, 4.9), d = NA, D = 0, 
            # .     max.p = 6, max.d = 3, max.q = 6, max.P = 0, max.Q = 0, max.order = 12, 
            # .     seasonal = TRUE, method = "ML", allowmean = TRUE, stepwise = FALSE, 
            # .     parallel = TRUE, num.cores = 4)
            # 2. forecast(ts_fm, h = h, xreg = fc_xreg)   # at line 26 of file <text>
            # 3. forecast.forecast_ARIMA(ts_fm, h = h, xreg = fc_xreg)
            # 4. predict(object, n.ahead = h)
            # 5. predict.Arima(object, n.ahead = h)
            # 6. eval.parent(xr)
            # 7. eval(expr, p)
            # 8. eval(expr, p)
            # 
            # This is because in the `predict.Arima`, we call `eval.parent(xr)`, where `xr <- object$call$xreg`.
            # Seems that if we assign NULL to a variable and pass it as xreg, it will not be exposed to the 
            # inner call environment and `predict.Arima` cannot find it. If I define a global variable `temp=NULL`
            # and pass it to the `auto.arima`, `eval.parent(xr)` works. So here, we have top explicitly pass `NULL`,
            # instead of some local variable that evaluates as a NULL.
            ts_fm <- auto.arima(tr_da_ts, ...) # forecast::auto.arima
        } else {
            ts_fm <- auto.arima(tr_da_ts, xreg = tr_xreg, ...) # forecast::auto.arima
        }
        print(ts_fm$nobs)
        # Forecast
        ts_fm$x <- tr_da_ts # https://stackoverflow.com/a/42464130/4307919
        if (is.null(xreg)) {
            ts_fc_res <- forecast(ts_fm, h = h)
        } else {
            ts_fc_res <- forecast(ts_fm, h = h, xreg = fc_xreg)
        }
        par(bg = 'white')
        tsdiag1(ts_fc_res$model, gof.lag = gof)
    }
    draw_arima_forecast_fig(da_ts, eotr, h, npts, frequency, ts_fc_res, main, xlab, ylab, ylim)
    ts_fc_res
}

#' Combine forecast results.
#'
#' @param forecast_obj A forecast object.
#' @param da_ts An original ts object.
#' @param eotr The end index of training data.
#' @param freq The frequency of the ts object.
#' @param stdout_len_lmt The length limit for printing out summary(forecast_obj) in stdout.
#' @export
comb_forecast_res <- function(forecast_obj, da_ts, eotr, freq, stdout_len_lmt = 5000) {
    out_len = sum(nchar(capture.output(summary(forecast_obj))))
    if (out_len <= stdout_len_lmt) {
        display(summary(forecast_obj))
    } else {
        display(paste("The \"summary(forecast_obj)\" call output length", out_len, "exceeds the output limit", stdout_len_lmt, "so that it is suppressed."))
    }
    fc_std_err <- (forecast_obj$upper[, 2] - forecast_obj$lower[, 2]) / 2 / qnorm(p = 0.975)
    display(forecast_obj$mean)
    display(fc_std_err)
    if (eotr + 1<=length(da_ts)) {
        actual_ts <- ts(da_ts[(eotr + 1):length(da_ts)], frequency = freq, start = time(da_ts)[eotr + 1])
        display(actual_ts)
        err_msg = sprintf(
            "The forecast$mean start time (%s) and actual values start time (%s) should be the same time", 
            paste(start(forecast_obj$mean), collapse = ","),
            paste(start(actual_ts), collapse = ",")
        )
        if (!all(start(forecast_obj$mean) == start(actual_ts))) {
            stop(err_msg)
        }
        multistep_ahead_forecast_tb <- cbind(forecast_obj$mean, fc_std_err, actual_ts)
        dimnames(multistep_ahead_forecast_tb)[[2]] <- c("Forecast", "Std. Error", "Actual")
    } else {
        multistep_ahead_forecast_tb <- cbind(forecast_obj$mean, fc_std_err)
        dimnames(multistep_ahead_forecast_tb)[[2]] <- c("Forecast", "Std. Error")
    }
    multistep_ahead_forecast_tb
}

#' Plot ACF and PACF for unit-roots testing.
#'
#' @param da_ts An original ts object.
#' @param freq The frequency of the ts object.
#' @param start The start of the diff(da_ts) object.
#' @param lag_max The lag.max param for acf().
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @export
plot_unit_root_figs <- function(da_ts, freq, start, lag_max, main, xlab, ylab) {
    par(bg = "white")
    diff_da_ts <- ts(diff(da_ts), frequency = freq, start = start)
    plot(da_ts, col = "black", type = "l", main = main, xlab = xlab, ylab = ylab)
    acf(da_ts, lag.max = 20)
    plot_time_fig(diff_da_ts, main = main, xlab = xlab, ylab = paste("diff(", ylab, ")", sep = ""))
    pacf(diff_da_ts, lag.max = 20)
}

#' Plot acf of a time-series.
#'
#' @param da A data series.
#' @param lag.max An acf param.
#' @param main The title of the figure.
#' @param w Plot width.
#' @param h Plot height.
#' @export
plot_acf <- function(da, lag.max = NULL, main=NULL, w=NULL, h=NULL, ...) {
    if (is.null(w) | is.null(h)) {
        par(bg = 'white')
    } else {
        par(bg = 'white', pin = c(w, h))
    }
    plot(acf(da, lag.max = lag.max, plot = F, ...), main = main)
}

#' Plot pacf and acf of a time-series.
#'
#' Don't pass in "`drop.lag.0` = T", b/c pacf doesn't accept this param.
#'
#' @param da A data series.
#' @param freq A frequency in taking seasonal difference.
#' @import TSA
#' @export
#' @examples 
#' plot_pacf_acf(diff(gdpdef), lag.max = lag.max, xlim = c(0, lag.max))
plot_pacf_acf <- function(da, freq=NA, ...) {
    par(mfrow = c(2, 1), bg = 'white')
    pacf(da, ...)
    TSA::acf(da, ...)
    pacf(diff(da), ...)
    TSA::acf(diff(da), ...)
    if (!is.na(freq)) {
        pacf(diff(da, freq), ...)
        TSA::acf(diff(da, freq), ...)
        pacf(diff(diff(da), freq), ...)
        TSA::acf(diff(diff(da), freq), ...)
    }
}

#' Perform and print eacf of a time-series.
#'
#' @param da A data series.
#' @param ar.max An eacf param.
#' @param ma.max An eacf param.
#' @import TSA IRdisplay
#' @export
perform_and_print_eacf <- function(da, ar.max, ma.max) {
    eacf_obj <- eacf(da, ar.max = ar.max, ma.max = ma.max)
    eacf_stats_tb = format(as.data.frame(eacf_obj$eacf), digits = 3)
    names(eacf_stats_tb) <- seq(from = 0, to = ma.max)
    display(eacf_stats_tb)
    display(eacf_obj$symbol)
    # pp67, asymptotic standard error of EACF
    display(2/sqrt(length(da)))
    c(eacf_obj, eacf_stats_tb)
}

#' Perform diagnosis on models (my version).
#'
#' @param object Model.
#' @param gof.lag The lag.max.
#' @param ... acf params.
#' @import TSA
#' @export
#' @examples
#' tsdiag1(ew_lrtn_m1, gof.lag = 36, `drop.lag.0` = FALSE)
#' tsdiag1(ew_lrtn_m1, gof.lag = 36, `drop.lag.0` = TRUE) 
tsdiag1 <- function(object, gof.lag = 30, ...) {
    # Specify the range of lags you're interested in
    lags <- 1:gof.lag
    p_values <- numeric(length(lags))

    # Compute the Ljung-Box test for each lag
    for (i in seq_along(lags)) {
        test_result <- Box.test(object$residuals, lag = lags[i], type = "Ljung-Box")
        p_values[i] <- test_result$p.value
    }

    # Plot
    par(mfrow = c(3, 1), bg = 'white')
    # Residual
    std_resi = object$residuals/sqrt(object$sigma2)
    plot(std_resi, type = 'n', main = "Standardized Residuals", xlab = "Time")
    segments(x0 = time(std_resi), y0 = 0, x1 = time(std_resi), y1 = std_resi)
    abline(h = 0, col = "blue", lty = 2, lwd = 1)
    # acf
    TSA::acf(as.numeric(object$residuals), lag.max = gof.lag, main = "ACF of Residuals", ...)
    # Ljung-Box
    plot(
        lags, p_values, type = "p", pch = 1, xlab = "Lag", ylab = "p value",
        main = "p values for Ljung-Box statistics", ylim = c(0, 1)
    )
    abline(h = 0.05, col = "blue", lty = 2, lwd = 1)  # Add a line for the common significance level of 0.05
}

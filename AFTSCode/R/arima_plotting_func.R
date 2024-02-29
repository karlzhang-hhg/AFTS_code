library(fUnitRoots)
library(TSA)
library(forecast)
library(IRdisplay)

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

#' Draw a ARIMA model and plot forecast results.
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
    par(bg = "white")
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
    lines(dummy_1st_fmean_ts)
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

#' Fit a ARIMA model and plot forecast results.
#'
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param order arima param.
#' @param seasonal arima param.
#' @param fixed arima param.
#' @param method arima param.
#' @param include.mean arima param.
#' @param transform.pars arima param.
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param ylim ylim.
#' @param ts_fc_res An object already fit.
#' @import forecast
#' @export
plot_arima_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency,
    order, seasonal, fixed, method, 
    include.mean, transform.pars,
    main=NULL, xlab=NULL, ylab=NULL, ylim = NULL,
    ts_fc_res=NULL
) {
    if (is.null(ts_fc_res)) {
        # arima model
        if (is.null(seasonal)) {seasonal = list(order = c(0L, 0L, 0L), period = NA)} # default value
        tr_da_ts <- ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
        if (is.null(transform.pars)) {
            ts_fm <- arima(tr_da_ts, order = order, fixed = fixed, 
                seasonal = seasonal, method = method, include.mean = include.mean
            )
        } else {
            ts_fm <- arima(tr_da_ts, order = order, fixed = fixed, 
                seasonal = seasonal, method = method, include.mean = include.mean, 
                transform.pars = transform.pars
            )
        }
        print(ts_fm$nobs)
        # Forecast
        ts_fm$x <- tr_da_ts # https://stackoverflow.com/a/42464130/4307919
        ts_fc_res <- forecast(ts_fm, h = h)
    }
    draw_arima_forecast_fig(da_ts, eotr, h, npts, frequency, ts_fc_res, main, xlab, ylab, ylim)
    ts_fc_res
}

#' Fit a ARIMA model using auto.arima and plot forecast results.
#'
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param xreg arima param (data.frame).
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @param ylim ylim.
#' @param ts_fc_res An object already fit.
#' @param ... auto.arima params.
#' @import forecast
#' @export
plot_auto_arima_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency,
    xreg=NULL, main=NULL, xlab=NULL, ylab=NULL, ylim = NULL, ts_fc_res = NULL,
    ...
) {
    if (is.null(ts_fc_res)) {
        if (is.null(xreg)) {
            tr_xreg <- NULL
            fc_xreg <- NULL
        } else {
            stopifnot("xreg should be of type matrix"=(is.matrix(xreg)))
            stopifnot("length(da_ts)!=dim(xreg)[1]"=(length(da_ts)==dim(xreg)[1]))
            tr_xreg <- xreg[1:eotr]
            fc_xreg <- xreg[(eotr+1):dim(xreg)[1]]
        }
        # arima model
        if (is.null(seasonal)) {seasonal = list(order = c(0L, 0L, 0L), period = NA)} # default value
        tr_da_ts <- ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
        ts_fm <- auto.arima(tr_da_ts, xreg = tr_xreg, ...)
        print(ts_fm$nobs)
        # Forecast
        ts_fm$x <- tr_da_ts # https://stackoverflow.com/a/42464130/4307919
        ts_fc_res <- forecast(ts_fm, h = h, xreg = fc_xreg)
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
#' @export
comb_forecast_res <- function(forecast_obj, da_ts, eotr, freq) {
    display(summary(forecast_obj))
    fc_std_err <- (forecast_obj$upper[, 2] - forecast_obj$lower[, 2]) / 2 / qnorm(p = 0.975)
    actual_ts <- ts(da_ts[(eotr + 1):length(da_ts)], frequency = freq, start = time(da_ts)[eotr + 1])
    display(forecast_obj$mean)
    display(fc_std_err)
    display(actual_ts)
    multistep_ahead_forecast_tb <- cbind(forecast_obj$mean, fc_std_err, actual_ts)
    dimnames(multistep_ahead_forecast_tb)[[2]] <- c("Forecast", "Std. Error", "Actual")
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
#' @param da An data series.
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
    plot(acf(da, lag.max = lag.max, plot = F, `drop.lag.0` = F, ...), main = main)
}

#' Plot pacf and acf of a time-series.
#'
#' @param da An data series.
#' @param lag.max An acf param.
#' @param main The title of the figure.
#' @param w Plot width.
#' @param h Plot height.
#' @export
plot_pacf_acf <- function(da, lag.max = NULL, main=NULL, w=NULL, h=NULL, ...) {
    if (is.null(w) | is.null(h)) {
        par(mfrow = c(2, 1), bg = 'white')
    } else {
        par(mfrow = c(2, 1), bg = 'white', pin = c(w, h))
    }
    plot(pacf(da, lag.max = lag.max, plot = F, `drop.lag.0` = F, ...), main = main)
    plot(acf(da, lag.max = lag.max, plot = F, `drop.lag.0` = F, ...), main = main)
}

#' Perform and print eacf of a time-series.
#'
#' @param da An data series.
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

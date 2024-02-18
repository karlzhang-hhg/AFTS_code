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
plot_time_fig <- function(ts, main=NULL, xlab=NULL, ylab=NULL) {
    par(bg = "white")
    plot(ts, type = 'l', main = main, xlab = xlab, ylab = ylab)
    points(ts, pch = '*')
}

#' Calculate \phi_0 from AR model.
#' 
#' @param ar_mod An AR model.
#' @param ord The order of AR model.
#' @param digits The number of digits to show float number.
#' @export 
cal_phi_0 <- function(ar_mod, ord, digits=6) {
    phi_0 = as.numeric((1-sum(ar_mod$coef[1:ord]))*ar_mod$coef['intercept'])
    print(format(phi_0, digits = digits))
    phi_0
}

#' Get \mu from ARIMA model.
#' 
#' @param arima_mod An ARIMA model.
#' @param digits The number of digits to show float number.
#' @export 
get_mu <- function(arima_mod, digits=6) {
    mu = as.numeric(arima_mod$coef['intercept'])
    print(format(mu, digits = digits))
    mu
}

#' Get \mu from ARIMA model.
#' 
#' @param da_ts A ts object.
#' @param eotr The end index of training data.
#' @param h The forecast horizon.
#' @param npts The number of points at the end of training data to plot.
#' @param frequency The frequency of the ts object.
#' @param order ARIMA param.
#' @param fixed ARIMA param.
#' @param method ARIMA param.
#' @param transform_pars ARIMA param.
#' @param main The title of the figure.
#' @param xlab xlab.
#' @param ylab ylab.
#' @import forecast
#' @export 
plot_forecast_fig <- function(
    da_ts, eotr, h, npts, frequency, 
    order, fixed, method, transform_pars,
    main, xlab, ylab, ylim = NULL
) {
    par(bg = "white")
    # arima model
    tr_da_ts = ts(da_ts[1:eotr], frequency = frequency, start = start(da_ts))
    if (is.null(transform_pars)) {
        ts_fm3 = arima(tr_da_ts, order = order, fixed = fixed, method = method)
    } else {
        ts_fm3 = arima(tr_da_ts, order = order, fixed = fixed, method = method, transform.pars = transform_pars)
    }
    print(ts_fm3$nobs)
    # Forecast
    ts_fm3$x <- tr_da_ts # https://stackoverflow.com/a/42464130/4307919
    ts_fc_res = forecast(ts_fm3, h = h)
    # Plot forecast
    if (is.null(npts)) {npts = eotr}
    xmin = time(da_ts)[eotr]-npts/frequency
    xmax = time(da_ts)[eotr]+(max(h, length(da_ts)-eotr)+1)/frequency
    cat(xmin, ";", xmax)
    
    plot(ts_fc_res, xlim = c(xmin, xmax), ylim = ylim, main = main, xlab = xlab, ylab = ylab)
    # Plot forecast mean
    dummy_1st_fmean_ts = ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$mean)), frequency = frequency, start = end(tr_da_ts))
    lines(dummy_1st_fmean_ts)
    points(dummy_1st_fmean_ts, pch = 1)
    # Plot confidence interval
    dummy_1st_flower_ts = ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$lower[,2])), frequency = frequency, start = end(tr_da_ts))
    dummy_1st_fupper_ts = ts(c(c(da_ts[eotr]), as.numeric(ts_fc_res$upper[,2])), frequency = frequency, start = end(tr_da_ts))
    lines(dummy_1st_flower_ts, lty=2)
    lines(dummy_1st_fupper_ts, lty=2)
    # Plot original data
    orig_plot_ts = ts(da_ts[(eotr-npts+1):length(da_ts)], frequency = frequency, start = time(da_ts)[eotr]-(npts-1)/frequency)
    lines(orig_plot_ts)
    points(orig_plot_ts, pch = 19)
    ts_fc_res
}

#' Get \mu from ARIMA model.
#' 
#' @param forecast_obj A forecast object.
#' @param da_ts An original ts object.
#' @param eotr The end index of training data.
#' @param freq The frequency of the ts object.
#' @export 
comb_forecast_res <- function(forecast_obj, da_ts, eotr, freq) {
    display(summary(forecast_obj))
    fc_std_err = (forecast_obj$upper[,2]-forecast_obj$lower[,2])/2/qnorm(p = 0.975)
    actual_ts = ts(da_ts[(eotr+1):length(da_ts)], frequency = freq, start = time(da_ts)[eotr+1])
    display(forecast_obj$mean); display(fc_std_err); display(actual_ts)
    multistep_ahead_forecast_tb = cbind(forecast_obj$mean, fc_std_err, actual_ts)
    dimnames(multistep_ahead_forecast_tb)[[2]] <- c("Forecast", "Std. Error", "Actual")
    multistep_ahead_forecast_tb
}

#' Get \mu from ARIMA model.
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
    diff_da_ts = ts(diff(da_ts), frequency = freq, start = start)
    plot(da_ts, col = 'black', type = 'l', main = main, xlab = xlab, ylab = ylab)
    acf(da_ts, lag.max = 20)
    plot_time_fig(diff_da_ts, main = main, xlab = xlab, ylab = paste("diff(", ylab, ")", sep = ""))
    pacf(diff_da_ts, lag.max = 20)
}

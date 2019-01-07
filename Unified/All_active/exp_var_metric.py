import numpy as np
import scipy.signal as signal


def calculate_spike_time_metrics(expt_trains, model_train, total_length, dt, sigma):
    """
    expt_trains is a list of nparrays, where each array contains the time indices of spike times (from all experimental trials)
    model_train is a single nparray which contains the time indices of spike times from the model
    total_length is the total number of time bins
    dt is the time step between each bin
    sigma is the length of the Gaussian convolution window (same units as dt)
    """
    trialtotrial_expvar_result = []
    for s in sigma:
        sigma_points = s / dt
        window_length = 10 * sigma_points
        gauss_func = signal.gaussian(window_length, sigma_points)
        binary_expt_trains = []
        for et in expt_trains:
            bt = np.zeros(total_length)
            bt[et] = 1
            binary_expt_trains.append(bt)
        convolved_expt_trains = [signal.fftconvolve(gauss_func, bt) for bt in binary_expt_trains]
        avg_convolved_expt_train = np.mean(convolved_expt_trains, axis=0)
        binary_model_train = np.zeros(total_length)
        if len(model_train) > 0:
            binary_model_train[model_train] = 1
        convolved_model_train = signal.fftconvolve(gauss_func, binary_model_train)

        f = np.array(convolved_expt_trains)
        trialtotrial_expvar_result.append(trial_expvar(f, convolved_model_train) / trial_expvar(f, avg_convolved_expt_train))

    return np.array(trialtotrial_expvar_result)


def trial_expvar(f, m):
    if f.shape[1] != len(m):
        print "Comparison trace has a different number of points than trial traces"
        return np.nan

    fbar = np.mean(f)
    mbar = np.mean(m)
    var_sum = np.sum((f - fbar) ** 2) + f.shape[0] * np.sum((m - mbar) ** 2)
    ev = (var_sum - (np.sum(((f - fbar) - (m - mbar)) ** 2))) / var_sum
    return ev

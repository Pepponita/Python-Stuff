import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import nicefunctions as fun

# loading the data
testrun = np.loadtxt(open("Testrun_511keV.Spe").readlines()[:-16], skiprows=12)
M1 = np.loadtxt(open("measurement1_time_offset.Spe").readlines()[:-16], skiprows=12)
M2 = np.loadtxt(open("measurement2_lifetime_in_aluminium.Spe").readlines()[:-16], skiprows=12)
M3 = np.loadtxt(open("measurement3_main.Spe").readlines()[:-16], skiprows=12)
timedata = np.loadtxt(open("timecalibration.Spe").readlines()[:-16], skiprows=12)

# normally plotting the data using fun.plot
fun.plot(M1, "Testrun", filename="Testrun")
fun.plot(M1, "Measurement 1", filename="Measurement_1_plot")
fun.plot(M2, "Measurement 2", filename="Measurement_2_plot")
fun.plot(M3, "Measurement 3", filename="Measurement_3_plot")
fun.plot(timedata, "Time calibration", filename="Time_calibration")

# Timecalibration part
# First we distinguish the different gaussians using fun.plotrange and plot them (will be commented out as this has no further use currently)

fun.plotrange(timedata,8860,8920,"Peak 8500","peak_8500")
"""
fun.plotrange(timedata,10160,10200,"Peak 10200","peak_10200")
fun.plotrange(timedata,11500,11550,"Peak 11500","peak_11500")
fun.plotrange(timedata,12790,12820,"Peak 12800","peak_12800")
fun.plotrange(timedata,14120,14150,"Peak 14100","peak_14100")
"""

# Then we use the gaussfit with information on range, amplitude, mu and sigma, we also get the mean values

timegaussfit_mean_1 = fun.custom_fit(timedata, 8884, 8900, [8892, 3, 20000], fun.gauss_function, "Peak 8500",
                                     "Timecalibrationgauss/TimeGauss_peak_8500.png", gauss_bool=True, xaxis="$\Delta$t ADC units")
timegaussfit_mean_2 = fun.custom_fit(timedata, 10170, 10190, [10182, 3, 2000], fun.gauss_function, "Peak 10200",
                                     "Timecalibrationgauss/TimeGauss_peak_10200", gauss_bool=True,xaxis="$\Delta$t ADC units" )
timegaussfit_mean_3 = fun.custom_fit(timedata, 11515, 11530, [11522, 3, 1500], fun.gauss_function, "Peak 11500",
                                     "Timecalibrationgauss/TimeGauss_peak_11500.png", gauss_bool=True,xaxis="$\Delta$t ADC units")
timegaussfit_mean_4 = fun.custom_fit(timedata, 12800, 12820, [12809, 3, 1400], fun.gauss_function, "Peak 12800",
                                     "Timecalibrationgauss/TimeGauss_peak_12800.png", gauss_bool=True,xaxis="$\Delta$t ADC units")
timegaussfit_mean_5 = fun.custom_fit(timedata, 14125, 14145, [14136, 3, 3500], fun.gauss_function, "Peak 14100",
                                     "Timecalibrationgauss/TimeGauss_peak_14100.png", gauss_bool=True,xaxis="$\Delta$t ADC units")
timecal_mean_lyst = [timegaussfit_mean_1, timegaussfit_mean_2, timegaussfit_mean_3, timegaussfit_mean_4,
                     timegaussfit_mean_5]

# Now we create the linear fit using known values for the delay
delay_lyst = [0, 4, 8, 12, 16]
timetranslation = fun.linear_fit(timecal_mean_lyst, delay_lyst,print_bool=True,xaxis="$\Delta$t ADC units",yaxis="$\Delta$t [ns]")

# fun.plot(M1,"Testrun",filename= "Testrun_with_timecal",timetrans=timetranslation)

# Measurement 1 part
# we can do a gaussian fit for the first measurement + residuals
fun.plot(M1, "Measurement 1", filename="_1_plot_adjusted_time", timetrans=timetranslation)
time_offset = fun.custom_fit(M1, 7500, 8500, [7800 * timetranslation, 20 * timetranslation, 80], fun.gauss_function,
                             "Measurement 1",
                             "Fits/Measurement_1_gaussfit.png", timetrans=timetranslation, linewidthh=2,
                             residual_bool=True, residual_file="Measurement_1", gauss_bool=True, print_bool=True,print_text="Time offset (ns) measurement 1: ")
print("Time offset (ns): " + str(time_offset))

# Measurement 2 part
#
fun.plot(M2, "Measurement 2", filename="_2_plot_adjusted_time", timetrans=timetranslation, linewidthh=1)
fun.custom_fit(M2, 7500, 9000, [8000 * timetranslation, 10 * timetranslation, 12, 50 * timetranslation],
               fun.gauss_exp_convolution, "Measurement 2", "Fits/Measurement_2_gauss_exp_convolutionfit.png",
               timetrans=timetranslation, linewidthh=2, residual_bool=True, residual_file="Measurement_2",print_bool=True, print_text="Lifetime (ns) measurement 2: ")


# measurement 3 part
fun.plot(M3, "Measurement 3", filename="_3_plot_adjusted_time", timetrans=timetranslation)
fun.custom_fit(M3,7500,10000,[25, 0.18, 200, 200, 20, 1],fun.final_gauss_exp_convolution,"Measurement 3","Fits/Measurement_3_final_gauss_exp_convolutionfit.png",timetrans=timetranslation,linewidthh=2,residual_bool=True,residual_file="Measurement_3",print_bool=True,print_text="Lifetime (ns) measurement 3: ",final_convolution_bool=True)
#fun.plot_test(M3, 0, 18000, fun.final_gauss_exp_convolution, [27.2, 0.18, 200, 200, 2, 0.1],
#              timetrans=timetranslation)

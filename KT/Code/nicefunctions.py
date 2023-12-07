import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.special as sp
import math

def plot(data,title,filename,xaxis="$\Delta$t[ns]",yaxis="Counts",linewidthh = 1, timetrans = 1):
    """
    This function plots a given set of count data to delta t
    params:
    data = dataset
    filename = name of file
    title = title of plot
    xaxis = name of xaxis
    yaxis = name of yaxis
    linewidth = width of line
    """
    x_range = np.linspace(0,len(data),len(data))*timetrans
    plt.plot(x_range, data, linewidth=linewidthh)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.title(title)
    plt.grid()
    filestring = "Plots/measurement_" + filename + ".png"
    plt.savefig(filestring)
    #plt.show()
    plt.clf()
    return None

def plotrange(data,start,stop,title,filename,xaxis="$\Delta$t ADC units",yaxis="Counts",linewidthh = 1):
    """
    This function plots a given set of count data to delta t in a given range
    params:
    data = dataset
    start = startpoint
    stop = endpoint
    filename = name of file
    title = title of plot
    xaxis = name of xaxis
    yaxis = name of yaxis
    linewidth = width of line
    """
    reduced_x_axis =  []
    reduced_data = []
    for i in range(len(data)):
        if start<=i<=stop:
            reduced_x_axis.append(i)
            reduced_data.append(data[i])
    plt.plot(reduced_x_axis, reduced_data, linewidth=linewidthh)
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.title(title)
    plt.grid()
    filestring = "Timecalibration/Timecal_" + filename + ".png"
    plt.savefig(filestring)
    # plt.show()
    plt.clf()
    return None

def gauss_function(x,mu,sigma,a):
    """
    Just a gauss function :3 with a amplitude
    """
    return (1/(sigma*np.sqrt(2*np.pi)))*a*np.exp(-((x-mu)**2)/(sigma**2))

def gauss_exp_convolution(x,mu,sigma,a,tau):
    """
    This function is the convolution of a gaussian and an exponential function
    """
    return ((a/2)*np.exp(((sigma**2)-2*tau*x+2*mu*tau)/(2*(tau**2)))*sp.erfc((tau*(mu-x)+(sigma**2))/(np.sqrt(2)*tau*sigma)))

def partial_final_gauss_exp_convolution(data,mu,sigma,a,tau):
    """
    For easier overview we seperate the last convolution into parts
    """
    res=[]
    for x in data:
        first = a/2
        second=np.exp(((sigma**2)-2*tau*x+2*mu*tau)/(2*(tau**2)))
        third=math.erfc((tau*(mu-x)+(sigma**2))/(np.sqrt(2)*tau*sigma))
        res.append(first*second*third)
    return np.array(res)

def final_gauss_exp_convolution(x,mu,sigma,a_1,a_2,tau_1,tau_2):
    """
    This function is the final convolution taking into account para and ortho positronium
    """
    one = partial_final_gauss_exp_convolution(x,mu,sigma,a_1,tau_1)
    two = partial_final_gauss_exp_convolution(x,mu,sigma,a_2,tau_2)
    return one+two
    #return ((1/(sigma*np.sqrt(2*np.pi)))*np.exp(-((x-mu)**2)/(2*(sigma**2)))*((a_1*np.exp(-x/tau_1))+(a_2*np.exp(-x/tau_2))))

def custom_fit(data,start,stop,initial_guess,function,title,filename,xaxis="$\Delta$ t [ns]",yaxis="Counts",linewidthh = 1,timetrans = 1, residual_bool=False, residual_file=None,gauss_bool = False,final_convolution_bool=False,print_bool=False,print_text=""):
    """
    This function plots a fit for a given function in a given range
    params:
    data = dataset
    start = startpoint of fit
    stop = endpoint of fit
    initial_guess = initial guess for the fit
    function = function which will be fitted
    filename = name and path of file
    title = title of plot
    xaxis = name of xaxis
    yaxis = name of yaxis
    linewidthh = width of line
    timetrans = time translation, gotten from the calibration
    residual_bool = boolean which determines if the residual will be plotted
    residual_file = name of the residual file
    gauss_bool = boolean which determines if the gauss parameters will be plotted
    final_convolution_bool = boolean which determines if the final convolution parameters will be plotted
    print_bool = boolean which determines if the fit parameters will be printed
    print_text = text which will be printed before the fit parameters

    Returns the mean of the fit

    Additionally, if residual = True, the residual will be calculated
    """
    #empty lists which will fill with data for the fit and plot parameters
    fit_x_axis = []
    fit_data = []
    reduced_x_axis=[]
    reduced_data=[]
    #filling the lists
    for i in range(len(data)):
        if start <= i <= stop: #if data is in gauss fit range add it
            fit_x_axis.append(i*timetrans)
            fit_data.append(data[i])
        if start-12<=i<=stop+12: #we widen the range for more plot visibility
            reduced_x_axis.append(i*timetrans)
            reduced_data.append(data[i])

    #Praying that the fit works
    #initial_guess = [initial_mu*timetrans, initial_sigma*timetrans,initial_amplitude]
    popt, pcov = opt.curve_fit(function, fit_x_axis, fit_data,initial_guess)

    #plotting the data and fit
    linear_delta_t_gauss = np.linspace(start,stop,len(fit_data))*timetrans
    fit_y_fit = function(linear_delta_t_gauss,*popt)
    plt.plot(reduced_x_axis, reduced_data, linewidth=1,label="Data")
    plt.plot(linear_delta_t_gauss, fit_y_fit, linestyle="--",color="red", linewidth = linewidthh,label="Fit",alpha=1)

    #if it is a gaussian we plot the mean and std
    if gauss_bool == True:
        plt.vlines(popt[0],ymin=0,ymax=gauss_function(popt[0],*popt),linestyle="--", color="g", label="Mean")
        plt.vlines((popt[0]+popt[1]), ymin=0, ymax=gauss_function(popt[0]+popt[1], *popt),linestyle="--", color="y", label="STD")

    #For the final measurement, the two exponentials and the data is plotted specially since it is a sum of two functions
    if final_convolution_bool == True:
        exp_gauss_1_fit = partial_final_gauss_exp_convolution(linear_delta_t_gauss,popt[0],popt[1],popt[2],popt[4])
        exp_gauss_2_fit = partial_final_gauss_exp_convolution(linear_delta_t_gauss,popt[0],popt[1],popt[3],popt[5])
        plt.plot(linear_delta_t_gauss,exp_gauss_1_fit,linestyle="--",color="orange",label="Exponential 1")
        plt.plot(linear_delta_t_gauss,exp_gauss_2_fit,linestyle="--",color="purple",label="Exponential 2")
        #Additionally the area under the data is highlighed
        #plt.fill_between(linear_delta_t_gauss,fit_y_fit,0,where=fit_y_fit>=0,facecolor="red",alpha=0.5,label="Area under curve")

    #showing additional info like title and axis
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.title(title)
    plt.legend()
    plt.grid()
    #filestring = "Timecalibrationgauss/TimeGauss_" + filename + ".png" old filestring used for the name of the file #old name, ignore
    plt.savefig(filename)
    # plt.show()
    plt.clf()

    #We also introduce the Residual and plot it
    if residual_bool == True:
        residual_lyst=[] #the empty residual lyst
        #we can reuse the restricted data and fitted data for this
        for i in range(len(fit_data)):
            if fit_data[i]==0: #to not have division by zero error we use this if statement
                temp_residual = 0
            else:
                temp_residual = (fit_data[i]-fit_y_fit[i])/(fit_data[i])
            residual_lyst.append(temp_residual)
        #now we can plot the residuals
        plt.plot(fit_x_axis,residual_lyst)
        plt.xlabel(xaxis)
        plt.ylabel("Residuals")
        plt.title("Residuals of the gaussian fit of " + title)
        plt.grid()
        residualstring = "Residual_"+residual_file
        filestring_residual = "Residuals/" + residualstring + ".png"
        plt.savefig(filestring_residual)
        #plt.show()
        plt.clf()

    #printing the fit parameters and their errors
    if print_bool == True:
        print("----------------------------------------")
        print(print_text + "\n")
        print("The values for the fitted parameters of the measurement are: ")
        print(popt)
        pcov_diag = np.sqrt(np.diag(pcov))
        print("The error on these parameters is given by: ")
        print(pcov_diag)
        print("----------------------------------------\n")

    return popt[0]

def linear_function(x,a,b):
    """
    Just a linear function :D
    """
    return a+b*x

def linear_fit(x_data,y_data,filename = "time_calibration",xaxis = "Channel",yaxis = "Delta t",plot_title = "channel vs time",print_bool = False):
    """
    This function performs a linear fit

    print_bool = boolean which determines if the fit parameters will be printed
    """
    #creating bset fit line and line data
    popt, pcov = opt.curve_fit(linear_function, x_data, y_data)
    linspace_x = np.linspace(x_data[0],x_data[-1],100)
    lin_y_fit = linear_function(linspace_x,*popt)
    #plotting
    plt.plot(linspace_x,lin_y_fit,linestyle="--",label="Linear fit")
    plt.scatter(x_data,y_data,color="r",label = "Mean")
    plt.xlabel(xaxis)
    plt.ylabel(yaxis)
    plt.title("Linear fit of "+plot_title)
    plt.legend()
    plt.grid()
    filestring = "Plots/Linear_fit_" + filename + ".png"
    #plt.show()
    plt.savefig(filestring)
    plt.clf()
    #additionally we return the translation of channel to time
    t_0=linear_function(x_data[0],*popt)
    t_1 = linear_function(x_data[-1],*popt)
    Deltat=t_1-t_0
    print("Delta t: "+str(Deltat))
    Deltax = x_data[-1]-x_data[0]
    print("Delta x: "+str(Deltax))
    TimeperChannel = Deltat/Deltax
    print("Time/Channel: "+str(TimeperChannel))

    if print_bool == True:
        print("----------------------------------------")
        print("Linear fit data:\n")
        print("The values for the fitted parameters of the measurement are: ")
        print(popt)
        pcov_diag = np.sqrt(np.diag(pcov))
        print("The error on these parameters is given by: ")
        print(pcov_diag)
        print("----------------------------------------\n")

    return TimeperChannel

def plot_test(data,start,stop,function,parameters,timetrans=1):
    """
    This function is used to test the fit functions
    """
    x_range = np.linspace(start*timetrans,stop*timetrans,len(data))
    y_range = function(x_range,*parameters)
    plt.plot(x_range,y_range,color="r",label="Fit")

    x_red = []
    dat_red = []
    for i in range(len(data)):
        if start<=i<=stop:
            x_red.append(i*timetrans)
            dat_red.append(data[i])
    plt.plot(x_red,dat_red,alpha=0.5,label="Data")
    plt.legend()
    plt.xlabel("time [ns]")
    plt.ylabel("Counts")
    plt.grid()
    plt.show()
    return None
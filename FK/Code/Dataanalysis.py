import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def get_data(data):
    """
    This function gets the data (returns 2 arrays, temperature and average resistance
    """
    temperature = data[:, 5]  # CHANGE TO 6 IF NOT TAFTER BUT TBASE
    resistance_pos = data[:, 2] / data[:, 1]
    resistance_neg = data[:, 4] / data[:, 3]
    average_resistance = (resistance_pos + resistance_neg) / 2
    return temperature, average_resistance

def datalimit(data,lowestTemp,highestTemp,number):
    """
    This function rescales the data according to the set boundaries lowest and highest temperature
    number: string, number of the measurement
    data: data to input
    """
    temperature = data[:, 5] #CHANGE TO 6 IF NOT TAFTER BUT TBASE
    resistance_pos = data[:, 2] / data[:, 1]
    resistance_neg = data[:, 4] / data[:, 3]
    average_resistance = (resistance_pos + resistance_neg) / 2
    new_temperature = []
    new_resistance = []

    for i in range(len(temperature)):
        temp = temperature[i]
        if lowestTemp<temp<highestTemp:
            new_temperature.append(temp)
            new_resistance.append(average_resistance[i])

    filestring = "Improved plots/rescaled_measurement_" + number + ".png"  # name of the file and where it goes
    titlestring = "Rescaled measurement " + number
    plt.scatter(new_temperature, new_resistance, s=3)
    plt.title(titlestring)
    plt.xlabel("Temperature [K]")
    plt.ylabel("Resistance [$\Omega$]")
    plt.grid()
    plt.savefig(filestring,bbox_inches="tight")
    plt.clf()
    return (new_temperature, new_resistance)

def plot(data, number,svalue = 0.8):
    """
    This function plots a measurement using the average of the positive and negative values
    for the resistance given by the measurement
    """
    #take data and make calculate average
    temperature = data[:, 5]#CHANGE TO 6 IF NOT TAFTER BUT TBASE
    resistance_pos = data[:, 2] / data[:, 1]
    resistance_neg = data[:, 4] / data[:, 3]
    average_resistance = (resistance_pos + resistance_neg) / 2
    #plot the data
    filestring = "Plots/measurement_" + number + ".png" #name of the file and where it goes
    titlestring = "Measurement " + number
    plt.scatter(temperature, average_resistance, svalue)
    plt.xlabel("Temperature [K]")
    plt.ylabel("Resistance [$\Omega$]")
    plt.title(titlestring)
    plt.grid()
    plt.savefig(filestring,bbox_inches="tight")
    #plt.show()
    plt.clf()
    return None

def linear_function(x,a,b):
    return a + b*x

def linearfit(data,lowestTemp,highestTemp,number,s=5,space=2,middlepoint = True):
    """
    This function linearly fits the transistion of the phaseshift and plots it
    """
    temperature = data[:, 5]#CHANGE TO 6 IF NOT TAFTER BUT TBASE
    resistance_pos = data[:, 2] / data[:, 1]
    resistance_neg = data[:, 4] / data[:, 3]
    average_resistance = (resistance_pos + resistance_neg) / 2
    new_temperature = []
    new_resistance = []
    fit_temperature = []
    fit_resistance = []
    #For loop for the resistance in the selected temperature range with added padding
    for i in range(len(temperature)):
        temp = temperature[i]
        if lowestTemp-space<temp<highestTemp+space:
            new_temperature.append(temp)
            new_resistance.append(average_resistance[i])
    #For loop for the same but now for the best fit line without padding
    for i in range(len(temperature)):
        temp = temperature[i]
        if lowestTemp<temp<highestTemp:
            fit_temperature.append(temp)
            fit_resistance.append(average_resistance[i])
    #best fit line
    popt, pcov = opt.curve_fit(linear_function,fit_temperature,fit_resistance)

    #Random stuff IGNORE
    #temptemp = linear_function(293,*popt)
    #print("Resistance at 293 K" + str(temptemp))
    #print(str(popt)+str(number))

    #name and scatter plot
    filestring = "Linear fits/Linear_fit_measurement_" + number + ".png"  # name of the file and where it goes
    titlestring = "Linear fit measurement " + number
    plt.scatter(new_temperature, new_resistance, s)
    #plotting the best fit line
    temp_range = np.linspace(lowestTemp,highestTemp,100)
    linear_resistance = linear_function(temp_range,*popt)
    plt.plot(temp_range,linear_resistance,color="red")
    #plotting the middle point
    middle_temp = (lowestTemp+highestTemp)/2
    print(middle_temp,number)
    if middlepoint == True:
        plt.plot(middle_temp,linear_function(middle_temp,*popt),color = "g", marker ="o")
        plt.errorbar(middle_temp,linear_function(middle_temp,*popt),0,0.7,color="black",capsize=5)
    #The rest of the plot
    plt.title(titlestring)
    plt.xlabel("Temperature [K]")
    plt.ylabel("Resistance [$\Omega$]")
    plt.grid()
    plt.savefig(filestring, bbox_inches="tight")
    #plt.show()
    plt.clf()



#data, no magnetic field
M1 = np.loadtxt("measurement_1_5mA",skiprows=7)
M2 = np.loadtxt("measurement_2_10mA",skiprows=7)
M3 = np.loadtxt("measurement_3_20mA",skiprows=7)
M4 = np.loadtxt("measurement_4_25mA.txt",skiprows=7) #corrupted file, data may be missing
M5 = np.loadtxt("measurement_5_30mA",skiprows=7)
#data, magnetic field
M6 = np.loadtxt("measurement_6_30mA_magneic_field_6_90",skiprows=7)
M7 = np.loadtxt("measurement_7_30_mgfield_6_90",skiprows=7)
M8 = np.loadtxt("measurement_8_with_mgfield_6_90",skiprows=7)
M9 = np.loadtxt("measurement_9_30mA_magneic_field_6_90",skiprows=7)
M10 = np.loadtxt("measurement_10_30mA_magneic_field_4_03",skiprows=7)
M11 = np.loadtxt("measurement_11_30mA_magneic_field_4_03",skiprows=7)
M12 = np.loadtxt("measurement_12_30mA_magneic_field_4_03",skiprows=7)

data_array = [M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12]
#overnight and std measurement
M_std = np.loadtxt("standard dev samle at equilibrium",skiprows=7)
M_overnight = np.loadtxt("overnight measrement 30mA",skiprows=7)


#Plotting of the data

"""
#We plot the file using a for loop and the function above
for i,data in enumerate (data_array,1):
    numberstr = str(i)
    plot(data,numberstr,5)

plot(M_std,"STD",3)
plot(M_overnight,"Overnight")
"""
#Limiting the data as needed


"""
#we now adjust the data to a better area for representation, this has to be done seperately for each value
datalimit(M2,75,90,"2")
datalimit(M3,130,140,"3")
datalimit(M4,98,102,"4")
datalimit(M5,98,102,"5")
#the magnetic field ones
datalimit(M9,103,107,"9")
datalimit(M12,80,85,"12")
datalimit(M_overnight,80,120,"overnight")
"""


#The Linear fit



linearfit(M1,94,96.5,"1",space=3)
linearfit(M3,94,97,"3",space=3)
linearfit(M4,93,93.7,"4")
linearfit(M5,93,93.8,"5")

linearfit(M_overnight,96,300,"overnight",middlepoint=False)

#STD measurement
std_raw_temp, std_raw_res = get_data(M_std)
std_temp = np.std(std_raw_temp)
std_res = np.std(std_raw_res)
print("STD temp: ", std_temp," STD res: ", std_res)

#!/usr/bin/env python3

import math
import os
import sys
import numpy as np
from scipy.interpolate import interp1d
# import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
# import palettable.scientific.sequential as colors
import palettable.cmocean.sequential as colors2
#from sklearn.metrics import r2_score

def read_logfile(logfile,do_eV=False):
    """
    Reads a Gaussian .log file and extracts the f and either
    E(eV) or E(lambda) values of each excitation.
    Note: to apply lambda shift in nm by make_spectrum(), 
          eV should be set to False
    """
    parse = False
    with open(logfile) as inp:
        lmbda = np.array([])
        eV = np.array([])
        f = np.array([])
        for line in inp:
            if "Excitation energies and oscillator strengths:" in line:
                parse = True
            elif "******************************" in line:
                parse = False
            if parse == True and "Excited State" in line:
                line.strip()
                dat = line.split()
                lmbda = np.append(lmbda,[float(dat[6])])
                eV = np.append(eV,[float(dat[4])])
                f = np.append(f,[float(dat[8][2:])])
    if do_eV == False:
        return lmbda, f
    else:
        return eV, f

def read_outfile(outfile,do_eV=False):
    """
    Reads an ORCA .out file and extracts the f and either
    E(eV) or E(lambda) values of each excitation.
    Note: to apply lambda shift in nm by make_spectrum(), 
          eV should be set to False
    """
    parse = False
    with open(outfile) as inp:
        lmbda = np.array([])
        eV = np.array([])
        f = np.array([])
        for line in inp:
            if "SPECTRUM VIA TRANSITION ELECTRIC DIPOLE" in line:
                parse = True
            elif "SPECTRUM VIA TRANSITION VELOCITY DIPOLE" in line: # TDDFT
                parse = False
            elif "CD SPECTRUM" in line: # CCSD
                parse = False
            if parse == True:
                line = line.strip()
                dat = line.split()
                if len(dat) >= 4:
                    try:
                        lmbda = np.append(lmbda,[float(dat[2])])
                        eV = np.append(eV,[1239.84193/float(dat[2])])
                        f = np.append(f,[float(dat[3])])
                    except ValueError:
                        pass
    # print(lmbda)
    # print(f)
    if do_eV == False:
        return lmbda, f
    else:
        return eV, f

def read_expt(csvfile,ref_pad=50):
    """
    Reads a .csv file and extracts the I vs E(lambda) 
    data from its two columns.
    The data is the transformed as follows:
        - I is interpolated to integer lambda values
          between min and max lambda (rounded to integer)
        - min(I) is shifted to zero
        - area under the spectrum is normalized to unity
        - normalization is done in the integer lambda range
    The integer lambda range should contain more points 
    than the input lambda range to avoid loss of detail.
    """
    E = np.array([])
    I = np.array([])
    with open(csvfile) as inp:
        for line in inp:
            if ';' in line and ',' in line:
                line = line.replace(',','.')
                line = line.replace(';',',')
            dat = line.split(',')
            E = np.append(E,float(dat[0]))
            I = np.append(I,float(dat[1])) 
    E_min = int(min(E)) + 1
    E_max = int(max(E))
    lambdarange = np.arange(E_min,E_max+1)
    f = interp1d(E, I, kind='cubic') # create an interpolation function
    I_ = f(lambdarange) # interpolate over the integer lambda values
    I_ = I_ - min(I_) # shift min(I_) to zero
    I_ = I_/np.trapz(I_,lambdarange) # normalize
    
    if ref_pad != 0: # add trailing zeroes to the data
        max_E = int(max(lambdarange))
        E_fill = np.arange(max_E+1,max_E+ref_pad+1)
        I_fill = np.zeros(len(E_fill))
        lambdarange = np.append(lambdarange,E_fill)
        I_ = np.append(I_,I_fill)

    return lambdarange, I_ # returning lambdarange should not be necessary

def save_expt(csvfile,ref_pad=50):
    lambdarange, I_ = read_expt(csvfile,ref_pad)
    outname = csvfile[:-4] + "_integerx.csv"
    if not os.path.exists(outname):
        k = open(outname, "w")
        for i in range(len(lambdarange)):
            newline = str(lambdarange[i]) + ',' + str(I_[i]) + '\n'
            k.write(newline)
        k.close()
    else:
        print(f'{outname} already exists, skipping reference csv transformation.')

def make_spectrum(energies,f,lambdarange,SIGMA,lambda_shift=1,int_range=[]):
    
    E = np.array([])
    try: # convert lambdarange to eV range
        eV_range = 1239.84193/lambdarange 
    except TypeError:
        eV_range = 1239.84193/np.array(lambdarange)
    for i in range(len(energies)): # apply the lambda-shift
        E = np.append(E,energies[i]*lambda_shift)
        # E[i] += lambda_shift
    E = 1239.84193/E # convert E(nm) to E(eV)
    
    Y = [0 for i in range(len(eV_range))]
    for i in range(len(Y)): # make spectrum in eV
        for n in range(len(E)):
            inexp = -0.5*((eV_range[i]-E[n])/SIGMA)**2 # all in eV
            eps = f[n]/(SIGMA*math.sqrt(2*math.pi))*math.exp(inexp)
            Y[i] += eps
    # Jacobi
    # hc = 1.986e-25
    # for i in range(len(eV_range)):
    #     Y[i] = (Y[i]*((eV_range[i])**2))/hc
    # X = np.flip(1239.84193/eV_range) # convert E(eV) to E(nm) and sort ascending
    # Y = np.flip(Y)
    X = 1239.84193/eV_range
    # fit to integer lambdas
    f = interp1d(X, Y, kind='cubic')
    I = f(lambdarange)
    # normalize
    if len(int_range) == 0:
        I = Y/np.trapz(I,lambdarange)
    else:
        splice_from = min(int_range)-min(lambdarange)
        splice_to = max(int_range)-max(lambdarange)
        if splice_to < 0: # the case where lambdarange is wider than int_range
            I = Y/np.trapz(I[splice_from:splice_to],int_range)
        else:
            int_lambdarange = lambdarange[splice_from:len(lambdarange)]
            # print(int_range[:len(int_lambdarange)])
            # print(int_lambdarange)
            I = Y/np.trapz(I[splice_from:len(lambdarange)],int_range[:len(int_lambdarange)])
    return lambdarange, I


def mae(predictions, targets):
    # Retrieving number of samples in dataset
    samples_num = len(predictions)
    # Summing absolute differences between predicted and expected values
    accumulated_error = 0.0
    for prediction, target in zip(predictions, targets):
        accumulated_error += np.abs(prediction - target)    
    # Calculating mean
    mae_error = (1.0 / samples_num) * accumulated_error    
    return mae_error

def mse(predictions, targets):
    # Retrieving number of samples in dataset
    samples_num = len(predictions)   
    # Summing square differences between predicted and expected values
    accumulated_error = 0.0
    for prediction, target in zip(predictions, targets):
        accumulated_error += (prediction - target)**2        
    # Calculating mean and dividing by 2
    mse_error = (1.0 / (2*samples_num)) * accumulated_error    
    return mse_error

def lg_mse(predictions, targets):
    return np.log(mse(predictions, targets))

def r_square(predictions, targets):
    """
    Calculates the square of the Pearson correlation coefficient.
    Correlations are useful to find patterns and relationships in 
    data but mostly useless to evaluate predictions.
    Return values may range from 0 to 1.
    Graphically, this can be understood as 
    “how close is the data to the line of best fit?”

    Parameters
    ----------
    predictions : numpy array
        Calculated values.
    targets : numpy array
        Reference values.

    Returns
    -------
    float
        The square of the Pearson correlation coefficient.

    """
    r_matrix = np.corrcoef(predictions, targets)
    r = r_matrix[0,1]
    return r**2

def R2(predictions, targets):
    """
    Calculates the coefficient of determination which 
    captures how well predictions match observations, 
    or how much of the variation in observed data is 
    explained by the predictions.
    Return values may range from -inf to 1. 
    Closer the data to the 1:1 line, 
    higher the coefficient of determination.
    Negative value means that using the mean of the 
    "targets" data yields a better prediction than using
    the "predictions" data.
    Also available as sklearn.metrics.r2_score()

    Parameters
    ----------
    predictions : numpy array
        Calculated values.
    targets : numpy array
        Reference values.

    Returns
    -------
    float
        The coefficient of determination.

    """
    #return r2_score(targets, predictions)
    #plot_predictions(predictions, targets)
    #plot_predictions_2(predictions, targets)
    #sys.exit(3)
    SSres = sum(map(lambda x: (x[0]-x[1])**2, zip(targets, predictions)))
    SStot = sum([(x-np.mean(targets))**2 for x in targets])
    return 1-(SSres/SStot)
    
def rmsle(predictions, targets):
    return np.sqrt(np.square(np.log(predictions + 1) - np.log(targets + 1)).mean())

def rmsd(predictions, targets):
    return np.sqrt(mse(predictions, targets))

def calc_error(X_ref,Y_ref,X,Y,errorfunc=mse):
    """
    Takes reference and predicted values 
    and calculates the error of the prediction.

    Parameters
    ----------
    X_ref : numpy array
        Reference X values.
    Y_ref : numpy array
        Reference Y values.
    X : numpy array
        Predicted X values.
    Y : numpy array
        Predicted Y values.
    errorfunc : function, optional
        Select the error metric. The default is mse.
        options: mae, mse, r_square, rmsle, lg_mse

    Returns
    -------
    float
        Returns the error value for a prediction.

    """
    T = {}
    P = {}
    for i in range(len(X_ref)):
        T[X_ref[i]] = Y_ref[i]
    for i in range(len(X)):
        P[X[i]] = Y[i]
    X_eval = np.intersect1d(X_ref,X)
    targets = np.array([T[i] for i in X_eval])
    predictions = np.array([P[i] for i in X_eval])
    return errorfunc(predictions, targets)

def get_errors_hmap(logfile,csvfile,lambdarange=[],ref_pad=50,errorfunc=mse,workdir='.',do_log=True):
    print("calculating errors...")
    molname = os.path.basename(logfile)[:-4]
    heatmap_name = molname+'_heatmap.svg'
    heatmap_name = os.path.join(workdir,heatmap_name)
    # print(f'workdir inside get_errors_hmap(): {workdir}')
    # print(f'heatmap name: {heatmap_name}')
    # print(f'error function passed as: {errorfunc}')
    # print(f'error function name: {errorfunc.__name__}')
    # sys.exit(3)
    errfn_name = errorfunc.__name__
    if len(lambdarange) == 0:
        lambdarange = np.arange(200,700)

    # SHIFTS = np.arange(70,131,0.5)/100 # 122 elements
    # SIGMAS = np.arange(10,45,0.4)/100 # 88 elements
    # SHIFTS = np.arange(70,131,0.75)/100 # 82 elements
    # SIGMAS = np.arange(10,45,0.4)/100 # 70 elements
    # SHIFTS = np.arange(70,131,0.65)/100 # 94 elements
    # SIGMAS = np.arange(10,45,0.4)/100 # 70 elements
    SHIFTS = np.linspace(0.7,1.3,90) # 90 elements
    SIGMAS = np.linspace(0.10,0.45,70) # 70 elements

    if logfile[-4:] == '.log': # Gaussian outputs
        E, f = read_logfile(logfile)
    elif logfile[-4:] == '.out': # ORCA outputs
        E, f = read_outfile(logfile)
    else:
        print(f"{logfile} is not recognized as a valid input for heatmap generation!")
        # sys.exit(2)
        return 2
    X_ref, Y_ref = read_expt(csvfile,ref_pad)

    int_range = np.arange(min(X_ref),max(X_ref)+1)
    lowest_err = 2e10
    highest_r2 = -2e10
    hmap_errors = np.array([], dtype=np.float64).reshape(0,len(SHIFTS))
    for n in SIGMAS:
        hmap_row = np.array([])
        for m in SHIFTS:
            X, Y = make_spectrum(E,f,lambdarange,n,m,int_range)
            current_err = calc_error(X_ref,Y_ref,X,Y,errorfunc)
            hmap_row = np.append(hmap_row,current_err)
            if current_err < lowest_err and errfn_name != "r_square" and errfn_name != "R2":
                lowest_sigma = n
                lowest_shift = m
                lowest_err = current_err
            elif current_err > highest_r2 and (errfn_name == "r_square" or errfn_name == "R2"):
                lowest_sigma = n
                lowest_shift = m
                highest_r2 = current_err
            #if (m > 0.95 and m < 1.05) and (n > 0.14 and n < 0.18):
            #    calc_error_debug(X_ref,Y_ref,X,Y)
            #    sys.exit(3) 
        if do_log == True and errfn_name != "r_square" and errfn_name != "R2":
            hmap_row = np.log(hmap_row) # use logarithmic scale
        hmap_errors = np.vstack([hmap_errors,hmap_row])

  
    fig, ax = plt.subplots()
    # im = ax.imshow(hmap_errors)
    # plt.imshow(hmap_errors, cmap=colors.Devon_20.mpl_colormap) # cmap="pink" gives good contrast but ugly color
    # plt.imshow(hmap_errors, cmap=colors.Tokyo_20.mpl_colormap)
    plt.imshow(hmap_errors, cmap=colors2.Haline_12.mpl_colormap)
    # We want to show all ticks...
    ax.set_xticks(np.arange(len(SHIFTS)))
    ax.set_yticks(np.arange(len(SIGMAS)))
    # ... and label them with the respective list entries
    # ax.set_xticklabels(SHIFTS)
    # ax.set_yticklabels(SIGMAS)
    ax.set_xticklabels(tic.FormatStrFormatter('%.2f').format_ticks(SHIFTS))
    ax.set_yticklabels(tic.FormatStrFormatter('%.2f').format_ticks(SIGMAS))
    # ax.contour(SHIFTS,SIGMAS, hmap_errors) #, levels=[40,80,100])
    ax.contour(np.arange(len(SHIFTS)),np.arange(len(SIGMAS)), hmap_errors, levels=15, linewidths=1.0, linestyles='dashed', colors = 'black')
    plt.xlabel('wavelength multiplication factor')
    plt.ylabel(r'$\Delta$ (eV)')
    n = 13  # Keeps every 10th label
    [l.set_visible(False) for (i,l) in enumerate(ax.xaxis.get_ticklabels()) if i % n != 0]
    [l.set_visible(False) for (i,l) in enumerate(ax.yaxis.get_ticklabels()) if i % n != 0]
    # plt.setp(ax.get_xticklabels()[::2], visible=False)
    # plt.setp(ax.get_yticklabels()[::2], visible=False)
    plt.colorbar().ax.tick_params(axis='y', direction='in')
    fig.tight_layout()
    plt.savefig(heatmap_name, bbox_inches='tight', pad_inches=0)
    # plt.clf()
    plt.close()
    # plt.show()
    # if errorfunc in [mse, mae, rmsle, rmsd]:
        # return lowest_sigma, lowest_shift, lowest_err
    # elif errorfunc == r_square:
        # return lowest_sigma, lowest_shift, highest_r2
    # else:
        # print("Please provide a valid error function!")
        # return 0,0,0
    if errfn_name == "r_square" or errfn_name == "R2":
        return lowest_sigma, lowest_shift, highest_r2
    else:
        return lowest_sigma, lowest_shift, lowest_err

def print_result(ref_csv,logfile,lambdarange,SIGMA,lambda_shift,ref_pad=50,dumpfile='foo.txt',writeopt=True,errorfunc=mse,workdir='.',logscale_hmap=True):
    if logfile[-4:] == '.log': # Gaussian outputs
        x,y = read_logfile(logfile)
    elif logfile[-4:] == '.out': # ORCA outputs
        x,y = read_outfile(logfile)
    else:
        print(f"{logfile} is not recognized as a valid input!")
        sys.exit(2)
    eX, eY = read_expt(ref_csv,ref_pad)
    max_eX = max(eX)
    int_range = np.arange(min(eX),max_eX+1)
    X,Y = make_spectrum(x,y,lambdarange,SIGMA,lambda_shift,int_range)
    e1 = calc_error(eX, eY, X, Y, errorfunc)
    print(f"{errorfunc.__name__} error calculated using sigma = {SIGMA} and no lambda_shift: {e1:.4e}")
    lowest_sigma, lowest_shift, lowest_err = get_errors_hmap(logfile,ref_csv,lambdarange,ref_pad,errorfunc,workdir,logscale_hmap)
    csvname = os.path.basename(logfile[:-4]) + "_opt_" + errorfunc.__name__ + ".csv"
    if type(dumpfile) == str and dumpfile != 'foo.txt':
        dumpfile = os.path.join(workdir,dumpfile)
        k = open(dumpfile, "a+")
        # dumpline = logfile + ',' + errorfunc.__name__ + ',' + str(lowest_sigma) + ',' + str(lowest_shift) + ',' + str(lowest_err) + '\n'
        dumpline = csvname + ',' + errorfunc.__name__ + ',' + str(lowest_err) + ',' + str(lowest_sigma) + ',' + str(lowest_shift) + '\n'
        k.write(dumpline)
        k.close()
    X_opt,Y_opt = make_spectrum(x,y,lambdarange,lowest_sigma,lowest_shift,int_range)
    if writeopt == True:
        # csvname = logfile[:-4] + "_opt_" + errorfunc.__name__ + ".csv"
        csvname = os.path.join(workdir,csvname)
        k = open(csvname, "w")
        for i in range(len(X_opt)):
            optline = str(X_opt[i]) + ',' + str(Y_opt[i]) + '\n'
            k.write(optline)
        k.close()
    print("optimized half-width: "+str(round(lowest_sigma,4))+" eV")
    print("optimized wavelength shift: "+str(round(lowest_shift,4))+" nm")
    print(f"{errorfunc.__name__} = {lowest_err:.4e}")
    # plt.plot(X,Y) # use these lines plot the spectrum with default parameters and experiment
    # plt.plot(X_opt,Y_opt,color='black')
    # plt.plot(eX,eY,linestyle='dashed',color='black')
    # plt.xlim(300,600)
    # plt.ylim(0,0.03)
    # plt.show()
    # plt.clf()
    # plt.close()

def plot_predictions(y_pred, y_test):
    x = range(len(y_pred))
    plt.scatter(x, y_pred, color='b')
    plt.scatter(x, y_test, color='r')
    plt.savefig("prediction_test.png", bbox_inches='tight', pad_inches=0)
    plt.clf()

def plot_predictions_2(y_pred, y_test):
    plt.scatter(y_pred, y_test)
    plt.savefig("prediction_test_xy.png", bbox_inches='tight', pad_inches=0)
    plt.clf()

def calc_error_debug(X_ref,Y_ref,X,Y,errorfunc=mse):
    T = {}
    P = {}
    for i in range(len(X_ref)):
        T[X_ref[i]] = Y_ref[i]
    for i in range(len(X)):
        P[X[i]] = Y[i]
    X_eval = np.intersect1d(X_ref,X)
    targets = np.array([T[i] for i in X_eval])
    predictions = np.array([P[i] for i in X_eval])
    print(predictions)
    print(targets)
    plot_predictions(predictions, targets)
    plot_predictions_2(predictions, targets)
    #return errorfunc(predictions, targets)


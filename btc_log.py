#!/usr/bin/env python
# coding: utf-8




#get_ipython().run_line_magic('pylab', 'notebook')
import scipy.stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import exists
from matplotlib.lines import Line2D
# https://www.canva.com/colors/color-wheel/
colors = ["green","blue"]







def mean_confidence_interval(pdf, confidence=0.95):

    num_sims,num_points = np.shape(pdf)
    mu = np.zeros(num_points)
    low = np.zeros(num_points)
    high = np.zeros(num_points)
    for i in range(num_points):
        data = pdf[:,i]
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
        mu[i] = m
        low[i] = m - h
        high[i] = m + h
        
    return mu,low,high

def mean_min_max(pdf, confidence=0.95):

    num_sims,num_points = np.shape(pdf)
    mu = np.zeros(num_points)
    low = np.zeros(num_points)
    high = np.zeros(num_points)
    for i in range(num_points):
        data = pdf[:,i]
        a = 1.0 * np.array(data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        mu[i] = m
        low[i] = min(data)
        high[i] = max(data)
        
    return mu,low,high


def weighted_avg_and_var(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = numpy.average(values, weights=weights)
    # Fast and numerically precise:
    variance = numpy.average((values-average)**2, weights=weights)
    return (average, np.sqrt(variance))

def extract_subvec(bx, pdf, idx):

    x_hat = bx[idx:idx + 20]
    y_hat = pdf[idx:idx + 20]
    yy_hat = y_hat[y_hat > 0]
    xx_hat = x_hat[y_hat > 0]
    return xx_hat, yy_hat

def pareto_fit(bx,pdf,increase):

    idx = np.argmax(pdf)
    idx = idx + increase
    xx_hat, yy_hat = extract_subvec(bx, pdf, idx)
    yfit, m, b1, r2, mse = power_law_fit(xx_hat, yy_hat)
    return m, b1

def power_law_fit(xhat, yhat):


    logx = np.log(xhat)
    logy = np.log(yhat)

    m, b = np.polyfit(logx, logy, 1)
    logyreg = b + m * logx
    sst = np.sum( np.power( logy-np.mean(logy), 2))
    ssreg = np.sum( np.power( logyreg - np.mean(logy), 2))
    r2 = ssreg/sst
    mse = np.sum( np.power(logy-logyreg, 2) ) / xhat.size

    yfit = np.exp(b) * np.power(xhat, m)
    
    return yfit, m, b, r2, mse


def create_bins(vals, num_bins, spacing = "log", x=[], weights=None, a_low=None, a_high=None, bin_edge = "center"):
    if a_low == None:
        a_low = 0.95 * np.min(vals)
    if a_high == None:
        a_high = np.max(vals)

    # Create bins 
    if x == []:
        if spacing == "linear":
            x = np.linspace(a_low,a_high,num_bins+1)
        elif spacing == "log":
            if min(a_low, a_high) > 0:
                x = np.logspace(np.log10(a_low), np.log10(a_high), num_bins+1)
            else:
                x = np.logspace(-2, 2, num_bins+1, endpoint=False)
                A = np.max(x)
                B = np.min(x)
                x = (x-A) * (a_low - a_high) / (B-A) + a_high
        else: 
            print("Unknown spacing type. Using Linear spacing")
            x = np.linspace(a_low,a_high,num_bins+1)
    return x

    
def create_pdf(vals, num_bins, spacing = "log", x = [], weights=None, a_low=None, a_high=None, bin_edge = "center"):
    """  create pdf of vals 

    Parameters
    ----------
        vals : array
           array of values to be binned
        num_bins : int
            Number of bins in the pdf
        spacing : string 
            spacing for the pdf, options are linear and log
        x : array
            array of bin edges
        weights :array
            weights corresponding to vals to be used to create a weighted pdf
        a_low : float
            lower value of bin range. If no value provided 0.95*min(vals) is used
        a_high : float
            upper value of bin range. If no value is provided max(vals) is used
        bin_edge: string
            which bin edge is returned. options are left, center, and right

    Returns
    -------
        bx : array
            bin edges or centers (x values of the pdf)
        pdf : array
            values of the pdf, normalized so the Riemann sum(pdf*dx) = 1.
    """

    # Pick bin range 
    if a_low == None:
        a_low = 0.95 * np.min(vals)
    if a_high == None:
        a_high = np.max(vals)

    # Create bins 
    if x == []:
        if spacing == "linear":
            x = np.linspace(a_low,a_high,num_bins+1)
        elif spacing == "log":
            if min(a_low, a_high) > 0:
                x = np.logspace(np.log10(a_low), np.log10(a_high), num_bins+1)
            else:
                x = np.logspace(-2, 2, num_bins+1, endpoint=False)
                A = np.max(x)
                B = np.min(x)
                x = (x-A) * (a_low - a_high) / (B-A) + a_high
        else: 
            print("Unknown spacing type. Using Linear spacing")
            x = np.linspace(a_low,a_high,num_bins+1)

    # Create PDF
    pdf, bin_edges = np.histogram(vals, bins=x, weights=weights, density=False)

    # Return arrays of the same size
    if bin_edge == "left":
        return bin_edges[:-1],pdf

    elif bin_edge == "right":
        return bin_edges[1:],pdf

    elif bin_edge == "center":
        bx = bin_edges[:-1] + 0.5*np.diff(bin_edges)
        return bx, pdf

    else: 
        print("Unknown bin edge type {0}. Returning left edges".format(bin_edge))
        return bin_edge[:-1],pdf


# # In[13]:


# ## flux weigthing
# #path = "/Users/jhyman/Projects/BES/exp_2_dfn/data/dfn_results/dfn_data_cm"

# # Need to find min and max of WHOLE dataset for plotting

# #a_min = 1e-02
# #a_max = 100
# #for i in range(30):
# #    data1 = load_partime_file(f"/home/msweeney2796/dfnWorks/examples/newBES/backbone_1.0e-09_2.0e+01_paper{i}/dfnTrans_full_dump/partime")
# #    if np.min(data1["total travel time"]) < a_min:
# #        a_min = np.min(data1["total travel time"])
# #    if np.max(data1["total travel time"]) > a_max:
# #        a_max = np.max(data1["total travel time"])
# #print(a_min,a_max)

# fig,ax = plt.subplots(figsize=(10,8)); 
# num_bins = 200
# increase = 10
# btc = np.zeros((10,num_bins))
# j = 1
# for i in range(10):
#     data = load_partime_file(f"/home/msweeney2796/dfnWorks/examples/BES_TPL/output_1.0e-09_2.0e+01_paper{i}/dfnTrans_dump/partime")
#     #data = load_partime_file(f"/lclscratch/msweeney/BES_TPL/output_big_1.0e-09_2.0e+01_paper{i}/dfnTrans_dump/partime")
#     bx,pdf = create_pdf(data["total travel time"],
#                         weights = data["flux weights"], 
#                         num_bins = num_bins)
#     btc[j-1,:] = pdf
#     j += 1
# vals,low,high = mean_min_max(btc)
# ax.loglog(bx,vals,colors[0],alpha=1,linewidth=2)
# ax.loglog(bx,low,colors[0],alpha=0.5,linewidth=2)
# ax.loglog(bx,high,colors[0],alpha=0.5,linewidth=2)
# ax.fill_between(bx,low,high,color=colors[0],alpha = 0.5)
# btc = np.zeros((10,num_bins))
# j = 1
# for i in range(10):
#     data = load_partime_file(f"/home/msweeney2796/dfnWorks/examples/BES_TPL/output_1.0e-09_5.0e+00_paper{i}/dfnTrans_dump/partime")
#     #data = load_partime_file(f"/lclscratch/msweeney/BES_TPL/output_big_1.0e-09_5.0e+00_paper{i}/dfnTrans_dump/partime")
#     bx,pdf = create_pdf(data["total travel time"],
#                         weights = data["flux weights"], 
#                         num_bins = num_bins)
#     btc[j-1,:] = pdf
#     j += 1
# vals,low,high = mean_min_max(btc)
# ax.loglog(bx,vals,colors[1],alpha=1,linewidth=2)
# ax.loglog(bx,low,colors[1],alpha=0.5,linewidth=2)
# ax.loglog(bx,high,colors[1],alpha=0.5,linewidth=2)
# ax.fill_between(bx,low,high,color=colors[1],alpha = 0.5)
# #ax.loglog(bx,pdf,colors[0])

# ax.set_xlabel('Time (years)',fontsize=16)
# ax.set_ylabel('PDF',fontsize=16)
# #ax.set_title(r'$\alpha = 3.5; \kappa = 100.0$',fontsize=16)
# ax.tick_params(axis='x', labelsize=16)
# ax.tick_params(axis='y', labelsize=16)
# #ax.set_xlim(2e-3,7e-2)
# #ax.set_ylim(1e-3,1e3)
# ax.grid(True)
# def bin_and_create_pdf(t_exit, nBins):
#     a_low = 0.9 * np.min(t_exit)
#     a_high = np.max(t_exit)
#     x = create_log_bins(nBins, a_low, a_high)
#     pdf, bin_edges = np.histogram(t_exit, bins=x, density=True)
#     bx = bin_edges[:-1] + np.diff(bin_edges)
#     return bx, pdf

# def create_log_bins(nBins, a_low, a_high):
#     if min(a_low, a_high) > 0:
#         import math
#         return np.logspace(math.log10(a_low), math.log10(a_high), nBins+1)
#     x = np.logspace(-8, 2, nBins+1, endpoint=False)
#     A = np.max(x)
#     B = np.min(x)
#     x = (x-A) * (a_low - a_high) / (B-A) + a_high
#     return x

# partime_0 = "/home/msweeney2796/dfnWorks/examples/BES_TPL/output_1.0e-09_0.0e+00_paper0/dfnTrans_dump/partime"
# #partime_0 = "/lclscratch/msweeney/BES_TPL/output_big_1.0e-09_0.0e+00_paper0/dfnTrans_dump/partime"

# times_0 = np.genfromtxt(partime_0,skip_header=1)
# bx,pdf = bin_and_create_pdf(times_0[:,2], 200)
# norm_time = bx[np.argmax(pdf)]
# plt.loglog(bx,pdf,'r')
# clines = [Line2D([0],[0], color='g',lw=1.5),
#             Line2D([0],[0],color='b',lw=1.5),
#             Line2D([0],[0],color='r',lw=1.5)]
# ax.legend(clines,['High Heterogeneity','Low Heterogeneity','Constant'],fontsize=16,loc='upper right')
# fig.savefig('new3.png')

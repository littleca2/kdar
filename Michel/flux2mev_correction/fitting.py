import numpy as np
from math import factorial
import ROOT
import sys
import matplotlib.pyplot as plt
from scipy import optimize
#from scipy.optimize import curve_fit
import math

MICHEL_E_ENDPOINT = 53.3


# From https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-for-a-dataset
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except (ValueError, msg):
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def GetData(h):
        nE = h.GetNbinsX()
        x, y = zip(*[(h.GetBinCenter(i+1), h.GetBinContent(i+1)) for i in range(h.GetNbinsX()) if h.GetBinContent(i+1) != 0])
        err = [np.sqrt(v) for v in y]
        return x,y,err


class MCMichelSpectrum():
    def __init__(self):
        x,y,err = np.loadtxt("/home/marzece/KDAR_Analysis/MichelAnalysis/AlexMichel/mc_michel_flux_spectrum.txt")
        ynew = savitzky_golay(y, 9,3)
        THRESHOLD = 10
        FLUX_MIN = 10e3 # This value was determined from by-eye inspection
        x,y,ynew, err = zip(*[(_x, _y, _yp, _err) for _x,_y, _yp, _err in zip(x,y,ynew,err) if _y > THRESHOLD and _x >= 10e3]) 
        # Have the default energy scale of the MC spectrum go to 1.0 as the max value
        self.xax = x/x[-1]
        # And normalize the yvalues such that they integrate to 1.0
        self.y = np.array(y)/np.trapz(y, x=self.xax)
        self.mc_err = np.array(err)/np.trapz(y, x=self.xax)

        # The prime values are just a smoothed version of the non-prime values
        self.yprime = np.array(ynew)/np.trapz(ynew, x=self.xax)
        self.mc_err_prime = err/np.trapz(ynew, x=self.xax)

    def get_val(self, x, ampl, scale, res, res_scale=0.5):
        y_vals = ampl*self.yprime
        xax = scale*self.xax
        bin_width = np.diff(xax)[0]

        # res represents the fractional error at the endpoint
        # The resolution scales like sqrt(E)
        res = np.power(self.xax, res_scale)*res
        result = np.zeros_like(x)
        sqrt2pi = np.sqrt(2*np.pi)
        for _x, y, sigma in zip(xax, y_vals, res):
            sigma = x*sigma
            coef = 1.0/(sigma*sqrt2pi)
            myGauss = coef*np.exp(-0.5*((x-_x)**2)/sigma**2)
            result += myGauss*y*bin_width
        result[x < xax[0]] = 0
        return result

    def do_fit(self, hist, full_errors=False):
        # Fit the given histogram with a Michel energy spectrum for the energy scale
        # and energy resolution (and amplitude).
        # Then return the fit values & errors and a TGraph for the best firt function
        E_vals, N_vals, err = GetData(hist)
        max_val = hist.GetBinCenter(hist.GetMaximumBin())
        magic_factor = 1.5 # Trial and error good guess value
        scale_guess = max_val*magic_factor
        # Fit parameters are the normalization, energy scale, and energy resolution
        guess = [hist.GetMaximum(), scale_guess , 0.04]
        low_b = [0, 0.1*scale_guess, 1e-5]
        upp_b = [5*hist.GetMaximum() , 10.0*scale_guess, 1.0]
        fit_vals, fit_err = curve_fit(self.get_val, E_vals, N_vals, p0=guess, sigma=err, bounds=(low_b, upp_b))
        if not full_errors:
            fit_err = np.sqrt(np.diag(fit_err))
        graph = ROOT.TGraph()
        for i, (x,y) in enumerate(zip(E_vals, self.get_val(np.array(E_vals), fit_vals[0], fit_vals[1], fit_vals[2]))):
            graph.SetPoint(i, x, y)
        # Re-order the values b/c I'm stupidly put them in the wrong order
        ret = [fit_vals[1], fit_vals[2], fit_vals[0]]
        ret_err = [fit_err[1], fit_err[2], fit_err[0]]
        return ret, ret_err, graph

def Res(E, R_ep):
        E_ep = MICHEL_E_ENDPOINT
        p = 0.02
        x = E/E_ep
        sigma = x*np.sqrt( ((R_ep**2 - p**2)/x) + p**2 )
        return sigma*(E_ep/E)

def Michel(x, par0, par1, par2):
        zScale = par0
        zRes = par1
        zAmpli = par2

        endpoint = MICHEL_E_ENDPOINT/zScale
        x = x/endpoint
        sums = 0
        NPOINTS = 500
        xax = 0.5/NPOINTS + np.linspace(0, 1, NPOINTS+1)[1:]
        for xAux in xax:
                yAux = (3.0 - 2.0*xAux)*xAux**2
                sigma = xAux*np.sqrt((zRes*zRes - 0.0004)/xAux + 0.0004)
                myGauss = np.exp(-0.5*((xAux - x)**2)/sigma**2 )/sigma
                sums += myGauss*yAux
        fitVal = sums*zAmpli*0.39894228*0.002
        return fitVal

def Michel_Flux(x, par0, par1, par2, par3):
        zScale = par0
        zRes = par1
        zAmpli = par2
        Ep = par3

        endpoint = Ep/zScale
        x = x/endpoint
        sums = 0
        NPOINTS = 500
        xax = 0.5/NPOINTS + np.linspace(0, 1, NPOINTS+1)[1:]
        for xAux in xax:
                yAux = (xAux**2)*(3.0 - 2.0*xAux)
                sigma = xAux*np.sqrt((zRes*zRes - 0.0004)/xAux + 0.0004)
                myGauss = np.exp(-0.5*((xAux - x)**2)/sigma**2 )/sigma
                sums += myGauss*yAux
        fitVal = sums*zAmpli*0.39894228*0.002
        return fitVal

def Michel_Conv(x, par0, par1, par2):
        #zScale = par0
        zRes = par0
        zAmpli = par1
        #Ep = par2
        
        endpoint = par2
        x = x/endpoint
        sums = 0
        NPOINTS = 500
        xax = 0.5/NPOINTS + np.linspace(0, 1, NPOINTS+1)[1:]
        for xAux in xax:
                yAux = (xAux**2)*(3.0 - 2.0*xAux)
                sigma = xAux*np.sqrt((zRes*zRes - 0.0004)/xAux + 0.0004)
                myGauss = np.exp(-0.5*((xAux - x)**2)/sigma**2 )/sigma
                sums += myGauss*yAux
                xAux += 0.002
        fitVal = sums*zAmpli*0.39894228*0.002
        return fitVal

def Michel_S(x, par0, par1, par2, par3): # Sensitivity to p1 parameter

        zScale = par0
        zRes = par1
        zAmpli = par2

        endpoint = MICHEL_E_ENDPOINT/zScale
        x = x/endpoint
        sums = 0
        NPOINTS = 500
        xax = 0.5/NPOINTS + np.linspace(0, 1, NPOINTS+1)[1:]
        for xAux in xax:
                yAux = (xAux**2)*(3.0 - 2.0*xAux)
                sigma = xAux*np.sqrt((zRes*zRes - par3*par3)/xAux + par3*par3)
                myGauss = np.exp(-0.5*((xAux - x)**2)/sigma**2 )/sigma
                sums += myGauss*yAux
                xAux += 0.002
        fitVal = sums*zAmpli*0.39894228*0.002
        return fitVal

class Muons:
    def __init__(self, nbins):
        mum_vals = np.loadtxt("/home/littleca/kdar/Michel/MC/MC_mu_minus_FV_edep_vals.nptxt")
        mup_vals = np.loadtxt("/home/littleca/kdar/Michel/MC/MC_mu_plus_FV_edep_vals.nptxt") 
        hmuminus = ROOT.TH1D("", "", nbins, 0, 100)
        hmuplus = ROOT.TH1D("", "", nbins, 0, 100)
        for v in mup_vals:
            hmuplus.Fill(v)

        for v in mum_vals:
            hmuminus.Fill(v)

        self.hcombined = ROOT.TH1D(hmuplus)
        self.hcombined.Add(hmuminus)

        self.muplus_dist = np.array([hmuplus[i] for i in range(hmuplus.GetNbinsX())])
        self.muminus_dist = np.array([hmuminus[i] for i in range(hmuminus.GetNbinsX())])
        #self.muplus_dist = np.array([hmuplus[i+1] for i in range(hmuplus.GetNbinsX()) if hmuplus.GetBinLowEdge(i+1) >= 20.0 and hmuplus.GetBinCenter(i+1) < 60])
        #self.muminus_dist = np.array([hmuminus[i+1] for i in range(hmuminus.GetNbinsX())if hmuminus.GetBinLowEdge(i+1) >= 20.0 and hmuminus.GetBinCenter(i+1) < 60])

def gaus(x, mu, sigma):
    """Returns the normalized gaussian value at "x" for gaussian with mean & width of mu and sigma."""
    coef = 1.0/(sigma*np.sqrt(2*np.pi))
    expo = -0.5 * ((x-mu)/sigma)**2
    return coef*np.exp(expo)

def get_resolution(E, sigma_ep=0.03, constant_term=0.0):
    """ Returns the RELATIVE energy resolution at the given energy E for
    the scaling & constant term specified at the endpoint energy """
    efrac = E/MICHEL_E_ENDPOINT
    return np.sqrt(sigma_ep**2/efrac + constant_term**2)

def apply_smearing(edep_spectrum, evals, a, c):
    sigmas = [get_resolution(x, a, c)*x for x in evals]
    bin_width = np.diff(evals)[0]
    smears = np.array([gaus(evals, x, s) for x,s in zip(evals, sigmas)])
    ret = np.sum([scale*smear*bin_width for scale, smear in zip(edep_spectrum, smears)],axis=0)
    return ret

def edep_fit(args, data_bins, binned_data, set_hist, muons_info) :
    scale = args[0]
    escale = math.exp(args[1])
    smear1 = args[2]
    smear2 = args[3]
    mixing_frac = args[4]
    bckg = args[5]
    #print(scale, escale,smear1, smear2, mixing_frac, bckg)

    mc_bin_centers = (np.linspace(0, 100, set_hist["n_bins"]*2) + 0.5)[:-1]
    mc_pred = muons_info.muplus_dist*mixing_frac + muons_info.muminus_dist*(1-mixing_frac)
    mc_pred = (mc_pred*scale) + bckg

    smeared_dist = apply_smearing(mc_pred, mc_bin_centers, smear1, smear2)


    mc_bin_centers_scaled = mc_bin_centers * escale
    mc_pred_w_data_binning_idx = np.digitize(mc_bin_centers_scaled, data_bins)
    mc_pred_w_data_binning = np.zeros(len(data_bins)-1)
    for i, bin_id in enumerate(mc_pred_w_data_binning_idx) :
        if bin_id == len(data_bins) :
            continue
        mc_pred_w_data_binning[bin_id-1] += mc_pred[i]


    chi2 = abs(np.sum((mc_pred_w_data_binning - binned_data)**2/mc_pred_w_data_binning))
    #print("chi2 = ", chi2)
    return chi2

def do_edep_fit(data_vals, set_hist):

    #data_dist = np.array([hdata[i+1] for i in range(hdata.GetNbinsX())])
    data_dist = data_vals
    # Re-bin the flux data
    data_bin_width = float(set_hist["max"] - set_hist["min"]) / float(set_hist["n_bins"])
    data_bin_centers = np.linspace( (set_hist["min"] + float(data_bin_width)/2), (set_hist["max"] - (float(data_bin_width)/2)), set_hist["n_bins"])
    data_bins = np.linspace(set_hist["min"], set_hist["max"], set_hist["n_bins"]+1)
    data_bin_idx = np.digitize(data_dist, data_bins)
    binned_data = np.zeros(len(data_bins)-1)
    for i, bin_id in enumerate(data_bin_idx):
        # Ignore overflow bin data
        if bin_id == len(data_bins):
            continue
        #binned_data[bin_id-1] += data_dist[i]	# Raw data
        binned_data[bin_id-1] += 1		# Histogram
 
    my_muons = Muons(set_hist["n_bins"]*2)

    meth_n = np.array([None, 'L-BFGS-B', 'Nelder-Mead', 'SLSQP', 'Powell', 'COBYLA'])#, 'TNC', 'BFGS'])
    chi_list = []
    suc_list = []
    #for i, name in enumerate(meth_n) :
    #    print(name)
        #bf_test = optimize.minimize(edep_fit, [np.log(960), 0.03, 0.03, 0.5, 0.0], bounds=[(np.log(500), np.log(1500)), (0.001, 1.0), (0.001, 1.0), (0.001, 1.0), (0.0, None)], method=name, args=(data_bin_centers, binned_data, set_hist, my_muons))
    #    bf_test = optimize.basinhopping(edep_fit, [1.0, np.log(960), 0.05, 0.05, 0.5, 10], stepsize=0.001, minimizer_kwargs={'bounds':[(0.0, None), (0, np.log(1e5)), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, None)], 'method':name, 'args':(data_bins, binned_data, set_hist, my_muons)})
    #    chi_list.append(bf_test.fun)
    #    suc_list.append(bf_test.success)
    #    print(bf_test)
    #    if bf_test.success and bf_test.x.all()>=0:
    #        if bf_test.fun < chi_list[i-1] and i > 0 :
    #            bf = bf_test
    #        elif i == 0 :
    #            bf = bf_test
    #print([(name, chi, suc) for name, chi, suc in sorted(zip(meth_n, chi_list, suc_list), key=lambda x : x[1])])
    bf = optimize.basinhopping(edep_fit, [1.0, np.log(960), 0.05, 0.05, 0.5, 10], stepsize=0.001, minimizer_kwargs={'bounds':[(0.0, None), (0, np.log(1e5)), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, None)], 'method':'Nelder-Mead', 'args':(data_bins, binned_data, set_hist, my_muons)})
    print(bf)
    #print(bf.hess_inv)
    #bf_err = np.sqrt(np.diag(bf.hess_inv.todense()))
    #bf_err = np.sqrt(np.diag(np.linalg.inv(bf.hess)))
    #print(bf_err)
    bf_err = [[0,0],[0,0]]    

    #my_bf = [21.9708, -0.652794, 0.026418, 0.0, 0.378431, 0]
    #bf.x = np.array(my_bf)

    bf_scale = bf.x[0]
    bf_escale = math.exp(bf.x[1])
    bf_smear1 = bf.x[2]
    bf_smear2 = bf.x[3]
    bf_mixing_frac = bf.x[4]
    bf_bckg = bf.x[5]

    gbf = ROOT.TGraph()

    mc_bin_centers = (np.linspace(0, 100, set_hist["n_bins"]*2) + 0.5)[:-1]
    mc_pred = my_muons.muplus_dist
    mc_pred = my_muons.muplus_dist*bf_mixing_frac + my_muons.muminus_dist*(1-bf_mixing_frac)
    mc_pred = (mc_pred*bf_scale) + bf_bckg

    mc_bin_centers_scaled = mc_bin_centers * bf_escale
    mc_pred_w_data_binning_idx = np.digitize(mc_bin_centers_scaled, data_bins)
    mc_pred_w_data_binning = np.zeros(len(data_bins-1))
    for i, bin_id in enumerate(mc_pred_w_data_binning_idx) :
        mc_pred_w_data_binning[bin_id-1] += mc_pred[i]

    for i, (x,v) in enumerate(zip(mc_bin_centers_scaled, mc_pred_w_data_binning)):
        gbf.SetPoint(i, x, v)

    return bf.x, bf_err, gbf

def do_michel_fit(hist, full_errors=False):
    # Fit the given histogram with a Michel energy spectrum for the energy scale
    # and energy resolution (and amplitude).
    # Then return the fit values & errors and a TGraph for the best firt function
    E_vals, N_vals, err = GetData(hist)
    max_val = hist.GetBinCenter(hist.GetMaximumBin())
    scale_guess = MICHEL_E_ENDPOINT/max_val
    #scale_guess = max_val
    # Fit parameters are the energy scale, energy resolution and normalization
    guess = [scale_guess , 0.04, hist.GetMaximum()]
    low_b = [0.1*scale_guess, 0.001, 0]
    upp_b = [10.0*scale_guess, 1.0, 2*hist.GetMaximum()]
    fit_vals, fit_err = curve_fit(Michel, E_vals, N_vals, p0=guess, sigma=err, bounds=(low_b, upp_b))
    if not full_errors:
        fit_err = np.sqrt(np.diag(fit_err))
    g = ROOT.TGraph()
    xax = np.linspace(hist.GetBinLowEdge(1), hist.GetBinLowEdge(hist.GetNbinsX()) + hist.GetBinWidth(hist.GetNbinsX()) , 1000)
    yvals = Michel(xax, fit_vals[0], fit_vals[1], fit_vals[2])
    for i, (x,y) in enumerate(zip(xax, yvals)):
        g.SetPoint(i, x, y)
    return fit_vals, fit_err, g

def do_michel_fit_conv(hist, full_errors=False):
    # Fit the given histogram with a Michel energy spectrum for the energy scale
    # and energy resolution (and amplitude).
    # Then return the fit values & errors and a TGraph for the best firt function
    E_vals, N_vals, err = GetData(hist)
    max_val = hist.GetBinCenter(hist.GetMaximumBin())
    scale_guess = max_val
    #scale_guess = max_val
    # Fit parameters are the energy scale, energy resolution and normalization
    low_b = [0.001, 0, 0.1*scale_guess]
    upp_b = [1.0, 2*hist.GetMaximum(), 10.0*scale_guess]

    guess = [0.04, hist.GetMaximum(), scale_guess]
    fit_vals, fit_err = curve_fit(Michel_Conv, E_vals, N_vals, p0=guess, sigma=err, bounds=(low_b, upp_b))
    if not full_errors:
        fit_err = np.sqrt(np.diag(fit_err))
    return fit_vals, fit_err

class PositionBasedHistograms():
    def __init__(self, name, title, nx, xlow, xhigh, ny, ylow, yhigh, nz, zlow, zhigh):
        name1 = name+"_h2" if name else ""
        name2 = name+"_h2_fit_scale" if name else ""
        name3 = name+"_h2_fit_resolution" if name else ""
        self.h2 = ROOT.TH2D( name, "%s;R^{2} [mm^{2}];Z [mm]" % title, nx, xlow, xhigh, ny, ylow, yhigh)
        self.scale_h2 = ROOT.TH2D(name2, "%s;R^{2} [mm^{2}];Z [mm]" % title, nx, xlow, xhigh, ny, ylow, yhigh)
        self.resolution_h2 = ROOT.TH2D(name3, "%s;R^{2} [mm^{2}];Z [mm]" % title, nx, xlow, xhigh, ny, ylow, yhigh)
        self.fit_vals = None
        self.hists = [ROOT.TH1D(name+"_h1_%i"%idx if name else "", "%s;R^{2} [mm^{2}];Z [mm]" % title, nz, zlow, zhigh) for idx in range(nx*ny)]
        self.nx = nx
        self.ny = ny
        self.lines = []
        self.graphs = []
    def Fill(self, x, y, z, weight=None):
        xbin = self.h2.GetXaxis().FindBin(x) - 1
        ybin = self.h2.GetYaxis().FindBin(y) - 1
        if(xbin < 0 or ybin < 0 or xbin >= self.nx or ybin >= self.ny):
            return
        if weight is  None:
            self.hists[ybin*self.nx + xbin].Fill(z)
        else:
            self.hists[ybin*self.nx + xbin].Fill(z, weight)

    def Draw(self, canvas, fitf=None):
        #canvas.cd()
        #self.h2.Draw("col")
        self.ll1 = ROOT.TLine(0, 1000, (1400/1850)**2, 1000)
        self.ll2 = ROOT.TLine((1400/1850)**2, 1000, (1400/1850)**2, -1000)
        self.ll3 = ROOT.TLine(0, -1000, (1400/1850)**2, -1000)
        canvas.Divide(self.nx, self.ny, 0.0, 0.0, 0)
        # Canvas/pad indices go from left to right, top to bottom.
        # TH2D indicies go from left to right, bottom to top
        for ix in range(self.nx):
            for iy in range(self.ny):
                canv_iy = (self.ny-1)-iy  # Reverse the vertical index to get the top pad
                num = canv_iy*self.nx + ix
                pad = canvas.cd(num+1)
                h = self.hists[iy*self.nx+ix]
                h.SetStats(0)
                h.SetTitle("")
                h.SetLineColor(4)
                h.Draw("HIST")
                if(fitf):
                    graph = fitf(h)
                    #_, _, graph = mc_michel.do_fit(h)
                    if(graph):
                        graph.SetLineColor(2)
                        graph.SetLineWidth(2)
                        self.graphs.append(graph)
                        graph.Draw("sameline")
                ll = ROOT.TLine(MICHEL_E_ENDPOINT, 0, MICHEL_E_ENDPOINT, pad.GetUymax())
                ll.SetLineColor(2)
                ll.SetLineWidth(2)
                self.lines.append(ll)
                ll.Draw()
        canvas.Update()

    def Fit(self):
        self.fit_vals = [do_michel_fit_conv(h) if h.GetEntries() >0 else None for h in self.hists]
        #self.fit_vals = [mc_michel.do_fit(h) for h in self.hists]

        # I want all the energy scale values relative to particular bin
        #scale_benchark = self.fit_vals[12*5 + 0][0][2];
        scale_benchark = MICHEL_E_ENDPOINT

        for ix in range(self.nx):
            for iy in range(self.ny):
                hist_idx = (iy+1)*self.nx + (ix+1)
                arr_idx = iy*self.nx + ix

                if(self.fit_vals[arr_idx] is None):
                    continue
                # The energy scale is reported by the fitter in weird way where
                # if energy scale is low, the value reported is high.
                # So take the inverse of the value so it becomes more sensible,
                # and then propagate that change to the error too
                scale_val = self.fit_vals[arr_idx][0][2]
                fractional_scale_err = self.fit_vals[arr_idx][1][2]/scale_val
                #scale_val = 1.0/scale_val
                scale_val /= scale_benchark;
                scale_err = fractional_scale_err*scale_val
                self.scale_h2.SetBinContent(ix+1, iy+1, scale_val)
                self.scale_h2.SetBinError(ix+1, iy+1, scale_err)

                self.resolution_h2.SetBinContent(ix+1, iy+1, self.fit_vals[arr_idx][0][0])
                self.resolution_h2.SetBinError(ix+1, iy+1, self.fit_vals[arr_idx][1][0])

    def DrawFit(self, c1, c2):
        if not self.fit_vals:
            print("Fit values not available!")
            return
        ll1 = ROOT.TLine(0, 1e3, 1.4e3**2, 1e3)
        ll2 = ROOT.TLine(1.4e3**2, 1e3, 1.4e3**2, -1e3)
        ll3 = ROOT.TLine(0, -1e3, 1.4e3**2, -1e3)
        [ll.SetLineColor(2) for ll in [ll1, ll2, ll3]]
        [ll.SetLineWidth(2) for ll in [ll1, ll2, ll3]]
        [self.lines.append(ll) for ll in [ll1, ll2, ll3]]
        c1.cd()
        self.scale_h2.SetStats(0)
        self.scale_h2.Draw("colztextE")
        [ll.Draw("same") for ll in [ll1, ll2, ll3]]
        c2.cd()
        self.resolution_h2.SetStats(0)
        self.resolution_h2.Draw("colztextE")
        [ll.Draw("same") for ll in [ll1, ll2, ll3]]

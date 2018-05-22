#!/usr/bin/env python
# ianh@astro.ox.ac.uk


import matplotlib
matplotlib.use('Agg')
import numpy
import pylab
import sys
from astropy.io import fits
from time import time,gmtime,strftime
from optparse import OptionParser


# ---------------------------------------------------------------------------------------
# Command line options
# ---------------------------------------------------------------------------------------


parser = OptionParser(usage='%prog [options] rms_map')
parser.add_option('--s_min',dest='s_min',help='Lower flux density limit for count determination [Jy] (default = 20e-6)',default=20e-6)
parser.add_option('--s_max',dest='s_max',help='Upper flux density limit for count determination [Jy] (default = 0.1)',default=0.1)
parser.add_option('--n_bins',dest='n_bins',help='Number of logarithmically-spaced bins between s_min and s_max (default = 25)',default=25)
parser.add_option('--sigma',dest='sigma',help='Sigma cut above which a source is considered to be detected (default = 5)',default=5.0)
parser.add_option('--n_iter',dest='n_iter',help='Number of iterations to simulate noisy detections (default = 1)',default=1)
parser.add_option('--s_cubed',dest='s_cubed',help='S-cubed query .result file (default = e207dfe12f12c123f003651be8c224a626e9907d.result)',
        default='e207dfe12f12c123f003651be8c224a626e9907d.result')
parser.add_option('--area',dest='area',help='Area of the S-cubed query [square degrees] (default = 100)',default=100.0)
parser.add_option('--col',dest='col',help='Flux density column to use from S-cubed query (default = itot_1400)',default='itot_1400')
parser.add_option('--doplot',dest='doplot',help='Plot the true and simulated source counts (default = True)',action='store_false',default=True)
parser.add_option('--pngname',dest='pngname',help='Name of output PNG file, output table will be written to a txt file of the same name',default='')


(options,args) = parser.parse_args()
s_min = float(options.s_min)
s_max = float(options.s_max)
n_bins = int(options.n_bins)
sigma = float(options.sigma)
n_iter = int(options.n_iter)
s_cubed = options.s_cubed
area = float(options.area)
col = options.col
doplot = options.doplot
pngname = options.pngname


if len(args) != 1:
    print 'Please specify a FITS RMS image.'
    sys.exit()
else:
    rms_map = args[0].rstrip('/')


if pngname == '':
    pngname = 'eddington_'+str(time()).split('.')[0]+'.png'


# ---------------------------------------------------------------------------------------
# Function definitions
# ---------------------------------------------------------------------------------------


def msg(message):
    tt = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    print tt+'   '+message


def centralize(bins):
    bincentres = 0.5*(bins[1:]+bins[:-1])
    return bincentres


def find_minimum(infits):
    inphdu = fits.open(infits)[0]
    if len(inphdu.data.shape) == 2:
        image = numpy.array(inphdu.data[:,:])
    elif len(inphdu.data.shape) == 3:
        image = numpy.array(inphdu.data[0,:,:])
    else:
        image = numpy.array(inphdu.data[0,0,:,:])
    image[image==0.0] = numpy.nan 
    return image.ravel()[numpy.nanargmin(image)]


def get_area(infits,rms):
    inphdu = fits.open(infits)[0]
    hdr = inphdu.header
    dx = hdr.get('CDELT1')
    dy = hdr.get('CDELT2')
    pixsize = abs(dx*dy)
    if len(inphdu.data.shape) == 2:
        image = numpy.array(inphdu.data[:,:])
    elif len(inphdu.data.shape) == 3:
        image = numpy.array(inphdu.data[0,:,:])
    else:
        image = numpy.array(inphdu.data[0,0,:,:])
    mask = (image < rms) & (~numpy.isnan(image))
    nonblanks = len(image[mask])
    a_sq_deg = nonblanks*pixsize
    a_sr = a_sq_deg/3283.0
    return a_sq_deg,a_sr,rms


def add_noise(infits,fluxes):
#    noise = numpy.zeros(len(fluxes))
    noise = []
    thresh = []
    inphdu = fits.open(infits)[0]
    if len(inphdu.data.shape) == 2:
        image = numpy.array(inphdu.data[:,:])
    elif len(inphdu.data.shape) == 3:
        image = numpy.array(inphdu.data[0,:,:])
    else:
        image = numpy.array(inphdu.data[0,0,:,:])
    image[image==0.0] = numpy.nan
    hdr = inphdu.header
    nx = hdr.get('NAXIS1')
    ny = hdr.get('NAXIS2')
    for i in range(0,len(fluxes)):
        rms = numpy.nan
        while numpy.isnan(rms):
            xx = numpy.random.randint(low=0,high=nx)
            yy = numpy.random.randint(low=0,high=ny)
            rms = image[yy,xx]
        noise.append(numpy.random.normal(loc=0.0,scale=rms))
        thresh.append(rms)
    noise = numpy.array(noise)
    thresh = numpy.array(thresh)
    noisy_fluxes = fluxes + noise
    return noise,thresh,noisy_fluxes


def getS3(result):
    f = open(result,'r')
    line = f.readline()
    cols = line.split(',')
    f.close()
    idx = numpy.where(numpy.array(cols)==col)[0][0]
    s3 = numpy.loadtxt(fname=result,skiprows=1,delimiter=',',usecols=idx)
    return 10.0**s3


def match(artist):
    return artist.__module__ == 'matplotlib.text'


# ---------------------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------------------


def main():
    msg('Reading '+s_cubed)
    msg('This might take a minute or two...')
    s3_true = getS3(s_cubed)
    msg('Catalogue contains '+str(len(s3_true))+' sources')

    min_rms = find_minimum(rms_map)
    msg('Minimum RMS noise '+str(min_rms))
    mask = s3_true > 0.5*min_rms
    s3_true = s3_true[mask]
    n_src = len(s3_true)
    msg('Catalogue reduces to '+str(n_src)+' sources above 0.5 times minimum RMS')

    my_bins = numpy.logspace(numpy.log10(s_min),numpy.log10(s_max),n_bins)

    true_counts, bin_edges = numpy.histogram(s3_true,bins=my_bins)
    true_counts = true_counts / area

    bin_centres = centralize(my_bins)
    bin_widths = bin_edges[1:] - bin_edges[:-1]

    true_euclidean = numpy.array((bin_centres**2.5) * (true_counts/bin_widths))

    noisy_euclidean = numpy.zeros(len(bin_centres))



    fEs = numpy.zeros(len(bin_centres))
    msg('Simulating noisy counts')

    for i in range(0,n_iter):
        msg('----- Iteration: '+str(i+1))

        noise,thresh,s3_noisy = add_noise(rms_map,s3_true)
        mask = s3_noisy > sigma*thresh
        noisy_detected = s3_noisy[mask]

        noisy_counts, bin_edges = numpy.histogram(noisy_detected,bins=my_bins)
        noisy_counts = noisy_counts / area

        temp_counts = numpy.array((bin_centres**2.5) * (noisy_counts/bin_widths))
        noisy_euclidean = numpy.add(noisy_euclidean,temp_counts)
        temp_fE = temp_counts / true_euclidean
        fEs = numpy.add(fEs,temp_fE)

    noisy_euclidean = noisy_euclidean / float(n_iter)
    fEs = fEs / float(n_iter)

    Ars = []
    Ads = []

    for i in range(0,len(bin_centres)):
        Ad,Ar,r = get_area(rms_map,bin_centres[i]/sigma)
        Ars.append(Ar)
        Ads.append(Ad)

    if doplot:
        fig = pylab.figure(figsize=(18,12))
        ax = fig.add_subplot(111,axisbg='#EEEEEE')
        ax.plot(bin_centres,true_euclidean,'.-',color='blue',zorder=100)
        ax.plot(bin_centres,noisy_euclidean,'.-',color='red',zorder=100)
        lim = sigma*min_rms
        ax.axvspan(lim,lim,0,1,color='black',alpha=0.5,zorder=100)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('log$_{10}$ S [Jy]')
        ax.set_ylabel('log$_{10}$ S$^{2.5}$ dN / dS [Jy$^{1.5}$ sr$^{-1}$]')
        ax.grid(b=True,which='minor',color='white',linestyle='-',lw=1)
        ax.grid(b=True,which='major',color='white',linestyle='-',lw=1)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x',which='both',bottom='off',top='off')
        ax.tick_params(axis='y',which='both',left='off',right='off')
        ax.legend(('True counts','Perturbed counts','Deepest limit of survey'),loc='lower right')
        for textobj in fig.findobj(match=match):
            textobj.set_fontsize(22)
        fig.savefig(pngname)
        msg('Rendered '+pngname)

    opfile = pngname.replace('.png','.txt')
    msg('Done. Here are the tabulated results, also written to '+opfile)

    print ''
    print '#------------------------------------------------------------------------------------------------------------------------'
    print '#    Bin_lower   Bin_upper   Bin_centre  True_raw    Det_raw     True_Euclid Det_Euclid  Eddington   Bin_area    Bin_area'
    print '#    Jy          Jy          Jy                                  Jy^1.5 sr-1 Jy^1.5 sr-1             sr          deg^2'
    print '#------------------------------------------------------------------------------------------------------------------------'

    for i in range(0,len(bin_centres)):
        s0 = round(bin_edges[i],7)
        s1 = round(bin_edges[i+1],7)
        sc = round(bin_centres[i],7)
        tc = int(true_counts[i])
        nc = int(noisy_counts[i])
        te = round(true_euclidean[i],4)
        ne = round(noisy_euclidean[i],4)
        fE = round(fEs[i],2)
        Ar = round(Ars[i],6)
        Ad = round(Ads[i],6)
        print '     %-12f%-12f%-12f%-12s%-12s%-12s%-12s%-12s%-12f%-12f' % (s0,s1,sc,tc,nc,te,ne,fE,Ar,Ad)
    print '#------------------------------------------------------------------------------------------------------------------------'
    print ''

    f = open(opfile,'w')
    print >>f,'#------------------------------------------------------------------------------------------------------------------------'
    print >>f,'#    Bin_lower   Bin_upper   Bin_centre  True_raw    Det_raw     True_Euclid Det_Euclid  Eddington   Bin_area    Bin_area'
    print >>f,'#    Jy          Jy          Jy                                  Jy^1.5 sr-1 Jy^1.5 sr-1             sr          deg^2'
    print >>f,'#------------------------------------------------------------------------------------------------------------------------'

    for i in range(0,len(bin_centres)):
        s0 = round(bin_edges[i],7)
        s1 = round(bin_edges[i+1],7)
        sc = round(bin_centres[i],7)
        tc = int(true_counts[i])
        nc = int(noisy_counts[i])
        te = round(true_euclidean[i],4)
        ne = round(noisy_euclidean[i],4)
        fE = round(fEs[i],2)
        Ar = round(Ars[i],6)
        Ad = round(Ads[i],6)
        print >>f,'     %-12f%-12f%-12f%-12s%-12s%-12s%-12s%-12s%-12f%-12f' % (s0,s1,sc,tc,nc,te,ne,fE,Ar,Ad)
    print >>f,'#------------------------------------------------------------------------------------------------------------------------'
    f.close()


if __name__ == "__main__":
    main()
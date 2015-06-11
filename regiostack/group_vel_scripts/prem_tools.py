import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

names = list(["RADIUS", "DEPTH", "DENSITY", "VPV", "VPH", "VSV", "VSH", "ETA",
              "QMU", "QKAPPA"])


def read_prem(fname):
    pd_prem = pd.read_csv(fname, sep=',', header=None, names=names)
    prem = pd_prem[names].values

    return prem

def plot_premVsv(prem, fname):

    plt.figure(figsize=(3, 5))
    plt.plot(prem[:, names.index("VSV")], -1.*prem[:, names.index("DEPTH")])

    plt.ylim([-600, 0])
    plt.xlim([3, 6])
    plt.xlabel('Vsv (km/s)')
    plt.ylabel('Depth (km)')
    plt.title('PREM')

    plt.savefig(fname, dpi=300)
    print fname

def write_Hermann_mod(prem, fname):

    # open the file and write the header information
    f = open(fname, 'w')
    f.write("MODEL.01\n")
    f.write("PREM\n")
    f.write("ISOTROPIC\n")
    f.write("KGS\n")
    f.write("SPHERICAL EARTH\n")
    f.write("1-D\n")
    f.write("CONSTANT VELOCITY\n")
    f.write("LINE08\n")
    f.write("LINE09\n")
    f.write("LINE10\n")
    f.write("LINE11\n")
    f.write("H(KM)   VP(KM/S)   VS(KM/S) RHO(GM/CC)     QP         QS       ETAP       ETAS      FREFP      FREFS\n")

    h = prem[:, names.index("DEPTH")]
    vp = prem[:, names.index("VPV")]
    vs = prem[:, names.index("VSV")]
    rho = prem[:, names.index("DENSITY")]
    q_mu = prem[:, names.index("QMU")]
    q_k = prem[:, names.index("QKAPPA")]
    qp, qs = Q_PS_from_Q_muk(q_mu, q_k, vp, vs)
    freq = 1.0

    prem_layers = depth_to_layers(h, vp, vs, rho, qp, qs)

    nrows, ncols = prem_layers.shape
    for i in xrange(nrows):
        f.write("%8.4f%10.4f%10.4f%10.4f%12.4f%12.4f%8.2f%8.2f%8.2f%8.2f\n"
                % (prem_layers[i, 0], prem_layers[i, 1], prem_layers[i, 2],
                prem_layers[i, 3], prem_layers[i, 4], prem_layers[i, 5], 0.0,
                0.0, freq, freq))

    # close the file
    f.close()
    print fname

def Q_PS_from_Q_muk(q_mu, q_k, vp, vs):
    npts = len(vp)
    QS = q_mu[:]
    QP = q_k[:]
    for i in xrange(npts):
        if q_mu[i] > 0:
            L = (4./3.) * (vs[i]/vp[i]) * (vs[i]/vp[i])
            QP[i] = 1./(L/q_mu[i] + (1-L)/q_k[i])
    return QP, QS

def depth_to_layers(h, vp, vs, rho, qp, qs):

    # get number of depths
    npts = len(vp)

    # set up arrays for layers
    new = np.empty((npts-1,6), dtype=float) 

    # reduce to layers by taking differences in depths and averages of values
    for i in xrange(npts-1):
        new[i, 0] = h[i+1] - h[i]
        new[i, 1] = (vp[i+1] + vp[i])/2.
        new[i, 2] = (vs[i+1] + vs[i])/2.
        new[i, 3] = (rho[i+1] + rho[i])/2.
        new[i, 4] = (qp[i+1] + qp[i])/2.
        new[i, 5] = (qs[i+1] + qs[i])/2.

    # remove 0 thickness layers (indicate discontinuities)
    zero_layers = np.arange(npts-1)[new[:, 0] < 1e-3]
    new_clean = np.delete(new, zero_layers, 0)

    return new_clean

def run_dispersion(mod_fname, pmin, pmax):

    do_Hermann_cleanup()

    # use Hermann computer codes to produce dispersion curve for model
    os.system('sprep96 -M %s -R -PMIN %.1f -PMAX %.1f' % (mod_fname, pmin, pmax))
    os.system('sdisp96')
    os.system('sregn96')
    os.system('sdpegn96 -R -U -PER -ASC')

    # save the dispersion file for PREM
    os.rename('SREGN.ASC', 'PREM_RU0.ASC')

    do_Hermann_cleanup()

    dsp = read_asc('PREM_RU0.ASC')
    return dsp

def do_Hermann_cleanup():

    Hermann_files = ["sdisp96.dat", "sdisp96.ray", "sregn96.egn", "SREGNU.PLT",
                     "SREGN.ASC"]

    # do cleanup
    for f in Hermann_files:
        if os.path.isfile(f):
            os.remove(f)

def read_asc(fname):
    """
    Reads an ascii-formatted theoretical dispersion curve from Hermann code
    """

    # set up list of names
    names = list(["RMODE", "NFREA", "PER", "FREQ", "CVEL", "UVEL", "ENERGY",
                  "GAMMA", "ELLIPT"])
    pd_dsp = pd.read_table(fname, sep="\s+", header=1, names=names)

    dsp_names = list(["PER", "UVEL"])
    dsp = pd_dsp[dsp_names].values

    return dsp


if __name__ == '__main__':

    do_disp = False
    fname = 'PREM_1s.csv'
    prem = read_prem(fname)
    
    plot_premVsv(prem, 'PREM_1s.png')
    write_Hermann_mod(prem, 'prem.mod')
    if do_disp:
        dsp = run_dispersion('prem.mod', 5.0, 250.0)
    else:
        dsp = read_asc('PREM_RU0.ASC')
    print dsp

import logging
import h5py
import numpy as np
import argparse
import glue
import matplotlib
from gwpy.table.lsctables import SnglBurstTable
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--single-ifo-trigs', type=str, required=True,
		help='HDF file cointaining coincident CBC triggers')
parser.add_argument('--template-file', type=str, required=True,
		help='HDF file containing template information for coincident triggers')
parser.add_argument('--ifo', type=str, required=True,
        help='IFO, H1 or L1')
parser.add_argument('--ranking-statistic', type=str, required=True, default='newsnr',
        choices=['snr','newsnr'],help='Ranking statistic, snr or newsnr')
parser.add_argument('--output-file', type=str, required=True,
        help='Full path to output file')
parser.add_argument('--N', type=int, required=True,
		help='Will plot the Nth loudest single ifo trigger')
args = parser.parse_args()

trigs_file = h5py.File(args.single_ifo_trigs,'r')
template_file = h5py.File(args.template_file,'r')
ifo = args.ifo

logging.info('Parsing HDF files')

# get SNR of single detector triggers and find index of loudest event
snr = trigs_file[ifo]['snr'][:]
highestSNR = np.argsort(snr)[-args.N]
chisq = trigs_file[ifo]['chisq'][:]
chisq_dof = trigs_file[ifo]['chisq_dof'][:]
reduced_chisq = chisq/(2*chisq_dof - 2)
# calculate newsnr from reduced chi squared and SNR
newsnr = np.array(snr)
index = np.where(reduced_chisq > 1.0)[0]
newsnr[index] *= (0.5*(1. + reduced_chisq[index]**3))**(-1./6)

highestnewSNR = np.argsort(newsnr)[-args.N]
end_times = trigs_file[ifo]['end_time'][:]
template_ids = trigs_file[ifo]['template_id'][:]

if args.ranking_statistic == 'snr':
    trig_time = end_times[highestSNR]
    trig_template = template_ids[highestSNR]
elif args.ranking_statistic == 'newsnr':
    trig_time = end_times[highestnewSNR]
    trig_template = template_ids[highestnewSNR]

mass1vec = template_file['mass1'][:]
mass2vec = template_file['mass2'][:]
mass1 = mass1vec[trig_template]
mass2 = mass2vec[trig_template]

hoft_chan = str(ifo) + ":GDS-CALIB_STRAIN"
start_time = trig_time - 40
end_time = trig_time + 40

logging.info('Fetching omicron triggers')
omicron_trigs = SnglBurstTable.fetch(hoft_chan,'omicron',start_time,end_time)
omicron_times = np.add(omicron_trigs.get_column('peak_time'),omicron_trigs.get_column('peak_time_ns')*10**-9)
omicron_snr = omicron_trigs.get_column('snr')
omicron_freq = omicron_trigs.get_column('peak_frequency')

from pycbc.workflow.segment import fromsegmentxml
import pycbc.pnutils

f_low = 30
f_high = pycbc.pnutils.f_SchwarzISCO(mass1+mass2)

inspiral_t, inspiral_f = pycbc.pnutils.get_inspiral_tf(trig_time, mass1, mass2, f_low, f_high)

logging.info('Plotting')

cm = plt.cm.get_cmap('jet')
plt.scatter(omicron_times,omicron_freq,c=omicron_snr,s=30,cmap=cm,linewidth=0)
plt.grid(b=True, which='both')
cbar = plt.colorbar()
cbar.set_label('Omicron trigger SNR')
plt.yscale('log')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(start_time,end_time)
if args.ranking_statistic == 'snr':
    plt.suptitle('CBC trigger SNR = ' + str(snr[highestSNR]) + 
                    ", newSNR = " + str(newsnr[highestSNR]),fontsize=12)
elif args.ranking_statistic == 'newsnr':
    plt.suptitle('CBC trigger SNR = ' + str(snr[highestnewSNR]) + 
                    ", newSNR = " +str(newsnr[highestnewSNR]),fontsize=12)
plt.title(str(mass1) + " - " + str(mass2) + " solar masses at GPS time " + str(trig_time),fontsize=12)
plt.hold(True)
plt.plot(inspiral_t,inspiral_f)
plt.savefig(args.output_file)


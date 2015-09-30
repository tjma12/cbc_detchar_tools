import logging
import h5py
import numpy as np
import argparse
import glue
import matplotlib
from gwpy.table.lsctables import SnglBurstTable
from gwpy.table.lsctables import SnglInspiralTable
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--single-ifo-trigs', type=str, required=True,
		help='XML file cointaining coincident CBC triggers or cache file')
parser.add_argument('--ifo', type=str, required=True,
        help='IFO, H1 or L1')
parser.add_argument('--ranking-statistic', type=str, required=True, default='newsnr',
        choices=['snr','newsnr'],help='Ranking statistic, snr or newsnr')
parser.add_argument('--output-file', type=str, required=True,
        help='Full path to output file')
parser.add_argument('--central-time',type=float,required=False,
		help='Central time to look at')
parser.add_argument('--window',type=float,required=False,default=40.0,
		help='Script will find loudest trigger within +/- window seconds')
parser.add_argument('--plot-window',type=float,required=True,default=40.0,
		help='Plot will display +/- plot-window seconds around the loudest trigger')
args = parser.parse_args()

events = SnglInspiralTable.read(args.single_ifo_trigs)
ifo = args.ifo

logging.info('Parsing XML files')

# get SNR of single detector triggers and find index of loudest event
snr = events.get_column('snr')
newsnr = [row.get_new_snr() for row in events]
end_time = events.get_column('end_time')
end_time_ns = events.get_column('end_time_ns')
end_times = np.add(end_time,end_time_ns*10**-9)

if args.ranking_statistic == 'snr':
    highest_idx = np.argsort(snr)[-1]
    trig_time = end_times[highest_idx]
elif args.ranking_statistic == 'newsnr':
    highest_idx = np.argsort(newsnr)[-1]
    trig_time = end_times[highest_idx]

mass1 = events[highest_idx].mass1
mass2 = events[highest_idx].mass2


hoft_chan = str(ifo) + ":GDS-CALIB_STRAIN"
start_time = trig_time - args.plot_window
end_time = trig_time + args.plot_window

logging.info('Fetching omicron triggers')
omicron_trigs = SnglBurstTable.fetch(hoft_chan,'omicron',start_time,end_time)
omicron_times = np.add(omicron_trigs.get_column('peak_time'),omicron_trigs.get_column('peak_time_ns')*10**-9)
omicron_snr = omicron_trigs.get_column('snr')
omicron_freq = omicron_trigs.get_column('peak_frequency')

window_idx = np.where((omicron_times > start_time) & (omicron_times < end_time))[0]
plot_snr = omicron_snr[window_idx]
plot_freq = omicron_freq[window_idx]
plot_times = omicron_times[window_idx]

logging.info('Loudest Omicron trig has SNR ' + str(max(omicron_snr)))


from pycbc.workflow.segment import fromsegmentxml
import pycbc.pnutils

f_low = 30
f_high = pycbc.pnutils.f_SchwarzISCO(mass1+mass2)

inspiral_t, inspiral_f = pycbc.pnutils.get_inspiral_tf(trig_time, mass1, mass2, f_low, f_high)

logging.info('Plotting')

cm = plt.cm.get_cmap('Reds')
plt.scatter(plot_times,plot_freq,c=plot_snr,s=30,cmap=cm,linewidth=0)
plt.grid(b=True, which='both')
cbar = plt.colorbar()
cbar.set_label('Omicron trigger SNR')
plt.yscale('log')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.xlim(start_time,end_time)
if args.ranking_statistic == 'snr':
    plt.suptitle('CBC trigger SNR = ' + str(snr[highest_idx]) + 
                    ", newSNR = " + str(newsnr[highest_idx]),fontsize=12)
elif args.ranking_statistic == 'newsnr':
    plt.suptitle('CBC trigger SNR = ' + str(snr[highest_idx]) + 
                    ", newSNR = " +str(newsnr[highest_idx]),fontsize=12)
plt.title(str(mass1) + " - " + str(mass2) + " solar masses at GPS time " + str(trig_time),fontsize=12)
plt.hold(True)
plt.plot(inspiral_t,inspiral_f)
plt.savefig(args.output_file)


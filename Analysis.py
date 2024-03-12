
MSDLen = 1

RMSDs = (i * RMSDs + np.array(RMSD_unif)) / (i + 1)
MSD = (i * MSD + np.array(GetMSDR(x_unif, y_unif, MSDLen))) / (i + 1)

RMSDs = np.power(np.array(RMSDs), 2)

fig, ax1 = plt.subplots()
fig2, ax1_b = plt.subplots()
fig3, ax1_c = plt.subplots()

ax2 = ax1.twinx()
ax2_c = ax1_c.twinx()

ax1_c.plot(timescale[MSDLen:],MSD)
ax2_c.plot(timescale,peptide_remaining,c = 'r')

ax1.set_yscale('log',base=10)
ax1.set_xscale('log',base=10)
ax1.plot((timescale),(distances))
ax1.plot((timescale),(timescale))
ax2.plot((timescale),peptide_remaining,c = 'r')
ax2.set_ylim(0,600)





ax1_b.plot(timescale,(distances))
ax2_b = ax1_b.twinx()
ax2_b.plot(timescale,peptide_remaining,c = 'r')

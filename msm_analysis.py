import matplotlib.pyplot as plt
import numpy as np
import pyemma
import mdtraj as md

s = "1"
gro = "res" + s + "_1us_no_solv.gro"
top = md.load(gro).topology
traj = md.load("res" + s + "_1us_no_solv.xtc", top=top)

# The phi[0], psi[1] and chi_1[2] indicies are given in tuples for each of the residues in the dictionary
# Keys are given as str of the structure number or as "tyr" for canonical Tyr. E.g. "1"
res_inds = {"1" : ((5,7,9,31),(7,9,31,33),(7,9,11,14)), "2" : ((5,7,9,26),(7,9,26,33),(7,9,11,14)), "3" : ((5,7,9,26),(7,9,26,35),(7,9,11,14)), "4" : ((5,7,9,27),(7,9,27,41),(7,9,11,14)), "5" : ((5,7,9,26),(7,9,26,32),(7,9,11,14)), "tyr": ((5,7,9,26),(7,9,26,28),(7,9,11,14))}

# Correct to 0-indexed indicies
indicies = []
for val in res_inds[s]:
    dihedral = (val - np.ones(4)).astype(int)
    indicies.append(dihedral)
phi_i, psi_i, chi_i = indicies

print([traj.top.atom(x) for x in chi_i]) 
print([traj.top.atom(x) for x in phi_i]) 
print([traj.top.atom(x) for x in psi_i]) 

dihedrals = md.compute_dihedrals(traj, [phi_i, psi_i])
np.shape(dihedrals)

phis = np.rad2deg(dihedrals[:,0]) + 180
psis = np.rad2deg(dihedrals[:,1]) + 180

state_trjs = [np.array([np.floor(phis[i]/10).astype(int) + 36*np.floor(psis[i]/10).astype(int) for i in range(len(psis))])]

its = pyemma.msm.its(state_trjs, lags=[1, 2, 5, 10, 20, 30, 40, 50], nits=5)

pyemma.plots.plot_implied_timescales(its, units='ps')
dtrajs_concatenated = np.concatenate(state_trjs)

msm = pyemma.msm.estimate_markov_model(state_trjs, lag=10, dt_traj="1 ps")

print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))
print(len(msm.active_set))

phi = np.array([np.floor(phis[i]/10).astype(int) for i in range(len(psis))])
psi = np.array([np.floor(psis[i]/10).astype(int) for i in range(len(psis))])

print("timescales" + str(msm.timescales()[:5]))

eigvec = msm.eigenvectors_left()
print('first eigenvector is one: {} (min={}, max={})'.format(
    np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))

font = {'color': 'black', 'weight': 'semibold', 'size': 20}
    
# Get the active states
active = msm.active_set

fig, axes = plt.subplots(1, 2, figsize=(16, 8), constrained_layout=True)
for i, ax in enumerate(axes.flat):
    
    # Initialize the eigenvector
    left_ev = np.zeros((36,36))
    
    # Assign left eigenvector of state to plot
    for j, m in enumerate(active):
        left_ev[m//36, m%36] = eigvec[i,j]
    
    left_ev = left_ev * 1000
    
    im = ax.pcolormesh(left_ev, cmap=plt.get_cmap('coolwarm'), vmin=-0.03, vmax=0.03)
    ax.set_xlabel(r'$\Phi$ Dihedral [discretized]', fontdict=font, labelpad=10)    
    ax.tick_params(axis='y', labelsize=20, direction='in', width=3, \
                    length=5, pad=10)
    ax.tick_params(axis='x', labelsize=20, direction='in', width=3, \
                    length=5, pad=10)
    for i in ["top","bottom","left","right"]:
        ax.spines[i].set_linewidth(3) 

#fig.subplots_adjust(right=0.85)
axes[0].set_ylabel(r'$\Psi$ Dihedral [discretized]', fontdict=font, labelpad=10)
cbar = plt.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=16, direction='out', width=3, length=5)
cbar.ax.set_yticklabels(np.round(np.arange(-0.03,0.04,0.01),2),fontsize=20)
cbar.outline.set_linewidth(3)
if which_res=="tyr":
    plt.suptitle("Canonical Tyr", fontsize=20)
    plt.savefig("left_ev_" + which_res + ".png")
else:
    plt.suptitle("Fluorinated pTyr " + which_res, fontsize=20) #fontproperties=font)
    plt.savefig("left_ev_" + which_res + ".png")



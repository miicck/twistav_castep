import matplotlib.pyplot as plt
import parse_castep as pc
import numpy as np
import sys
import os

from ast import literal_eval
from scipy.optimize import curve_fit

def out(msg, tag):
	active_tags = []
	active_tags.append("twist_av_energies")
	#active_tags.append("raw_energies")
	active_tags.append("errors")
	if tag not in active_tags:
		return;
	print msg;

def bash(cmd):
	return os.popen(cmd).read()

def get_reblocked_energy(direc):
	for rbl in [-1, 256, 128, 64, 32, 16, 8, 4, 2, 1]:
		cmd = 'printf "%d\n5\n%d\n%d\n" -1 '+str(rbl)+' -1 | reblock'
		res = bash("cd "+direc+"; "+cmd)
		for line in res.split("\n"):
			if "Total energy (using Ewald)" in line:
				ret = [float(w) for w in line.split(":")[-1].split()]
				if rbl != -1:
					return [ret, False]
				else:
					return [ret, True]
	return [None, None]

def twist_av_model(ddft, e_dmc, c):
	return e_dmc + c * ddft

def fit_twist_av(qmc_energies, dft_energies, dft_fine_energy, qmc_errors):
	if dft_fine_energy == None:
		dft_fine_energy = np.mean(dft_energies)
	ddft = np.array(dft_energies) - dft_fine_energy

	f = open("E_v_twist.dat","w")
	for i in range(0,len(qmc_energies)):
		f.write(str(i+1)+","+str(qmc_energies[i])+","+str(qmc_errors[i])+str(",")+str(ddft[i])+"\n")
	f.close();
	qmc_e = None
	for line in bash("tahelper").split("\n"):
		if "Final twist-averaged energy" in line:
			qmc_e = float(line.split(":")[-1])	
			continue
		if "Error in twist-averaged energy" in line:
			dqmc_e = float(line.split(":")[-1])
	bash("rm E_v_twist.dat")

	return [qmc_e, dqmc_e]

def get_tau(direc):
	return float(bash("cd "+direc+"; cat input | grep dtdmc | cut -d: -f2"))

def tau_extrapolated_energy(qmc_e_tau):
	e1  = qmc_e_tau[0][0][0]
	de1 = qmc_e_tau[0][0][1]
	t1  = qmc_e_tau[0][1]
	e2  = qmc_e_tau[1][0][0]
	de2 = qmc_e_tau[1][0][1]
	t2  = qmc_e_tau[1][1]
	de1 = de1 * (1+t1/(t2-t1))
	de2 = -de2 * t1/(t2-t1)
	return [e1 - t1*(e2-e1)/(t2-t1), np.sqrt(de1**2+de2**2)]

def get_cell_data():

	cell_data = []
	cell_ns = []
	for d in os.listdir("output"):
		if not os.path.isdir("output/"+d):
			continue
		if not d.endswith("cell"):
			continue
		n = int(d.split("_")[0])
		cell_ns.append(n)

	cell_ns.sort()

	for n_cells in cell_ns:
		cell_dir = "output/"+str(n_cells)+"_cell"
		out(str(n_cells)+" cells:", "twist_av_energies")

		dft_twist_es  = []
		qmc_twist_es  = []
		qmc_twist_des = []
		for twist_dir in os.listdir(cell_dir):
			if not twist_dir.startswith("twist"):
				continue
			twist_dir = cell_dir + "/" + twist_dir
			out("","raw_energies")
			out(twist_dir,"raw_energies")

			dft_dir = twist_dir + "/dft"
			cast = None
			for f in os.listdir(dft_dir):
				if f.endswith(".castep"):
					cast = pc.parse_castep(dft_dir+"/"+f)
					break
			dft_e = cast["0K energies"][-1]/cast["ions in cell"]
			out("DFT energy: "+str(dft_e)+" eV","raw_energies")

			qmc_dir = twist_dir + "/qmc"
			qmc_e_tau = []
			for tau_dir in os.listdir(qmc_dir):
				tau_dir = qmc_dir + "/" + tau_dir
				qmc_e, reblock_success = get_reblocked_energy(tau_dir)
				if qmc_e == None:
					out ("Error, could not get reblocked energy!","errors")
					continue
				qmc_t = get_tau(tau_dir)
				out("QMC energy (tau = "+str(qmc_t)+"): "+str(qmc_e[0])+" +/- "+str(qmc_e[1])+" eV", "raw_energies")
				if not reblock_success:
					out("(Reblocking failed to identify optimum block length)", "raw_energies")
				qmc_e_tau.append([qmc_e,qmc_t])
			if len(qmc_e_tau) < 2:
				out("Too few tau points to carry out zero-timestep extrapolation!", "errors")
				continue

			qmc_e_tau_0, dqmc_e_tau_0 = tau_extrapolated_energy(qmc_e_tau)
			out("DMC energy (tau = 0): "+str(qmc_e_tau_0)+" +/- "+str(dqmc_e_tau_0), "raw_energies")
			dft_twist_es.append(dft_e)
			qmc_twist_es.append(qmc_e_tau_0)
			qmc_twist_des.append(dqmc_e_tau_0)

		if len(qmc_twist_es) == 0:
			out("Too few twists to carry out infinite-size extrapolation!", "errors")
			continue
		twist_av = fit_twist_av(qmc_twist_es, dft_twist_es, None, qmc_twist_des)

		out("", "raw_energies")
		out("Twist av DFT energy: "+str(np.mean(dft_twist_es))+r" \sigma = "+str(np.std(dft_twist_es)), "twist_av_energies")
		out("Twist av QMC energy: "+str(twist_av[0])+" +/- "+str(twist_av[1]), "twist_av_energies")
		out("", "twist_av_energies")

		dft_mean_diff = list(np.abs( np.array(dft_twist_es) - np.mean(dft_twist_es) ))
		dft_best_index = dft_mean_diff.index(np.min(dft_mean_diff))

		cell_data.append([n_cells, twist_av[0], twist_av[1], np.mean(qmc_twist_es), qmc_twist_es[dft_best_index]])

	cell_data.sort()
	f = open("output/analysis_results","w")
	f.write(str(cell_data))
	return cell_data

def inf_size_func(ninv, e_inf, b):
	return e_inf - b*ninv	

def inf_size_extrapolation(ns, qmc_es, qmc_des):
	ninv = 1.0/ns
	return curve_fit(inf_size_func, ninv, qmc_es)

# Program start
	
min_n = 0           # The minimum supercell size (in prim. cells) to analyse
recalculate = False # Do we re-analyse, or attempt to read the results of previous analysis

# Parse arguments
for a in sys.argv:
	if a.startswith("start="):
		min_n = int(a.split("=")[-1])
	if a == "recalculate":
		recalculate = True

if recalculate:
	# Recalculate analysis results from scratch
	cell_data = get_cell_data()
else:
	# Read in analysis results
	cell_data = literal_eval(open("output/analysis_results").read())

cell_data_new = []
for cd in cell_data:
	if cd[0] < min_n:
		continue
	cell_data_new.append(cd)
cell_data = cell_data_new

# Plot things
ns, qmc_es, qmc_des, qmc_es_simple, qmc_es_dft_best = np.array(zip(*cell_data))
inf_size_par, inf_size_covar = inf_size_extrapolation(ns, qmc_es, qmc_des)
print inf_size_par

plt.plot(1.0/ns, qmc_es, label="QMC twist-fit")
plt.fill_between(1.0/ns, qmc_es + qmc_des/2, qmc_es - qmc_des/2, alpha=0.2)
plt.plot(1.0/ns, qmc_es_simple, label="QMC twist-average")
plt.plot(1.0/ns, qmc_es_dft_best, label="QMC best-dft")
ngrid = np.linspace(0, np.max(1.0/ns), 50)
ex_lab = "Infinie size extrapolation (E_inf = "+str(inf_size_par[0]) 
ex_lab +=" +/- " +str(np.sqrt(inf_size_par[1]))+")"
plt.plot(ngrid, inf_size_func(ngrid, *inf_size_par), label=ex_lab)
plt.xlabel("1/N where N = number of primative cells in supercell")
plt.ylabel("DMC energy (eV/atom)")
plt.legend()
plt.show()

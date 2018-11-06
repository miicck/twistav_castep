import matplotlib.pyplot as plt
import parse_castep as pc
import numpy as np
import sys
import os

from ast import literal_eval
from scipy.optimize import curve_fit

# Utility for enabling/disabling various sources of output
def out(msg, tag):
	active_tags = []
	active_tags.append("twist_av_energies")
	#active_tags.append("raw_energies")
	active_tags.append("errors")
	if tag not in active_tags:
		return;
	print msg;

# Run cmd in the shell and return the result
def bash(cmd):
	return os.popen(cmd).read()

# Get the reblocked energy of a CASINO calulation in
# the given directory, uses the "reblock" utility.
# Returns [[qmc_e, qmc_de], reblock]
# where qmc_e   = energy
#       qmc_de  = error in energy
#	reblock = true if reblocking was successful
def get_reblocked_energy(direc):
	# Look first for preferred block length, then longest
	for rbl in [-1, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]:
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

# The twist averaging fitting function
def twist_av_model(ddft, e_dmc, c):
	return e_dmc + c * ddft

# Carries out twist averaging of CASINO calculations, using the "tahelper" utility.
def fit_twist_av(qmc_energies, dft_energies, dft_fine_energy, qmc_errors):

	# If no fine DFT energy was provided, use the mean
	# of the dft_energies given for particular twists.
	if dft_fine_energy == None:
		dft_fine_energy = np.mean(dft_energies)

	# Calculate the difference of each dft energy
	# from the fine DFT energy
	ddft = np.array(dft_energies) - dft_fine_energy

	# Write a file that "tahelper" can read.
	f = open("E_v_twist.dat","w")
	for i in range(0,len(qmc_energies)):
		f.write(str(i+1)+","+str(qmc_energies[i])+","+str(qmc_errors[i])+str(",")+str(ddft[i])+"\n")
	f.close();

	# Run "tahelper" and parse the result for the twist
	# averaged energy and associated error.
	qmc_e = None
	for line in bash("tahelper").split("\n"):
		if "Final twist-averaged energy" in line:
			qmc_e = float(line.split(":")[-1])	
			continue
		if "Error in twist-averaged energy" in line:
			dqmc_e = float(line.split(":")[-1])
	bash("rm E_v_twist.dat")

	return [qmc_e, dqmc_e]

# Get the DMC timestep used in a CASINO calculation in the given directory
# (assumes the input file is still present and correct)
def get_tau(direc):
	return float(bash("cd "+direc+"; cat input | grep dtdmc | cut -d: -f2"))

# Extrapolate for the zero-timestep qmc_energy
# from data in the form [[[qmc_e, qmc_de], tau], ...]
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

# Traverse the output directory and gather data
def get_cell_data():

	cell_data = []
	cell_ns = []

	# Get the number of cells in each supercell
	# (by parsing the filename)
	for d in os.listdir("output"):
		if not os.path.isdir("output/"+d):
			continue
		if not d.endswith("cell"):
			continue
		n = int(d.split("_")[0])
		cell_ns.append(n)

	# Sort so we work in a nice order
	cell_ns.sort()

	# Run over each supercell
	for n_cells in cell_ns:

		# Construct the directory of this supercell
		cell_dir = "output/"+str(n_cells)+"_cell"
		out(str(n_cells)+" cells:", "twist_av_energies")

		dft_twist_es  = []
		qmc_twist_es  = []
		qmc_twist_des = []

		# Run over each twist for this supercell
		for twist_dir in os.listdir(cell_dir):

			# Check that this file is a twist directory
			if not twist_dir.startswith("twist"):
				continue
			twist_dir = cell_dir + "/" + twist_dir
			if not os.path.isdir(twist_dir):
				continue

			out("","raw_energies")
			out(twist_dir,"raw_energies")

			# Open and parse the .castep file
			# for the DFT calulation of this twist
			dft_dir = twist_dir + "/dft"
			cast = None
			for f in os.listdir(dft_dir):
				if f.endswith(".castep"):
					cast = pc.parse_castep(dft_dir+"/"+f)
					break

			# Get the DFT energy of this twist
			dft_e = cast["0K energies"][-1]/cast["ions in cell"]
			out("DFT energy: "+str(dft_e)+" eV","raw_energies")

			# Run over the different tau values for which a
			# DMC calculation was ran for this twist
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

			# Need at least two points to carry out zero-timestep extrapolation
			if len(qmc_e_tau) < 2:
				out("Too few tau points to carry out zero-timestep extrapolation!", "errors")
				continue

			# Extrapolate this twist to zero timestep
			qmc_e_tau_0, dqmc_e_tau_0 = tau_extrapolated_energy(qmc_e_tau)
			out("DMC energy (tau = 0): "+str(qmc_e_tau_0)+" +/- "+str(dqmc_e_tau_0), "raw_energies")

			# Record results for this twist
			dft_twist_es.append(dft_e)
			qmc_twist_es.append(qmc_e_tau_0)
			qmc_twist_des.append(dqmc_e_tau_0)

		
		# Carry out twist averaging
		if len(qmc_twist_es) == 0:
			out("Too few twists to carry out twist averaging!", "errors")
			continue
		twist_av = fit_twist_av(qmc_twist_es, dft_twist_es, None, qmc_twist_des)

		# Output some things
		out("", "raw_energies")
		out("Twist av DFT energy: "+str(np.mean(dft_twist_es))+r" \sigma = "+str(np.std(dft_twist_es)), "twist_av_energies")
		out("Twist av QMC energy: "+str(twist_av[0])+" +/- "+str(twist_av[1]), "twist_av_energies")
		out("", "twist_av_energies")

		# Accumulate the resulting E_{twist-av DMC}(N)
		cell_data.append([n_cells, twist_av[0], twist_av[1]])

	# Write resulting cell_data to analysis_results file
	# First sort by n
	cell_data.sort()
	f = open("output/analysis_results","w")
	f.write(str(cell_data))
	return cell_data

# Function to fit infinite size extrapolation to
def inf_size_func(ninv, e_inf, b):
	return e_inf - b*ninv	

# Carry out infinite size extrapolation
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
	try:
		cell_data = literal_eval(open("output/analysis_results").read())
	except:
		cell_data = get_cell_data()

# Discard supercells with < min_n primitive cells
cell_data_new = []
for cd in cell_data:
	if cd[0] < min_n:
		continue
	cell_data_new.append(cd)
cell_data = cell_data_new

# Plot things
ns, qmc_es, qmc_des = np.array(zip(*cell_data))
par_inf, cov_inf    = inf_size_extrapolation(ns, qmc_es, qmc_des)

# PLot infinite size extrapolation
plt.plot(1.0/ns, qmc_es, label="QMC twist-averaged, zero-timestep energy")
plt.fill_between(1.0/ns, qmc_es + qmc_des/2, qmc_es - qmc_des/2, alpha=0.2)
ngrid = np.linspace(0, np.max(1.0/ns), 5)
ex_lab = "Infinie size extrapolation (E_inf = "+str(par_inf[0]) 
ex_lab +=" +/- " +str(np.sqrt(cov_inf[0][0]))+")"
plt.plot(ngrid, inf_size_func(ngrid, *par_inf), label=ex_lab)
plt.xlabel("1/N where N = number of primative cells in supercell")
plt.ylabel("DMC energy (eV/atom)")
plt.legend()
plt.show()

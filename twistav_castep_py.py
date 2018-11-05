# A script which carries out twist-averaged zero-time-step-extrapolated
# infinite size extrapolation of crystalline systems in qmc

import os
import numpy as np
import parse_castep as pc
import sys

# Get the number of electrons in the system
def get_electrons():
	if len(sys.argv) < 2:
		print "Error: I require one argument, the number of electrons in the primitive cell!"
		quit()
	return int(sys.argv[1])

# Run cmd in a bash terminal and
# return the output
def bash(cmd):
	return os.popen(cmd).read()

# Print an emphasied line
def print_header(header):
	line = ""
	for i in range(0,len(header)):
		line += "-"
	print ""
	print line
	print header
	print line

# Get the optimium supercell and kpoints for
# a particular lattice and number of cells
def get_supercell(sc_lat, ncell):
	
	global supercell_cmd
	print "Calling supercell utility..."
	sc_cmd = 'printf "'+sc_lat + str(ncell) + "\nG\n\nn\n" + '"' + " | " + supercell_cmd
	sc_res = bash(sc_cmd)

	sc_lines = sc_res.split("\n")
	sc_matrix_text = None

	kv_lines = []
	kvs = []
	for i, scl in enumerate(sc_lines):
		if "Best geometry with" in scl:
			sc_matrix_text = sc_lines[i+3] + "\n" + sc_lines[i+4] + "\n" + sc_lines[i+5]
		if "For k-vector offset" in scl:
			kv_lines.append(i)
			kvs.append([scl.split(":")[-1],[]])
	kv_lines.append(len(sc_lines))

	for i in range(0,len(kv_lines)-1):
		for j in range(kv_lines[i]+2, kv_lines[i+1]-1):
			kvs[i][1].append(sc_lines[j])

	return [sc_matrix_text, kvs]

# Make the dft input .cell file with the given kpoints and
# (directory quilified) filename, w.r.t the template .cell file
def make_dft_input_cell(kpoints, filename, template_cellfile):
	to_write =  "# ===== Written by twistav_castep_py.py ===== "
	to_write += "\n%block kpoints_list"
	for k in kpoints:
		to_write += "\n"+k
	to_write += "\n%endblock kpoints_list"
	to_write += "\n# =========================================== "
	f = open(filename, "w")
	f.write(to_write + "\n\n" + open(template_cellfile).read())

# Make the qmc input for a system with the given supercell matrix, number of cells
# electrons per cell and (directory-qualified) filename
def make_qmc_input(supercell_matrix_text, ncells, nue, nde, tau, filename):
	to_write =  "# ====== Written by twistav_castep_py.py ===== "
	to_write += "\ndtdmc : " + str(tau)
	to_write += "\nneu   : " + str(nue)
	to_write += "\nned   : " + str(nde)
	to_write += "\nperiodic : T"
	to_write += "\natom_basis_type : blip"
	to_write += "\n%block scell_matrix\n"
	to_write += supercell_matrix_text
	to_write += "\n%endblock scell_matrix\n"
	to_write += "# ============================================ "
	f = open(filename,"w")
	f.write(to_write+"\n\n"+open("input").read())

def get_reblocked_energy(direc):
	cmd = 'printf "%d\n5\n%d\n" -1 -1 | reblock'
	res = bash("cd "+direc+"; "+cmd)
	for line in res.split("\n"):
		if "Total energy (using Ewald)" in line:
			return [float(w) for w in line.split(":")[-1].split()]


# ======================
# ==== MAIN PROGRAM ====
# ======================

electrons     = get_electrons()
max_cells     = 50
tau           = 0.03
supercell_cmd = "supercell"
castep_cmd    = "nice -15 mpirun castep.mpi"
castep2casino = "castep2casino"
casino_cmd    = "nice -15 runqmc"

# Read the input cell file
cell_file = None
template_cell_file = None
system_name = None
for f in os.listdir("."):
	if f.endswith(".cell"):
		cell_file = pc.parse_cell(f)
		template_cell_file = f
		system_name = f.split("/")[-1].split(".")[0]

# Construct the lattice vector string to pass to supercell utility
lat = [cell_file["lattice a"], cell_file["lattice b"], cell_file["lattice c"]]
sc_lat = ""
for l in lat:
	sc_lat += str(l[0]) + " " + str(l[1]) + " " + str(l[2]) + "\n"	

# Run over each cell size
for ncell in range(1, max_cells+1):
	
	print_header("Building supercell from "+str(ncell)+" primative cell(s)")

	# Create the directory for this number of cells, skip if exists
	cell_dir = "output/"+str(ncell)+"_cell_in_progress"
	cell_dir_complete = "output/"+str(ncell)+"_cell"
	if os.path.exists(cell_dir) or os.path.exists(cell_dir_complete):
		print cell_dir + " exists! skipping..."
		continue
	os.makedirs(cell_dir)

	# Work out number of electrons in cell (minimize spin polarization)
	nue = ncell * electrons / 2
	nde = ncell * electrons - nue
	print "Up electrons   : "+str(nue)
	print "Down electrons : "+str(nde)

	# Get optimum supercell and kpoint-sets (twists) that give
	# real-valued wavefunctions
	sc_matrix, kvs = get_supercell(sc_lat, ncell)
	
	print "Optimal supercell:"
	print sc_matrix
	print "Twists resulting in real wvfns : "+str(len(kvs))
	print "Kpoints per twist              : "+str(len(kvs[0][1]))


	# Run over all twists 
	for twist_i in range(0,len(kvs)):

		# Print header
		twist_header="Running twist "+str(twist_i+1)+"/"+len(kvs)+ " ("+str(len(kvs[twist_i][1]))+" kpoint(s))"
		twist_sep = ""
		for tsi in range(0, len(twist_header)):
			twist_sep += "-"
		print twist_sep
		print twist_header

		# Create twist directory
		twist_dir = cell_dir+"/twist_"+str(twist_i)

		# Create dft directory/copy all files from base dir into it
		dft_dir = twist_dir + "/dft"
		os.makedirs(dft_dir)
		bash("cp * "+dft_dir+" >/dev/null 2>&1")

		# Make the input .cell file	
		make_dft_input_cell(kvs[twist_i][1],dft_dir+"/"+system_name+".cell",template_cell_file)

		# Run CASTEP to generate orbitals
		print "DFT..."
		cast_stout = bash("cd "+dft_dir+";"+castep_cmd+" "+system_name)
		if len(cast_stout.strip()) > 0:
			print cast_stout
		cast = pc.parse_castep(dft_dir+"/"+system_name+".castep")
		print "DFT energy: "+str(cast["0K energies"][-1]/cast["ions in cell"])+" eV/atom"

		# Convert orbitals into CASINO plane-wave format, and
		# then into a blip basis
		bash("cd "+dft_dir+";"+castep2casino+" "+system_name)
		bash("cd "+dft_dir+"; mv "+system_name+".casino pwfn.data")
		bash("cd "+dft_dir+'; printf "1\n0\nn\nn\n" | blip')

		# Run over time-steps (tau) to allow extrapolation to zero timestep
		for tau_mult in [1,4]:

			# Create qmc directory/copy all files from base dir into it
			qmc_dir = twist_dir + "/qmc/tau_" + str(tau_mult)
			os.makedirs(qmc_dir)
			bash("cp * "+dft_dir+" "+qmc_dir+" >/dev/null 2>&1")
			bash("cp "+dft_dir+"/bwfn.data "+qmc_dir) # Copy blip orbitals to qmc directory

			# Make the qmc input file
			make_qmc_input(sc_matrix, ncell, nue, nde, tau*float(tau_mult), qmc_dir+"/input")

			# Run CASINO
			print "QMC (tau = tau_min x " + str(tau_mult) + " = "+str(tau*float(tau_mult))+") ..."
			qmc_stout = bash("cd "+qmc_dir+"; "+casino_cmd)
			if len(qmc_stout.strip()) > 0:
				print qmc_stout
			try:
				qmc_energy, qmc_error = get_reblocked_energy(qmc_dir)
				print "QMC energy: "+str(qmc_energy)+" +/- "+str(qmc_error)+" ev/atom"
			except:
				print "QMC reblocking failed, best estimate of energy unknown."

	print cell_dir + " complete, renaming to " + cell_dir_complete
	bash("mv "+cell_dir+" "+cell_dir_complete)

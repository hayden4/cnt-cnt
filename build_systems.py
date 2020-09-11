import math
import sys
import os.path
import shutil

C_MASS = 12.011000
C_C_PAIRWISE = [0.06400, 4.0100]
C_C_BONDLENGTH = 1.417000
C_C_BOND = [C_C_BONDLENGTH, 470.836100, -627.617900, 1327.634500]


# Creates an armchair (n,n) CNT of specified dimensions
class CarbonNanotube:
	
	# s - C-C bond separation
	# n - horizontal vector size ( the 'n' from (n,n) )
	# L - length of CNT
	def __init__(self, s, n, L):
		self.s = float(s)
		self.num_atoms = n+1
		self.L = float(L)
		self.atoms = []
		self.bonds = []
	
	# Constructs the atoms and bonds in the tube
	def buildCNT(self):
		dt = (2 * math.pi) / float(self.num_atoms)
		radius = self.s * math.sqrt(3 / (2 * (1 - math.cos(dt))))
		for i in range(self.num_atoms):
			self.atoms.append([radius * math.cos(i * dt), radius * math.sin(i * dt), 0])
		z = 0
		while z < self.L:
			# sublayer 1
			z += self.s / 2
			if z > self.L: break
			for i in range(self.num_atoms-1):
				self.atoms.append([radius * math.cos(i*dt + (dt/2)), radius * math.sin(i*dt + (dt/2)), z])
				self.bonds.append([len(self.atoms) - self.num_atoms    , len(self.atoms)])
				self.bonds.append([len(self.atoms) - self.num_atoms + 1, len(self.atoms)])
			self.atoms.append([radius * math.cos(-dt/2), radius * math.sin(-dt/2), z])
			self.bonds.append([len(self.atoms) - 2*self.num_atoms+1, len(self.atoms)])
			self.bonds.append([len(self.atoms) -   self.num_atoms  , len(self.atoms)])

			# sublayer 2
			z += self.s
			if z > self.L: break
			for i in range(self.num_atoms):
				self.atoms.append([radius * math.cos(i*dt + (dt/2)), radius * math.sin(i*dt + (dt/2)), z])
				self.bonds.append([len(self.atoms) - self.num_atoms, len(self.atoms)])

			# sublayer 3
			z += self.s / 2
			if z > self.L: break
			self.atoms.append([radius, 0, z])
			self.bonds.append([len(self.atoms)-1             , len(self.atoms)])
			self.bonds.append([len(self.atoms)-self.num_atoms, len(self.atoms)])
			for i in range(1, self.num_atoms):
				self.atoms.append([radius * math.cos(i*dt), radius * math.sin(i*dt), z])
				self.bonds.append([len(self.atoms)-self.num_atoms-1, len(self.atoms)])
				self.bonds.append([len(self.atoms)-self.num_atoms  , len(self.atoms)])

			# sublayer 4
			z += self.s
			if z > self.L: break
			for i in range(self.num_atoms):
				self.atoms.append([radius * math.cos(i*dt), radius * math.sin(i*dt), z])
				self.bonds.append([len(self.atoms)-self.num_atoms, len(self.atoms)])

		# Translate so center of tube is at the origin
		com = [0,0,0]
		for i in range(3):
			for pos in self.atoms:
				com[i] += pos[i]
			com[i] /= len(self.atoms)

		for j in range(3): 
			for i in range(len(self.atoms)):
				self.atoms[i][j] -= com[j]
	

	# Translates the tube a distance r in the y-direction
	def translate(self, r):
		for i in range(len(self.atoms)):
			self.atoms[i][1] += r
	

	# Rotates the CNT so the central axis matches the spherical unit
	# vector (1, theta, phi)
	# Angles expected to be in radians
	def rotate(self, theta, phi):
		for i in range(len(self.atoms)):
			pos = self.atoms[i]
			self.atoms[i] = [
				pos[0] * math.cos(theta) + pos[2] * math.sin(theta),
				pos[1],
				-pos[0] * math.sin(theta) + pos[2] * math.cos(theta)
			]

		for i in range(len(self.atoms)):
			pos = self.atoms[i]
			self.atoms[i] = [
				pos[0] * math.cos(phi) - pos[1] * math.sin(phi),
				pos[0] * math.sin(phi) + pos[1] * math.cos(phi),
				pos[2]
			]
		



def cnt_system(theta, phi):
	# CNT Configuration
	# s = C-C bond length
	# n = number of carbon atoms in ring - 1
	# L = length of CNT
	n = 10
	L = 150
	s = C_C_BONDLENGTH

	# Build CNT 1
	cnt1 = CarbonNanotube(s, n, L)
	cnt1.buildCNT()

	# Build CNT 2
	cnt2 = CarbonNanotube(s, n, L)
	cnt2.buildCNT()

	# Rotate CNT 2
	cnt2.rotate(theta, phi)

	return (cnt1, cnt2)


def write_system(filename, cnt1, cnt2):
	# Combine tubes
	atoms = [a for a in cnt1.atoms]
	bonds = [b for b in cnt1.bonds]

	atomid_offset = len(atoms)
	for i in range(len(cnt2.atoms)):
		atoms.append(cnt2.atoms[i])
	
	for i in range(len(cnt2.bonds)):
		bonds.append([cnt2.bonds[i][0]+atomid_offset, cnt2.bonds[i][1]+atomid_offset])
	
	# Write file
	outfile = open(filename, 'w')
	outfile.write("""LAMMPS Description

	{a} atoms
	{b} bonds
	1 atom types
	1 bond types
	{xl} {xh} xlo xhi
	{yl} {yh} ylo yhi
	{zl} {zh} zlo zhi

Masses

1 {m}

Atoms

""".format(a=len(atoms), b=len(bonds), xl=-100, xh=100, yl=-100, yh=100, zl=-100, zh=100, m=C_MASS))
	
	for i in range(len(atoms)):
		outfile.write("{i} {m} 1 {q} {x} {y} {z}\n".format(i=i+1, m=1 if i < atomid_offset else 2, q=0, x=atoms[i][0], y=atoms[i][1], z=atoms[i][2]))
	
	outfile.write("\nBonds\n\n")

	for i in range(len(bonds)):
		outfile.write("{i} 1 {a} {b}\n".format(i=i+1, a=bonds[i][0], b=bonds[i][1]))
	
	outfile.write("""
Pair Coeffs

1 {a} {b}

Bond Coeffs

1 {c} {d} {e} {f}
""".format(a=C_C_PAIRWISE[0], b=C_C_PAIRWISE[1], c=C_C_BOND[0], d=C_C_BOND[1], e=C_C_BOND[2], f=C_C_BOND[3]))
	
	
	outfile.close()


def main():
	if len(sys.argv) != 3:
		print("Requires argument: theta, phi in degrees")
	
	theta = math.pi * float(sys.argv[1]) / 180.0
	phi = math.pi * float(sys.argv[2]) / 180.0

	print("Theta: {a}".format(a=theta))
	print("Phi: {b}".format(b=phi))

	offset = 25.0 * math.sin(theta) * math.sin(phi) / 2
	r_high = 20.0 + offset
	r_low  = 10.0 + offset
	num_steps = 100

	step_size = (r_high - r_low) / num_steps

	(cnt1, cnt2) = cnt_system(theta, phi)

	if not os.path.exists("sims"):
		os.makedirs("sims")
	
	r = r_low
	cnt2.translate(r)

	for i in range(num_steps):
		r += step_size
		dir_name = "{a:2.3f}".format(a=r)
		dir_path = "sims/{a}".format(a=dir_name)

		print("Generating system at separation {a} in path: {b}".format(a=dir_name, b=dir_path))

		if not os.path.exists(dir_path):
			os.makedirs(dir_path)

		# write run file
		shutil.copy2("run.in", dir_path)
		# write pbs file
		shutil.copy2("potential_test.pbs", dir_path + "/potential_{a}.pbs".format(a=dir_name))
		# translate CNT 2 to correct position
		cnt2.translate(step_size)
		# write data file
		write_system(dir_path + "/cnt.lammps", cnt1, cnt2)

	print("Done generating systems.")


if __name__ == '__main__':
	main()

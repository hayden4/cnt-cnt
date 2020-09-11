from os import listdir
from os.path import isdir, join

def readfile(filename):
	infile = open(filename, 'r')

	NO_STATE = 0
	READ_DATA = 1
	STATE = NO_STATE

	potential = 0.0

	for line in infile:
		if line.startswith("checked_atoms:"):
			STATE = READ_DATA
			continue

		if STATE == READ_DATA:
			items = line.split()
			if (items[4] == '1'):
				potential = float(items[5])
				STATE = NO_STATE
	
	infile.close()

	return potential


def readfixfile(filename):
	infile = open(filename, 'r')

	data = []

	for line in infile:
		if line.startswith("#"):
			continue

		items = line.split()
		data.append(float(items[1]))
	
	infile.close()

	return data[0]


def torque(c, x, f):
	r = [x[i] - c[i] for i in range(3)]
	t = [0 for i in range(3)]
	t[0] = r[1] * f[2] - r[2] * f[1]
	t[1] = r[2] * f[0] - r[0] * f[2]
	t[2] = r[0] * f[1] - r[1] * f[0]
	return t


def readDumpFile(filename):
	print("Reading dump file: " + filename)
	infile = open(filename, 'r')

	tube1 = [[] for i in range(2)]
	tube2 = [[] for i in range(2)]

	NO_STATE = 0
	READ_ATOMS = 1

	state = NO_STATE

	for line in infile:
		line = line.strip()
		if line.startswith("ITEM: ATOMS"):
			state = READ_ATOMS
			continue

		if line.startswith("ITEM: TIMESTEP"):
			if state == READ_ATOMS:
				break
			else:
				state = NO_STATE

		if state == READ_ATOMS:
			items = line.split()
			type = int(items[1])
			x = [float(i) for i in items[2:5]]
			f = [float(i) for i in items[5:8]]

			if type == 1:
				tube1[0].append(x)
				tube1[1].append(f)
			else:
				tube2[0].append(x)
				tube2[1].append(f)
		

	infile.close()

	c1 = [0 for i in range(3)]
	for x in tube1[0]:
		for i in range(3):
			c1[i] += x[i]
	
	c2 = [0 for i in range(3)]
	for x in tube2[0]:
		for i in range(3):
			c2[i] += x[i]

	for i in range(3):
		c1[i] /= len(tube1[0])
		c2[i] /= len(tube2[0])
	
	print("c1: " + str(c1))
	print("c2: " + str(c2))

	print("num atoms tube 1: " + str(len(tube1[1])))
	print("num atoms tube 2: " + str(len(tube2[1])))

	nf1 = [0 for i in range(3)]
	for f in tube1[1]:
		for i in range(3):
			nf1[i] += f[i]

	nf2 = [0 for i in range(3)]
	for f in tube2[1]:
		for i in range(3):
			nf2[i] += f[i]
	
	t1 = [0 for i in range(3)]
	for i in range(len(tube1[0])):
		t = torque(c1, tube1[0][i], tube1[1][i])
		for j in range(3):
			t1[j] += t[j]
	
	t2 = [0 for i in range(3)]
	for i in range(len(tube2[0])):
		t = torque(c2, tube2[0][i], tube2[1][i])
		for j in range(3):
			t2[j] += t[j]
	
	results = []
	for i in range(3):
		results.append(nf1[i])
	for i in range(3):
		results.append(nf2[i])
	for i in range(3):
		results.append(t1[i])
	for i in range(3):
		results.append(t2[i])

	return results



def main():

	dirs = [f for f in listdir("sims") if isdir(join("sims", f))]
	#data = sorted([(float(dir), readfile(join("sims", dir, "post_run.data"))) for dir in dirs], key=lambda x: x[0])
	cnt_data = sorted(
		[(float(dir),
		 readfixfile(join("sims", dir, "cnt_interaction_force.data")),
		 readDumpFile(join("sims", dir, "cnt.dump")))
		 for dir in dirs], key=lambda x: x[0])

	#outfile = open("potential_table.data", 'w')
	#for (r, pot) in data:
	#	outfile.write("{a:2.3f}\t{b}\n".format(a=r, b=pot))
	#	
	#outfile.close()

	outfile = open("cnt_potential_table.data", 'w')
	outfile.write("# sep	energy	f1x	f1y	f1z	f2x	f2y	f2z	t1x	t1y	t1z	t2x	t2y	t2z\n")
	for (r, pot, forces) in cnt_data:
		line = "{a:2.3f}\t{b}".format(a=r, b=pot)
		for x in forces:
			line += "\t{c}".format(c=x)
		line += "\n"
		outfile.write(line)
	
	outfile.close()



if __name__ == '__main__':
	main()

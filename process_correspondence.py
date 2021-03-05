import sys,string
def process_correspondence(correspond_file):
	correspondence = {}
	genomes = {}
	with open(correspond_file, "r") as correspondFile:
		for line in correspondFile:
			if not line:
				continue
			lineSplit = line.strip().split()
			if len(lineSplit)<2:
				continue
			genome = lineSplit[1]
			spectrumfile = lineSplit[0]
			if genome in correspondence:
				correspondence[genome].append(spectrumfile)
			else:
				correspondence[genome] = [spectrumfile]
	return correspondence

	

if __name__ == "__main__":
	correspond_file = sys.argv[1]
	correspondence = process_correspondence(correspond_file)
	with open(correspond_file+".NRPminer_processed","w") as foratted_correspondence:
		for genome in correspondence:
			foratted_correspondence.write(genome + "\t" + "\t".join(correspondence[genome]) + "\n")


from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import argparse
import re
import time


parser = argparse.ArgumentParser(description="Get the sequences you want from a fasta file.")
parser.add_argument("fasta", help="fasta file that has a corresponding .tax file in the same directory")
parser.add_argument("--taxa", help="give the name of taxa that you desire... not case-sensitive, replace spaces with underscores. Takes comma separated lists with no spaces. Put a '~' before the first item to take all taxa EXCEPT the given ones. Put a '+' at the beginning and it will take all descendents in that group. You can do both modes at once.")
parser.add_argument("--out", default="new_fasta.fasta", help="name for output file")
parser.add_argument("--taxfile", help="only necessary if your .tax file is in a different directory than your fasta file")
parser.add_argument("--names", help="extract sequences with a given sequence name")
args = parser.parse_args()
Entrez.email = "jason.vailionis@uconn.edu"

if args.taxfile is None:
	taxfile = "".join(args.fasta.split(".")[:-1]) + ".tax"
else:
	taxfile = args.taxfile

if args.taxa is None:
	no_taxa = True
else:
	no_taxa = False
	taxa_list = re.sub("[+~]", "", args.taxa).split(",")
	taxa_list = [re.sub("_", " ", item).lower() for item in taxa_list]

if args.names is None:
	no_names = True
else:
	no_names = False
	names_list = re.sub("[~]", "", args.names).split(",")
	names_list = [re.sub("_", " ", item).lower() for item in names_list]

if no_taxa and no_names:
	raise Exception("No positional arguments given")

fasta = list(SeqIO.parse(args.fasta, "fasta"))

def taxa_is_true(taxa):
	if no_taxa is True:
		return no_taxa
	taxa = taxa.lower()
	if "+" in args.taxa:
		if taxa == "bacteria":
			taxid = list('2')
		elif taxa == "actinobacteria":
			taxid = list('201174')
		else:
			handle = Entrez.esearch(db="Taxonomy", term=taxa)
			record = handle.read().replace("\n","")
			taxid = re.findall("<Id>([0-9]*)</Id>", record)
			handle.close()
		if len(taxid) > 1:
			print("An Entrez search with the name '" + taxa + "' gave more than one result. Taxonomy cannot be determined so this sequence will be included. I don't know any cases in which this occurs except 'bacteria' and 'actinobacteria', which are handled separately in this program, but who knows..")
			return True
		print(taxid)
		handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
		record = handle.read().replace("\n","")
		lineage = re.findall("<Lineage>([A-Za-z; /-]*)</Lineage>", record)
		handle.close()
		time.sleep(0.4)     # entrez only allows three queries per second
		if len(lineage) is 0:     # this happens when the taxa is N/A or cellular organisms
			return accept(taxa in taxa_list, "~" in args.taxa)
		else:
			print(lineage)
			lineage = [l.lower() for l in lineage[0].split("; ")]
			lineage.append(taxa)
			for i in taxa_list:
				if i in lineage:
					return accept(True, "~" in args.taxa)
			return accept(False, "~" in args.taxa)
	else:
		return accept(taxa in taxa_list, "~" in args.taxa)

def name_is_true(seqname):
	if no_names is True:
		return no_names
	else:
		for name in names_list:
			if name in seqname:
				return accept(True, "~" in args.names)
		return accept(False, "~" in args.names)

def accept(match, exclude):
	if match is True and exclude is False:
		return True
	if match is True and exclude is True:
		return False
	if match is False and exclude is False:
		return False
	if match is False and exclude is True:
		return True

with open(args.out, "w") as o:
	with open("".join(args.out.split(".")[:-1]) + ".tax", "w") as t:
		with open(taxfile, "r") as f:
			tax_parse = enumerate([line.rstrip("\n").split("\t") for line in f])
		for i, info in tax_parse:     # i is the index, info[0] is the seqid, info[1] is the taxonomy
			if taxa_is_true(info[1]) and name_is_true(info[0]):
				o.write(">" + fasta[i].id + "\n" + str(fasta[i].seq) + "\n")
				t.write(info[0].lstrip(">") + "\t" + info[1] + "\n")

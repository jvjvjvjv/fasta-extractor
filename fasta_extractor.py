from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import re

parser = argparse.ArgumentParser(description="Get the sequences you want from a fasta file.")
parser.add_argument("fasta", help="fasta file that has a corresponding .tax file in the same directory")
parser.add_argument("--taxa", help="give the name of taxa that you desire... not case-sensitive, replace spaces with underscores. Takes comma separated lists with no spaces. Put a '~' before the first item to take all taxa EXCEPT the given ones")
parser.add_argument("--out", default="new_fasta.fasta", help="name for output file")
parser.add_argument("--taxfile", help="only necessary if your .tax file is in a different directory than your fasta file")
# make it so you can extract a sequence based on seqname ie give it "I6732" and it will return everything that matches, or even just by index
args = parser.parse_args()

if args.taxa.startswith("~"):
	exclude = True
else:
	exclude = False

if args.taxfile is None:
	taxfile = "".join(args.fasta.split(".")[:-1]) + ".tax"
else:
	taxfile = args.taxfile

taxlist = args.taxa.lstrip("~").split(",")
taxlist = [re.sub("_", " ", item).lower() for item in taxlist]
fasta = list(SeqIO.parse(args.fasta, "fasta"))

with open(args.out, "w") as o:
	with open(taxfile, "r") as f:
		tax_parse = enumerate([line.rstrip("\n").split("\t") for line in f])
		for i, info in tax_parse:
			if (info[1].lower() not in taxlist) == exclude:
				o.write(">" + fasta[i].id + "\n")
				o.write(str(fasta[i].seq) + "\n")

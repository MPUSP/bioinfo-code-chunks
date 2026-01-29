import argparse
import os
from Bio import Blast  # pyright: ignore[reportMissingImports]
from Bio import SeqIO  # pyright: ignore[reportMissingImports]

# get input parameters command line arguments
parser = argparse.ArgumentParser(description="Run BLASTP search")
parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
parser.add_argument("-o", "--output", required=True, help="Output directory")
parser.add_argument("-u", "--email", required=True, help="Email for NCBI")
parser.add_argument(
    "-f", "--filter", required=False, help="Entrez query filter", default=None
)
parser.add_argument(
    "--matrix-name",
    type=str,
    default="BLOSUM45",
    help="Substitution matrix when running BLASTP",
)
parser.add_argument(
    "--gapcosts",
    type=str,
    default="13 3",
    help="Gap opening and extension costs formatted as two integers",
)
parser.add_argument(
    "--expect",
    type=float,
    default=0.001,
    help="E-value cutoff for filtering BLAST hits",
)
parser.add_argument(
    "--hitlist-size",
    type=int,
    default=1000,
    help="Maximum number of BLAST hits to return",
)
args = parser.parse_args()
print(f"Started blastp for input: {args.input}")

# set email for NCBI
Blast.email = args.email

# read the provided FASTA file
fasta_input = SeqIO.read(args.input, "fasta")
print(f"Input sequence has length: {len(fasta_input.seq)}")

# perform BLASTP search against the refseq_protein database optionally using a filter in Entrez format
print("Submitting BLASTP request to NCBI...")
result_stream = Blast.qblast(
    program="blastp",
    database="refseq_protein",
    sequence=str(fasta_input.seq),
    matrix_name=args.matrix_name,
    gapcosts=args.gapcosts,
    expect=args.expect,
    hitlist_size=args.hitlist_size,
    format_type="XML",
    entrez_query=args.filter,
)

# make sure output directory exists
os.makedirs(args.output, exist_ok=True)

# save the full BLAST result to a file
with open(f"{args.output}/results.xml", "wb") as out_stream:
    out_stream.write(result_stream.read())
result_stream.close()
print(f"BLASTP result saved to: {args.output}/results.xml")

# parse the BLAST result
with open(f"{args.output}/results.xml", "rb") as result:
    blast_record = Blast.read(result)
print(f"Total hits found: {len(blast_record)}")

# print the top 3 hits
print("Top 3 BLASTP hits:\n")
for hit in blast_record[0:3]:
    for hsp in hit:
        print("****Alignment****")
        print("protein:", hit.target.name)
        print("description:", hit.target.description)
        print("length:", len(hit.target))
        print("e value:", hsp.annotations["evalue"])
        print("identity:", hsp.annotations["identity"])
        print(hsp)

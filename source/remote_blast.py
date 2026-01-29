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
args = parser.parse_args()
print(f"Started blastp for input: {args.input}")

# set email for NCBI
Blast.email = args.email

# read the provided FASTA file
fasta_input = SeqIO.read(args.input, "fasta")
print(f"Input sequence has length: {len(fasta_input.seq)}")

# perform BLASTP search against the refseq_protein database, excluding genus Streptococcus
print("Submitting BLASTP request to NCBI...")
result_stream = Blast.qblast(
    program="blastp",
    database="refseq_protein",
    sequence=str(fasta_input.seq),
    matrix_name="BLOSUM45",
    gapcosts="13 3",
    expect=0.001,
    hitlist_size=1000,
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
result = open(f"{args.output}/results.xml", "rb")
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

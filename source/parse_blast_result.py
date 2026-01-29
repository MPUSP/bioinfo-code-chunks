import argparse
from itertools import count
from Bio import Blast  # pyright: ignore[reportMissingImports]
from Bio import SeqIO  # pyright: ignore[reportMissingImports]

# get input parameters command line arguments
parser = argparse.ArgumentParser(description="prepare flags query from BLASTP results")
parser.add_argument("-i", "--input", required=True, help="Input xml file")
parser.add_argument(
    "-o", "--output", required=True, help="Output text file with WP protein IDs"
)
parser.add_argument(
    "-m",
    "--max_hits",
    required=False,
    help="Maximum number of hits to retrieve",
    default=None,
)
parser.add_argument(
    "-f", "--filter", required=False, help="Entrez query filter", default=None
)
args = parser.parse_args()
print(f"Prepare protein IDs for Flags2 from input: {args.input}")

# parse the BLAST result
result = open(f"{args.input}", "rb")
blast_record = Blast.read(result)
print(f"Total hits found: {len(blast_record)}")

# extract seq records
filtered_result = {}
wp_number = []
count = 0
for hit in blast_record:
    for hsp in hit:
        if args.filter:
            if args.filter not in hit.target.description:
                continue
        if not hit.target.description.startswith("MULTISPECIES"):
            filtered_result[hit.target.name] = {
                "description": hit.target.description,
                "sequence": str(
                    hsp.target.seq[hsp.coordinates[0][0] : hsp.coordinates[0][1]]
                ),
                "length": len(hsp.target.seq),
                "e_value": hsp.annotations["evalue"],
                "identity": hsp.annotations["identity"],
                "alignment": hsp,
            }
            wp_number.append(hit.target.name)
            count += 1
            break
    if args.max_hits:
        if count >= int(args.max_hits):
            break
print(f"Total hits after filtering: {len(filtered_result)}")

# save the filtered results to fasta file
with open(args.output, "w") as out_stream:
    for k, v in filtered_result.items():
        out_stream.write(
            f">{k} {v['description']} length={v['length']} evalue={v['e_value']} identity={v['identity']}\n"
        )
        out_stream.write(f"{v['sequence']}\n")
print(f"Filtered BLASTP results saved to: {args.output}")

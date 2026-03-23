# imports
import csv
import json
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
import argparse


def main():

    # parse arguments
    parser = argparse.ArgumentParser(
        description="Prepare JSON files for ENA submission"
    )
    parser.add_argument(
        "--tsv", type=Path, help="Path to TSV file with run/experiment metadata"
    )
    parser.add_argument(
        "--raw_prefix", type=Path, help="Prefix path to raw fastq files", default=Path("./")
    )
    parser.add_argument(
        "--out_prefix",
        type=Path,
        default="json",
        help="Prefix path to output JSON files",
    )
    args = parser.parse_args()
    tsv_path = args.tsv
    raw_prefix = args.raw_prefix
    out_prefix = args.out_prefix

    # import TSV table with same layout as normal paired end submission, except 2 extra columns for UMI
    template = {
        "study": "",
        "sample": "",
        "name": "",
        "platform": "",
        "instrument": "",
        "insert_size": "",
        "library_name": "",
        "library_source": "",
        "library_selection": "",
        "library_strategy": "",
        "fastq": []
    }
    samples = defaultdict()
    with tsv_path.open() as f:
        content = [l for l in f.readlines() if not l.startswith("FileType")]
        reader = csv.DictReader(content, delimiter="\t")
        for row in reader:
            samples[row["sample"]] = row

    # make a list of JSON records based on the template
    records = []
    for row in samples.values():
        fastq_files = []
        if row["forward_file_name"]:
            fastq_files.append({
                "value": str(raw_prefix / row["forward_file_name"]),
                "attributes": {
                    "read_type": row["library_layout"].lower()
                }
            })
        if row["reverse_file_name"]:
            fastq_files.append({
                "value": str(raw_prefix / row["reverse_file_name"]),
                "attributes": {
                    "read_type": row["library_layout"].lower()
                }
            })
        if row["umi_file_name"]:
            fastq_files.append({
                "value": str(raw_prefix / row["umi_file_name"]),
                "attributes": {
                    "read_type": "umi_barcode"
                }
            })
        rec = deepcopy(template)
        rec["study"] = row["study"]
        rec["sample"] = row["sample"]
        rec["name"] = row["library_name"]
        rec["platform"] = row["platform"]
        rec["insert_size"] = row["insert_size"]
        rec["instrument"] = row["instrument_model"]
        rec["library_name"] = row["library_name"]
        rec["library_source"] = row["library_source"]
        rec["library_selection"] = row["library_selection"]
        rec["library_strategy"] = row["library_strategy"]
        rec["fastq"] = fastq_files
        records.append(rec)

    # export 1 json file per sample with all associated runs
    out_prefix.mkdir(exist_ok=True)
    for i, rec in enumerate(records, start=1):
        out_path = out_prefix / f"{i:02d}_{rec['sample']}.json"
        with out_path.open("w") as f:
            json.dump(rec, f, indent=2)
    print(f"Successfuly created {len(records)} JSON files in {out_prefix}.")


if __name__ == "__main__":
    main()

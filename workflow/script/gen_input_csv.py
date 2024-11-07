from pathlib import Path
from argparse import ArgumentParser
import pandas as pd


def get_args():
    parser = ArgumentParser(
        "Generate input csv files for snp_gatk workflow."
        "example: python gen_input_csv.py -I data/reads -O "
        "config --fq1 .1.fastq.gz --fq2 .2.fastq.gz"
    )
    parser.add_argument("-I", "--input_dir", type=str, help="dir path of fastq reads")
    parser.add_argument(
        "-O",
        "--output_dir",
        default="config",
        type=str,
        help="should always be config, do not modify unless you"
        " know what you are doing",
    )
    parser.add_argument(
        "--fq1", type=str, default=".1.fastq.gz", help="suffix for fastq1"
    )
    parser.add_argument(
        "--fq2", type=str, default=".2.fastq.gz", help="suffix for fastq2"
    )
    return parser.parse_args()


def main():
    args = get_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = [
        str(f.name).replace(args.fq1, "") for f in input_dir.glob(f"*{args.fq1}")
    ]
    sample_df = pd.DataFrame(samples, columns=["sample"])
    sample_df.to_csv(output_dir / "samples.tsv", index=False, sep="\t")

    fq1_paths = [str(f) for f in input_dir.glob(f"*{args.fq1}")]
    fq2_paths = [str(f) for f in input_dir.glob(f"*{args.fq2}")]

    units_df = pd.DataFrame(
        dict(
            sample=samples,
            platform=["ILLUMINA"] * len(samples),
            fq1=fq1_paths,
            fq2=fq2_paths,
        )
    )
    units_df.to_csv(output_dir / "units.tsv", index=False, sep="\t")


if __name__ == "__main__":
    main()

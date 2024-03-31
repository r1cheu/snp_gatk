from pathlib import Path
from argparse import ArgumentParser
import platform
import pandas as pd


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-I", "--input_dir", type=str) 
    parser.add_argument("-O", "--output_dir", type=str)
    parser.add_argument("-s", '--suffix', type=str, default=".1.fastq.gz")
    return parser.parse_args()


def main():
    args = get_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = [str(f.name).replace(args.suffix, '') for f in input_dir.glob(f"*{args.suffix}")]
    sample_df = pd.DataFrame(samples, columns=["sample"])
    sample_df.to_csv(output_dir / "samples.tsv", index=False, sep="\t")  
    
    fq1_paths = [str(f) for f in input_dir.glob(f"*{args.suffix}")] 
    fq2_paths = [str(f).replace(".1.fastq.gz", ".2.fastq.gz") for f in fq1_paths]

    units_df = pd.DataFrame(dict(sample=samples, platform=["ILLUMINA"] * len(samples), fq1=fq1_paths, fq2=fq2_paths))
    units_df.to_csv(output_dir / "units.tsv", index=False, sep="\t")

if __name__ == "__main__":
    main()
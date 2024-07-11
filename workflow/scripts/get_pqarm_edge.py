import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType("rb"))
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("-e", "--edge", default=20_000, type=int)
    args = ap.parse_args()

    df = pl.read_csv(args.input, has_header=False, new_columns=["chr", "start", "stop"], separator="\t")

    pl.concat([
        df.select(
            "chr",
            start=pl.col("start") - args.edge,
            stop=pl.col("start"),
            arm=pl.lit("p-arm")
        ),
        df.select(
            "chr",
            start=pl.col("stop"),
            stop=pl.col("stop") + args.edge,
            arm=pl.lit("q-arm")
        )
    ]).sort(by="chr").write_csv(args.output, include_header=False, separator="\t")

if __name__ == "__main__":
    raise SystemExit(main())
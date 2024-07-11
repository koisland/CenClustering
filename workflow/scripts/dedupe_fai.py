import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=str)
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ap.add_argument("-l", "--len", type=int, default=20_000, help="Size of region.")
    ap.add_argument("--heuristic", type=str, default="min", help="Dupe choosing heurstic. Min or max abs diff.")
    
    args = ap.parse_args()

    df = pl.read_csv(args.input, has_header=False, separator="\t", new_columns=["chr", "length", "1", "2", "3"])
    
    df = (
        df.with_columns(
            pl.col("chr").str.splitn(":", 2)
        )
        .unnest("chr")
        .rename({"field_0": "chr", "field_1": "coords"})
        .with_columns(diff=(pl.col("length") - args.len).abs())
        .filter(
            pl.when(args.heuristic == "min")
            .then(pl.col("diff") == pl.col("diff").min().over(["chr"]))
            .otherwise(pl.col("diff") == pl.col("diff").max().over(["chr"]))
        )
        .with_columns(
            chr=pl.when(pl.col("coords").is_null()).then(pl.col("chr")).otherwise(pl.col("chr") + ":" + pl.col("coords"))
        )
        .drop("diff")
    )
    df.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
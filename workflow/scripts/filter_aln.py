import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=str)
    ap.add_argument("-o", "--output", type=argparse.FileType("wt"), default=sys.stdout)
    ap.add_argument("--aln_length_thr", type=int, default=10_000)
    ap.add_argument("--perc_human_thr", type=float, default=0.01)
    args = ap.parse_args()

    df = pl.read_csv(args.input, has_header=True, separator="\t")
    df_filtered = (
        df.with_columns(
            reference_aln_length=pl.col("reference_end") - pl.col("reference_start"),
            query_aln_length=pl.col("query_end") - pl.col("query_start"),
        )
        .filter(
            (pl.col("reference_aln_length") >= args.aln_length_thr)
            & (pl.col("query_aln_length") >= args.aln_length_thr)
            & (
                pl.col("query_aln_length").sub(pl.col("reference_aln_length")).abs()
                <= pl.col("reference_length").mul(args.perc_human_thr)
            )
        )
        .drop("reference_aln_length", "query_aln_length")
        .select(
            [
                "query_name",
                "query_start",
                "query_end",
                "query_length",
                "#reference_name",
                "reference_start",
                "reference_end",
                "reference_length",
                "strand",
                "perID_by_matches",
                "perID_by_events",
                "perID_by_all",
                "matches",
                "mismatches",
                "deletion_events",
                "insertion_events",
                "deletions",
                "insertions",
            ]
        )
        .with_columns(
            pl.col("query_name").str.strip_chars_end("+-"),
            query_aln_ort=pl.col("query_name").str.extract("(\+|-)")
        )
    )
    df_filtered.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())

import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType("rb"))
    ap.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    args = ap.parse_args()

    df = pl.read_csv(
        args.input,
        has_header=False,
        new_columns=[
            "idx",
            "div",
            "deldiv",
            "insdiv",
            "contig",
            "start",
            "end",
            "left",
            "C",
            "type",
            "rClass",
            "right",
            "x",
            "y",
            "z",
            "other",
        ],
        separator=" ",
    )

    with pl.Config(tbl_cols=16):
        df = (
            df.filter(pl.col("type") != "ALR/Alpha")
            .filter(~pl.col("rClass").str.contains_any(["LINE", "SINE", "LTR"]))
            .select("contig", "start", "end", pl.col("type") + "#" + pl.col("rClass"))
        )
        df.write_csv(args.output, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())

import re
import argparse
import pandas as pd

def filter_hits(blast_results_file: str, filter_keyword: str, sample=None) -> pd.DataFrame:
    """Parse BLAST results file and generate a human-readable .tsv file
    containing viral contigs and relevant columns.

    blast_results_file: path to (tab-delimited) BLAST results file
    filter_keyword: keyword used to extract hits (ex. "Viruses")

    Returns a DataFrame containing filtered BLAST results."""

    results = pd.read_table(blast_results_file)
    keep_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "slen",
        "sscinames",
        "scomnames",
        "sblastnames",
        "sskingdoms",
        "stitle",
    ]
    results = results[keep_cols]
    results = results.drop_duplicates(subset=["qseqid"])
    results = results[
        (results["sscinames"] == filter_keyword)
        | (results["scomnames"] == filter_keyword)
        | (results["sblastnames"] == filter_keyword)
        | (results["sskingdoms"] == filter_keyword)
    ]
    results.sort_values(
        ["sscinames", "pident", "length"], ascending=[True, False, False], inplace=True
    )
    results = results[
        ["qseqid", "sscinames", "pident", "length", "qlen", "slen", "sseqid", "stitle"]
    ]
    contig_num = re.compile(
        "_[0-9]+"
    )  # get contig numbers. Expects headers like ">NODE_562_length.."
    results["qseqid"] = results["qseqid"].map(lambda x: contig_num.findall(x)[0][1:])
    access = re.compile("\|[^a-z]+\.[0-9]+")  # get NCBI accession IDs
    results["sseqid"] = results["sseqid"].map(lambda x: access.findall(x)[0][1:])
    
    if isinstance(sample, str):
        results['sample'] = sample
        results = results[['sample'] + list(results.columns[:-1])]
    
    return results


def summarize_hits(hits):
    summary = {
        "": [],
        "Virus": [],
        "Contigs": [],
        "Length": [],
        "Identity": [],
        "mxlen": [],
        "mxid": [],
    }
    for n in list(dict.fromkeys(hits["sscinames"])):
        onevirus = hits[hits["sscinames"] == n]  # DataFrame for each virus
        summary["Virus"].append(n)
        summary["Contigs"].append(onevirus.shape[0])
        summary["Length"].append(
            str(onevirus["qlen"].min()) + " - " + str(onevirus["qlen"].max())
        )
        summary["Identity"].append(
            str(onevirus["pident"].min()) + " - " + str(onevirus["pident"].max())
        )
        summary[""].append("")
        summary["mxlen"].append(onevirus["qlen"].max())
        summary["mxid"].append(onevirus["pident"].max())
    summary = pd.DataFrame(summary, columns=list(summary))
    summary.sort_values(
        ["Contigs", "mxlen", "mxid"], ascending=[False, False, False], inplace=True
    )
    summary.drop(labels=["mxlen", "mxid"], axis=1, inplace=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("blast_results", help="tab-delimited BLAST results")
    parser.add_argument("output", help="name of output file to write")
    parser.add_argument("-s", "--sample_name", type=str, help="specify a sample name to include in the output")
    args = parser.parse_args()

    header = [
        "Contig",
        "Virus",
        "Identity",
        "Hit Length",
        "Query Length",
        "Subject Length",
        "Accession",
        "Title",
    ]

    sample_name = None
    if args.sample_name:
        sample_name = args.sample_name
        header = ['Sample'] + header

    all_hits = filter_hits(args.blast_results, "Viruses", sample=sample_name)
    all_hits.to_csv(args.output, index=False, header=header, sep="\t")

    # sum_hits = summarize_hits(all_hits)
    # sum_hits.to_csv(sys.argv[2]+'_summary.csv', index=False, header=list(sum_hits))

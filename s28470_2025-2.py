#!/usr/bin/env python3
"""
Extended NCBI GenBank Retriever
"""

from Bio import Entrez, SeqIO
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        print(f"Organism: {records[0]['ScientificName']} (TaxID: {taxid})")
        term = f"txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y")
        result = Entrez.read(handle)
        self.env = result["WebEnv"]
        self.key = result["QueryKey"]
        return int(result["Count"])

    def fetch_records(self, min_len, max_len, max_records=100):
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                               retstart=0, retmax=min(max_records, 500),
                               webenv=self.env, query_key=self.key)
        text = handle.read()
        records = SeqIO.parse(StringIO(text), "genbank")
        return [{"accession": r.id, "length": len(r.seq), "description": r.description}
                for r in records if min_len <= len(r.seq) <= max_len]

def save_csv(data, taxid):
    pd.DataFrame(data).to_csv(f"taxid_{taxid}_filtered.csv", index=False)

def plot_lengths(data, taxid):
    sorted_data = sorted(data, key=lambda x: -x["length"])
    x = [d["accession"] for d in sorted_data]
    y = [d["length"] for d in sorted_data]
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, marker='o')
    plt.xticks(rotation=90)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(f"taxid_{taxid}_lengths.png")
    plt.close()

def main():
    email = input("Enter your email: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter TaxID: ")
    min_len = int(input("Min sequence length: "))
    max_len = int(input("Max sequence length: "))
    retriever = NCBIRetriever(email, api_key)
    count = retriever.search_taxid(taxid)
    print(f"Found {count} records.")
    data = retriever.fetch_records(min_len, max_len)
    if not data:
        print("No sequences matched length filter.")
        return
    save_csv(data, taxid)
    plot_lengths(data, taxid)
    print("CSV and PNG files saved.")

if __name__ == "__main__":
    main()
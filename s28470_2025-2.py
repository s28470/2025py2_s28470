from Bio import Entrez, SeqIO
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt

class GenBankRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key

    def search(self, taxid):
        tax_info = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid, retmode="xml"))
        print(f"{tax_info[0]['ScientificName']} (TaxID: {taxid})")
        search = Entrez.read(Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y"))
        self.webenv = search["WebEnv"]
        self.query_key = search["QueryKey"]
        return int(search["Count"])

    def fetch_filtered(self, min_len, max_len, limit=100):
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                               retstart=0, retmax=limit,
                               webenv=self.webenv, query_key=self.query_key)
        records = SeqIO.parse(StringIO(handle.read()), "genbank")
        return [
            {"accession": r.id, "length": len(r.seq), "description": r.description}
            for r in records if min_len <= len(r.seq) <= max_len
        ]

def save_csv(data, taxid):
    pd.DataFrame(data).to_csv(f"taxid_{taxid}_filtered.csv", index=False)

def plot_lengths(data, taxid):
    sorted_data = sorted(data, key=lambda x: -x["length"])
    x = [d["accession"] for d in sorted_data]
    y = [d["length"] for d in sorted_data]
    plt.figure(figsize=(10, 5))
    plt.plot(x, y, marker='o')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"taxid_{taxid}_lengths.png")
    plt.close()

def main():
    email = input("Email: ")
    api_key = input("NCBI API key: ")
    taxid = input("TaxID: ")
    min_len = int(input("Min sequence length: "))
    max_len = int(input("Max sequence length: "))

    retriever = GenBankRetriever(email, api_key)
    total = retriever.search(taxid)
    print(f"Found {total} records.")

    data = retriever.fetch_filtered(min_len, max_len)
    if not data:
        print("No sequences matched.")
        return

    save_csv(data, taxid)
    plot_lengths(data, taxid)
    print("Saved CSV and PNG.")

if __name__ == "__main__":
    main()
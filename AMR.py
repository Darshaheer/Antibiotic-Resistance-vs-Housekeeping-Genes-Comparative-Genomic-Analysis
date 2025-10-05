# Loading Libraries
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import pandas as pd 
from collections import Counter

# ClustalW Path
Clustal_path = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"

#Entrez Email
Entrez.email = "shaheerdar085@gmail.com"

#Sequences to Download
ARG_list = [
    "MG653169.1", "OM574584.1", "MG460317.1",
    "AY466395.1", "JF901756.1"
]
HK_list = [
    "KP670789.1", "MK564156.1", "EU896799.1",
    "EU892694.1", "AB035920.1"
]

# ----------------------------------1. Sequence Downloading ------------------------------------------
with Entrez.efetch(db= "nucleotide", id= ARG_list, rettype="gb", retmode="text") as file:
    records = list(SeqIO.parse(file, "genbank"))

SeqIO.write(records, "My_ARG.gb", "genbank")
print(f"{len(records)} sequences have been saved to the file")

with Entrez.efetch(db= "nucleotide", id= HK_list, rettype="gb", retmode="text") as file:
    record = list(SeqIO.parse(file, "genbank"))

SeqIO.write(record, "My_HKgenes.gb", "genbank")
print(f"{len(record)} sequences have been saved to the file")

# ----------------------------------2. GC % ------------------------------------------
ARG_record = list(SeqIO.parse("My_ARG.gb", "genbank"))
ARG_gc = [gc_fraction(rec.seq) * 100 for rec in ARG_record]
print(f"The GC content of Resistance genes is:\n")
print(ARG_gc)

HK_record = list(SeqIO.parse("My_HKgenes.gb", "genbank"))
HK_gc = [gc_fraction(recs.seq) * 100 for recs in HK_record]
print(f"The GC content of Housekeeping genes is:\n")
print(HK_gc)

#Saving Results
with open("GC Summary.txt", "w") as file:
    file.write("The GC Contents for Both Resistant and Housekeeping Genes is:\n")
    file.write("The Resistant Gene:\n")
    for idx, gc in enumerate(ARG_gc, 1):
        file.write(f"ARG_{idx} : {gc:.2f}%\n")
    file.write("The Housekeeping Genes:\n")
    for idx, gc in enumerate(HK_gc, 1):
        file.write(f"HK_{idx} : {gc:.2f}%")
    
print("The GC Summary file has been saved")

# ----------------------------------3. Multiple Sequence Alignment ------------------------------------------
SeqIO.convert("My_ARG.gb", "genbank", "My_ARG.fasta", "fasta")
SeqIO.convert("My_HKgenes.gb", "genbank", "My_HKgenes.fasta", "fasta")

clustalw_cline = ClustalwCommandline(Clustal_path, infile="My_ARG.fasta")
stdout, stderr = clustalw_cline()
print("ARG Alignment Complete")
ARG_alignment = AlignIO.read("My_ARG.aln", "clustal")
print("ARG Alignment Length:", ARG_alignment.get_alignment_length())

clustalw_cline = ClustalwCommandline(Clustal_path, infile="My_HKgenes.fasta")
stdout, stderr = clustalw_cline()
print("Housekeeping genes Alignment Complete")
HK_alignment = AlignIO.read("My_HKgenes.aln", "clustal")
print("HK genes Alignment Length:", HK_alignment.get_alignment_length())

AlignIO.write(ARG_alignment, "My_ARG_aligned.fasta", "fasta")
AlignIO.write(HK_alignment, "My_HKgenes_aligned.fasta", "fasta")

print("Both ARG and HK alignments saved for further analysis!")

# ----------------------------------4. Conservation Scores ------------------------------------------
ARG_alignment = AlignIO.read("My_ARG_aligned.fasta", "fasta")
HK_alignment = AlignIO.read("My_HKgenes_aligned.fasta", "fasta")

def conservation_score(alignment):
    scores = []
    aln_length = alignment.get_alignment_length()
    for i in range(aln_length):  
        column = alignment[:, i]  
        most_common = max(set(column), key=column.count)
        score = column.count(most_common) / len(column)  
        scores.append(score)
    return scores

ARG_scores = conservation_score(ARG_alignment)
HK_scores = conservation_score(HK_alignment)

# Visualization
plt.figure(figsize=(12,6))
plt.plot(ARG_scores, label="ARG", color="blue")
plt.plot(HK_scores, label="HK", color="orange")
plt.title("Conservation Across Genes (Resistant vs Housekeeping)")
plt.xlabel("Position in Alignment")
plt.ylabel("Conservation Score (0–1)")
plt.legend()
plt.tight_layout()
plt.show()

# Saving Results
with open("Conservation_summary.txt", "w") as f:
    f.write("Conservation Scores\n\n")
    f.write("ARG:\n")
    f.write(", ".join([f"{s:.2f}" for s in ARG_scores]) + "\n\n")
    f.write("HK:\n")
    f.write(", ".join([f"{n:.2f}" for n in HK_scores]) + "\n")

print("Conservation summary saved to Conservation_summary.txt")

# ----------------------------------5. Nucleotide Bias ------------------------------------------
def nucleotide_bias(fasta_file, gene_type, gc_values):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = []
    for i, record in enumerate(records):
        seq = str(record.seq).upper()
        length = len(seq)
        if length == 0:
            continue
        A = seq.count("A")
        T = seq.count("T")
        G = seq.count("G")
        C = seq.count("C")
        data.append({
            "Gene_ID": record.id,
            "Type": gene_type,
            "A%": (A / length) * 100,
            "T%": (T / length) * 100,
            "G%": (G / length) * 100,
            "C%": (C / length) * 100,
            "GC%": gc_values[i]  
        })
    return pd.DataFrame(data)


ARG_bias = nucleotide_bias("My_ARG.fasta", "Resistant", ARG_gc)
HK_bias = nucleotide_bias("My_HKgenes.fasta", "Housekeeping", HK_gc)

# Combine both tables
combined_bias = pd.concat([ARG_bias, HK_bias], ignore_index=True)

# Save to Excel
combined_bias.to_excel("Gene_Nucleotide_Composition.xlsx", index=False)

print("Nucleotide bias + GC% summary saved as 'Gene_Nucleotide_Composition.xlsx'")

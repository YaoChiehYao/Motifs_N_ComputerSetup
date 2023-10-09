"""
Author : Yaochieh Yao
This code is for assignment 2, Motifs and Computer Set-up of the
BINF6400 Genomic course. The homework requires reading a bacteria
genome as input with contigs, finding the motif score in 13 bases
of each sequence over a threshold of 7.25 to decide the start codon's
position. After that, find all possible open reading frames that end
with a stop codon, write a new fasta file, and evaluate the orfs in
the research question.
"""


def read_fasta(fasta_file):
    seqs = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.replace('>', '').replace('|test', '').strip()
            else:
                if header not in seqs:
                    seqs[header] = ''
                seqs[header] += line.strip()
    return seqs


def score_motif(motif):
    """
    This function takes in a 13bp sequence and returns the score
    """
    # letter to number matching for index
    base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    motif_score = [[.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # A
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],  # T
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],  # C
                   [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]]  # G
    indexes = [base_idx[base] for base in motif if base in base_idx]

    scores = []
    for idx, value in enumerate(indexes):
        if 0 <= value < len(motif_score) and 0 <= idx < len(motif_score[value]):
            scores.append(motif_score[value][idx])
        else:
            # Handle the case where the index is out of range
            scores.append(0)

    return sum(scores)


def find_ORFs(sequence):
    """
    This function continues with passed motifs to work on using
    the starting positions to find possible orfs with a stop
    codon in the sequence.
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for i in range(0, len(sequence)-2, 3):
        codon = sequence[i:i+3]

        if codon in stop_codons:
            orfs.append(sequence[0:i+3])

    return orfs


def scanSeq(sequence):
    """
    Takes in a sequence and returns an array of sequences,
    an array of corresponding ORF start positions, and an
    array of corresponding ORF lengths
    """
    motifs = [sequence[i:i+13]
              for i in range(len(sequence)-12) if score_motif(sequence[i:i+13]) > 7.25]
    start_pos = [i for i in range(len(sequence)-12)
                 if score_motif(sequence[i:i+13]) > 7.25]

    orfs = []

    for motif, start in zip(motifs, start_pos):
        start_indices = [j for j in range(
            len(motif)) if motif[j:j+3] in {'ATG', 'GTG'}]

        for index in start_indices:
            starting = start + index
            identified_orfs = find_ORFs(sequence[starting:])

            if identified_orfs:
                orfs.extend(identified_orfs)

    len_orf = [len(orf) for orf in orfs]

    return start_pos, len_orf, orfs


def write_fasta(filename, header_orf_dict):
    """
    The function outputs a file with a FASTA sequence of orfs
    labeled by the order of original file
    """
    with open(filename, "w") as w:
        for header, orfs in header_orf_dict.items():
            for orf in orfs:
                w.write(f">{header}\n{orf}\n")
    return filename


if __name__ == '__main__':
    """Main Business Logic"""
    fasta_seqs = read_fasta('spaceSeq.fa')
    headers = [header for header in fasta_seqs.keys()]
    sequences = [seq for seq in fasta_seqs.values()]
    header_orf_dict = {}
    for i in range(len(sequences)):
        start_pos, len_orf, orfs = scanSeq(sequences[i])
        header_orf_dict[headers[i]] = orfs
    write_fasta('output.fa', header_orf_dict)

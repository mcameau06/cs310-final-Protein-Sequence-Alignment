import argparse
from algorithms import NeedlemanWunsch, SmithWaterman
from Bio import SeqIO 
from itertools import combinations


def read_sequences(filepath):
    sequences = {}
    # each file may have multiple sequences
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def alignment_stats(a1, a2):
    matches = 0
    mismatches = 0
    gaps = 0
    for c1, c2 in zip(a1, a2):
        if c1 == '-' or c2 == '-':
            gaps += 1
        elif c1 == c2:
            matches += 1
        else:
            mismatches += 1
    total = len(a1)
    print(f"Matches:    {matches}/{total}")
    print(f"Mismatches: {mismatches}")
    print(f"Gaps:       {gaps}")


def run_alignment(args):
    try:
        seqs1 = read_sequences(args.file1)
        seqs2 = read_sequences(args.file2)

    except FileNotFoundError:
        print("A FASTA file isn't found")

    all_seqs = {**seqs1, **seqs2}

    # get alignments for all combinations of chains of sequences 

    for (name1, seq1), (name2, seq2) in combinations(all_seqs.items(), 2):
        print(f"\n{name1} vs {name2}")
        if args.algorithm == "nw":      
            aligner = NeedlemanWunsch(seq1, seq2, args.match, args.mismatch, args.gap)
        else:
            aligner = SmithWaterman(seq1, seq2, args.match, args.mismatch, args.gap)

        try:
            a1, a2 = aligner.align()

            print(f"Seq 1: {a1}")
            print(f"Seq 2: {a2}")
            alignment_stats(a1, a2)

        except Exception as e:
            print(f"Error performing alignment, {e}")

if __name__ == "__main__":
        
    parser = argparse.ArgumentParser(description="Sequence Alignment")

    parser.add_argument("--file1",type=str, help="first FASTA file")
    parser.add_argument("--file2", type=str, help="second FASTA file")
    parser.add_argument("--algorithm",choices=["nw", "sw"], help="nw or sw")
    parser.add_argument("--match", type=int, default=1, help="match score")
    parser.add_argument("--mismatch", type=int, default=-1, help="mismatch penalty")
    parser.add_argument("--gap", type=int, default=-2, help= "gap penalty")

    args = parser.parse_args()
    
    run_alignment(args)

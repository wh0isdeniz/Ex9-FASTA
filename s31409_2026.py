import random
import csv

def validate_positive_int(prompt: str, min_val: int = 1, max_val: int = 100_000) -> int:
    """Gets an integer from the user within [min_val, max_val].
    Repeats the prompt on invalid input instead of crashing."""
    while True:
        raw = input(prompt)
        try:
            value = int(raw)
            if min_val <= value <= max_val:
                return value
            else:
                print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")
        except ValueError:
            print(f"Error: value must be an integer in the range [{min_val}, {max_val}].")


# SEQUENCE GENERATION

def generate_sequence(length: int) -> str:
    """Returns a random DNA sequence of the specified length.
    Uses only standard nucleotides: A, C, G, T."""
    nucleotides = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(nucleotides) for _ in range(length))


# STATISTICS
def calculate_stats(sequence: str) -> dict:
    """Returns a dictionary of sequence statistics.
    Keys: 'A', 'C', 'G', 'T' (float, %), 'GC' (float, %).
    Must be called on the raw sequence BEFORE insert_name()."""
    length = len(sequence)
    counts = {nuc: sequence.count(nuc) for nuc in 'ACGT'}
    stats = {nuc: round(counts[nuc] / length * 100, 2) for nuc in 'ACGT'}
    stats['GC'] = round((counts['G'] + counts['C']) / length * 100, 2)
    return stats


# NAME INSERTION
def insert_name(sequence: str, name: str) -> str:
    """Inserts a name at a random position in the sequence.
    Name is written in lowercase so it is visually distinguishable
    from the uppercase nucleotides. The name does NOT affect statistics
    or declared sequence length — call calculate_stats() before this."""
    position = random.randint(0, len(sequence))  # 0..len inclusive
    return sequence[:position] + name.lower() + sequence[position:]

# FASTA FORMATTING

def format_fasta(seq_id: str, description: str, sequence: str, line_width: int = 80) -> str:
    """Returns a correctly formatted FASTA record as a string.
    Header: '>ID description' (description omitted if empty).
    Sequence is wrapped at line_width characters (last line may be shorter)."""
    # Build header line
    if description:
        header = f">{seq_id} {description}"
    else:
        header = f">{seq_id}"

    # Break sequence into fixed-width lines
    lines = [sequence[i:i + line_width] for i in range(0, len(sequence), line_width)]

    return header + '\n' + '\n'.join(lines) + '\n'


# OPTIONAL 1: MOTIF SEARCH

def find_motif(sequence: str, motif: str) -> list:
    """Searches for all occurrences of a motif in the sequence.
    Returns a list of 1-based positions (biological convention).
    Overlapping occurrences are counted."""
    positions = []
    start = 0
    # Use find() in a loop to catch overlapping matches
    while True:
        pos = sequence.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # convert to 1-based index
        start = pos + 1  # +1 allows overlapping matches
    return positions


# OPTIONAL 2: COMPLEMENTARY & REVERSE COMPLEMENTARY SEQUENCES

# Complement mapping: each base pairs with its Watson-Crick partner
COMPLEMENT_MAP = str.maketrans('ACGTacgt', 'TGCAtgca')

def complement(sequence: str) -> str:
    """Returns the complementary DNA strand (5'→3' direction preserved).
    A↔T, C↔G (Watson-Crick base pairing rules)."""
    return sequence.translate(COMPLEMENT_MAP)

def reverse_complement(sequence: str) -> str:
    """Returns the reverse complementary strand.
    This represents the antiparallel partner strand read 5'→3'."""
    return complement(sequence)[::-1]


# OPTIONAL 3: IN SILICO TRANSCRIPTION (DNA → mRNA)


def transcribe(sequence: str) -> str:
    """Generates the mRNA sequence from the coding DNA strand.
    Replaces thymine (T) with uracil (U) — the only difference
    between DNA and mRNA in terms of nucleotide composition."""
    return sequence.replace('T', 'U').replace('t', 'u')


# OPTIONAL 4: GC SLIDING WINDOW ANALYSIS

def sliding_window_gc(sequence: str, window_size: int = 10, step: int = 1) -> list:
    """Calculates GC content in a sliding window across the sequence.
    Returns a list of dicts with 'start_position' (1-based) and 'gc_content' (%).
    Only pure nucleotide characters are counted (lowercase name chars are skipped)."""
    results = []
    # Work only on uppercase nucleotide characters
    clean_seq = ''.join(c for c in sequence if c in 'ACGT')
    n = len(clean_seq)

    for i in range(0, n - window_size + 1, step):
        window = clean_seq[i:i + window_size]
        gc = round((window.count('G') + window.count('C')) / window_size * 100, 2)
        results.append({'start_position': i + 1, 'gc_content': gc})  # 1-based
    return results


def save_sliding_window_csv(results: list, filename: str) -> None:
    """Saves sliding window GC analysis results to a CSV file.
    Columns: start_position, gc_content."""
    with open(filename, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['start_position', 'gc_content'])
        writer.writeheader()
        writer.writerows(results)
    print(f"GC sliding window data saved to: {filename}")


# OPTIONAL 5: TRANSLATION (mRNA → PROTEIN)


CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate(mrna_sequence: str) -> str:
    """Translates an mRNA sequence into a protein (amino acid) sequence.
    Uses the standard genetic code. Translation stops at the first stop codon (*).
    Incomplete codons at the end are ignored.
    Returns single-letter amino acid codes."""

    clean = ''.join(c for c in mrna_sequence if c in 'ACGU')
    protein = []
    for i in range(0, len(clean) - 2, 3):  # step 3 = one codon at a time
        codon = clean[i:i + 3]
        amino_acid = CODON_TABLE.get(codon, '?')  # '?' for unknown codons
        if amino_acid == '*':  # stop codon reached
            break
        protein.append(amino_acid)
    return ''.join(protein)


# MAIN

def main():
    """Main program flow: collects user input, generates sequence,
    calculates statistics, saves FASTA file, and runs optional analyses."""

    print("=== Random DNA Sequence Generator ===\n")

    length = validate_positive_int("Enter sequence length: ")

    while True:
        seq_id = input("Enter sequence ID: ").strip()
        if seq_id and ' ' not in seq_id:
            break
        print("Error: ID cannot be empty or contain whitespace.")

    description = input("Enter a description of the sequence (optional): ").strip()
    name = input("Enter your name: ").strip()

    sequence = generate_sequence(length)
    stats = calculate_stats(sequence)
    sequence_with_name = insert_name(sequence, name)

    fasta_str = format_fasta(seq_id, description, sequence_with_name)
    fasta_filename = f"{seq_id}.fasta"
    with open(fasta_filename, 'w') as f:
        f.write(fasta_str)
    print(f"\nSequence saved to file: {fasta_filename}")

    print(f"\nSequence statistics (n={length}):")
    for nuc in 'ACGT':
        print(f"  {nuc}: {stats[nuc]:.2f}%")
    print(f"  GC-content: {stats['GC']:.2f}%")

    # --- Optional 1: Motif search ---
    print("\n--- Motif Search ---")
    motif = input("Enter a motif to search for (or press Enter to skip): ").strip().upper()
    if motif:
        positions = find_motif(sequence, motif)
        if positions:
            print(f"Motif '{motif}' found at positions (1-based): {positions}")
        else:
            print(f"Motif '{motif}' not found in the sequence.")

    # --- Optional 2: Complementary sequences (saved to FASTA) ---
    print("\n--- Complementary Sequences ---")
    comp_seq = complement(sequence)
    rev_comp_seq = reverse_complement(sequence)
    with open(fasta_filename, 'a') as f:
        f.write(format_fasta(f"{seq_id}_complement", "Complementary strand", comp_seq))
        f.write(format_fasta(f"{seq_id}_revcomp", "Reverse complementary strand", rev_comp_seq))
    print(f"Complement and reverse complement appended to {fasta_filename}")

    # --- Optional 3: Transcription ---
    print("\n--- Transcription (DNA → mRNA) ---")
    mrna = transcribe(sequence)
    with open(fasta_filename, 'a') as f:
        f.write(format_fasta(f"{seq_id}_mRNA", "mRNA transcript", mrna))
    print(f"mRNA sequence appended to {fasta_filename}")

    # --- Optional 4: GC Sliding Window ---
    print("\n--- GC Sliding Window Analysis ---")
    window_size = validate_positive_int("Enter window size (e.g. 10): ", min_val=1, max_val=length)
    gc_results = sliding_window_gc(sequence, window_size=window_size)
    csv_filename = f"{seq_id}_gc_window.csv"
    save_sliding_window_csv(gc_results, csv_filename)

    # --- Optional 5: Translation ---
    print("\n--- Translation (mRNA → Protein) ---")
    protein = translate(mrna)
    if protein:
        print(f"Protein sequence ({len(protein)} aa): {protein}")
        with open(fasta_filename, 'a') as f:
            f.write(format_fasta(f"{seq_id}_protein", "Translated protein sequence", protein))
        print(f"Protein sequence appended to {fasta_filename}")
    else:
        print("No protein could be translated (no valid codons found).")

    print("\nDone!")


if __name__ == "__main__":
    main()
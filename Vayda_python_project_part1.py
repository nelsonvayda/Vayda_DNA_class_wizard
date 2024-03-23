#!/usr/bin/env python
# coding: utf-8

# create class seq
class seq:
    # Call instance attributes
    def __init__(self, name, sequence, organism, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print(self.name)
        print(self.type)
        print(self.organism)
        print(self.sequence)

    # define the length function
    def length(self):
        x = len(self.sequence)
        print(x)

    # define the fasta_out function
    def fasta_out(self):
        with open(f"{self.name}.fa", "w") as f:
            f.write(
                ">"
                + self.name
                + "_"
                + self.organism
                + "_"
                + self.type
                + "\n"
                + self.sequence
            )


# codon translation table
standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


# amino acid molecular weights
aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# create protein class
class Protein(seq):
    def __init__(self, name, sequence, organism, type, size):
        self.size = size

        super().__init__(name, sequence, organism, type)

    # define fasta out function
    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()

    # define molecular weight function
    def mol_weight(self):
        total_weight = 0
        for aa in self.sequence:
            total_weight += aa_mol_weights.get(aa, 0)
        return total_weight


# create nucleotide class
class nucleotide(seq):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    # define gc_content function
    def gc_content(self):
        total_length = len(self.sequence)
        gc_amount = self.sequence.count("G") + self.sequence.count("C")
        gc_percent = 100 * gc_amount / total_length
        print(gc_percent)


# create DNA class
class DNA(nucleotide):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    def transcribe(self):
        transcript = self.sequence.replace("T", "U")
        print(transcript)

    # define reverse complement method
    @classmethod
    def reverse_complement(cls, dna_sequence=None, obj=None):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        # If a direct sequence is not provided, but an object is, use the object's sequence
        if dna_sequence is None and obj is not None:
            dna_sequence = obj.sequence
        return "".join([complement[base] for base in reversed(dna_sequence)])

    # define coding frames method
    def six_frames(self):
        frames = []
        # Generate three forward frames
        for i in range(3):
            frames.append(self.sequence[i:])
        # Generate reverse complement
        rev_comp_sequence = self.reverse_complement(dna_sequence=self.sequence)
        # Generate three reverse frames
        for i in range(3):
            frames.append(rev_comp_sequence[i:])
        return frames


# create RNA class
class RNA(nucleotide):
    def __init__(self, name, sequence, organism, type):

        super().__init__(name, sequence, organism, type)

    # define find start function
    def start(self):
        print(self.sequence.find("AUG"))

    # define translate function
    def translate(self):
        protein_sequence = ""
        sequence = self.sequence.replace("T", "U")  # Replace T with U for RNA
        start_index = sequence.find("AUG")

        if start_index == -1:
            return "Start codon AUG not found."

        for i in range(start_index, len(sequence), 3):
            codon = sequence[i : i + 3]
            if len(codon) != 3:
                break  # Stop if the remaining sequence is less than 3 nucleotides
            amino_acid = standard_code.get(codon, "")
            if amino_acid == "*":
                break  # Stop translation at a stop codon
            protein_sequence += amino_acid

        return protein_sequence


uidA = DNA(
    name="uidA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="bacteria",
    type="DNA",
)

uidA.fasta_out()

print(uidA.six_frames())

print(DNA.reverse_complement(obj=uidA))

uidA.transcribe()


uidA_RNA = RNA(
    name="uidA_RNA",
    sequence="CGCAUGUUACGUCCUGUAGAAACCCCAACCCGUGAAAUCAAAAAA",
    organism="bacteria",
    type="RNA",
)

uidA_RNA.fasta_out()

uidA_RNA.translate()


uidA_protein = Protein(
    name="uidA_protein",
    sequence="MLRPVETPTREIKK",
    organism="bacteria",
    type="Protein",
    size="42",
)

uidA_protein.fasta_out()

uidA_protein.mol_weight()

import pysam
import numpy as np
from Bio import SeqIO
import pandas as pd
import pytest

class Encoder:
    """
    The purpose of this class is to develop one-hot encoding for 
    nucleotide sequences and cigar strings of fasta and any number
    of bam files within bed windows, all extracted from respective files. 
    """

    def __init__(self):
        self.bed_window_list = []
        self.one_hot_dict = {nt:pos for pos, nt in enumerate("ACGT")}
        self.transition_dict = {}
        self.edge_list = [('start','M'), ('M','M'), ('M','N'), ('N','N'), ('N','M'), ('M','end')]
        self.edge_dict = {transition:i for i, transition in enumerate(self.edge_list)}
        self.transition_matrices = []
        
    def extract_exon_sequences_from_fasta(self, fasta_file, bed_file):
        """
        This def extracts nucleotide strings within each of the bed windows.

        Input: fasta and bed files
        Output: dictionary of nucleotide sequences with keys of read name,window position, strand
        """
        
        # Parse FASTA file
        fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

        # Read BED file and extract sequences for each region
        exon_sequences = {}
        with open(bed_file, "r") as bed:
            for line in bed:
                # Skip empty lines or lines that don't have enough columns
                if line.startswith("track") or len(line.strip()) == 0:
                    continue
                
                parts = line.strip().split()
                seq_id, start, end, exon_id, quality, strand = parts[:6]
                start, end = int(start), int(end)
                print(seq_id, strand)
                self.bed_window_list.append((seq_id,start,end))
                # Retrieve the sequence and slice according to the BED coordinates
                if seq_id in fasta_sequences:
                    exon_seq = fasta_sequences[seq_id].seq[start:end]
                    exon_position_key = f"{exon_id},{str(start)},{str(end)},{strand}"
                    exon_sequences[exon_position_key] = exon_seq

        return exon_sequences
        

    def bed_encoder(self, seq, reverse_flag):
        """
        This def develops one-hot encoding for the nucleotide sequences within each bed window. 

        Input:nucleotide sequences, strand boolean flag
        Output: one-hot encoded nucleotide sequence numpy matrix
        """
        transition_matrix = np.zeros((8, len(seq)), dtype=int)
    
        x = 4 if reverse_flag else 0
        for i, nt in enumerate(seq):
            transition_matrix[self.one_hot_dict[nt] + x, i] = 1
    
        return transition_matrix
        

    def generate_cigar_transitions(self, cigar_tuples, reverse_flag):
        """
        Generates a list of transitions from a CIGAR tuple list.
        Each transition tuple contains the previous and current CIGAR operations,
        with 'start' and 'end' markers at the beginning and end of the list.

        Input: list of cigar tuples from each bam alignment, flag boolean 
        Output: list of edges which depict the transition path of each CIGAR string
        """
        # Dictionary to map CIGAR operation codes to symbols
        cigar_dict = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}

        # Convert CIGAR tuples into a flat list of states
        state_list = []
        for enum, length in cigar_tuples:
            state_char = cigar_dict[enum]
            state_list.extend([state_char] * length)

        # Initialize the transition list with the start marker
        transitions = []
        # Generate the transition list as tuples of (previous_state, current_state)
        previous_state = state_list[0]
        
        for current_state in state_list[1:]:
            transitions.append((previous_state, current_state))
            previous_state = current_state

        #adding start and end edges
        finished_transition_list = [('start', state_list[0])]
        if reverse_flag: 
            finished_transition_list.extend(transitions[::-1])
            finished_transition_list.append((previous_state, 'end'))
            return finished_transition_list
        else:
            finished_transition_list.extend(transitions)
            finished_transition_list.append(((previous_state, 'end')))
            return finished_transition_list
    
    def overlap_detector(self, bam_list):
        """
        This def determines if bam alignments exist within the bed windows contained in the bed_window_list.
        Should handle bam alignments which are subsets of bed windows, and bed windows which are subsets of bam alignments. 

        Input: list of bam files
        Output: dictionary of state transition paths recorded as lists of edge tuples
        """
    #iterate through the list of N bam files
        for bam_file in bam_list:
            bam_file = pysam.AlignmentFile(bam_file, "rb")

            # Iterate through each read in the BAM file
            for read in bam_file:
                read_chrom = read.reference_name
                read_start = read.reference_start
                read_end = read.reference_end
                #flag for reverse
                reverse_flag = read.is_reverse
            
                #associate bed windows and bam alignments
                for chrom, start, end in self.bed_window_list:
                    if chrom == read_chrom and start <= read_start <= end and start <= read_end <= end:
                        cigar_tuples_list = read.cigartuples
                        if reverse_flag:
                            self.transition_dict[str(read).split("\t")[0]+'_reverse'] = self.generate_cigar_transitions(cigar_tuples_list, reverse_flag)
                        else:
                            self.transition_dict[str(read).split("\t")[0]] = self.generate_cigar_transitions(cigar_tuples_list, reverse_flag)

    def transitioner(self):
        """
        This def develops one-hot encoded numpy matrices from CIGAR state transitions paths.

        Input:dictionary of transition paths
        Output:one-hot encoded matrix of 
        """
        for key in self.transition_dict.keys():
            transition_list = self.transition_dict[key]
            print(transition_list)
            # Create a new matrix for each alignment
            transition_matrix = np.zeros((12, len(transition_list)), dtype=int)

            print(transition_matrix.shape)        
            # Encode transitions for this alignment
            x = 0 if 'reverse' not in key else 6
            for pos, transition in enumerate(transition_list):
                if transition in self.edge_dict:
                    transition_matrix[self.edge_dict[transition] + x, pos] = 1
    
            #transition_tensor = torch.tensor(transition_matrix)
            #transition_matrices.append(transition_tensor)
    
            df = pd.DataFrame(transition_matrix, index=self.edge_list*2, columns=[f"Pos_{j+1}" for j in range(transition_matrix.shape[1])])
            print(df)
            csv_filename = f"transition_matrix{key}.csv"
    
    def driver(self):
        """Driver function to manage class functon calls.
        
        Input:
        self.extract_exon_sequences_from_fasta("test_seq.fa", "test_reads.bed")
        self.bed_encoder(seq, reverse_flag)
        self.overlap_detector(['test_reads.bam', 'test_reads.bam'])
        self.transitioner()
        
        Output: function calls
        """
        
        bed_seq_dict = self.extract_exon_sequences_from_fasta("test_seq.fa", "test_reads.bed")
        for key in bed_seq_dict:
            seq = bed_seq_dict[key]
            print(seq)
            reverse_flag = False
            if '-' in key:
                reverse_flag = True
                #for reverse alignments, should the requence be the RC? 
                seq = seq.lower().replace('a','T').replace('t','A').replace('g','C').replace('c','G')[::-1]
                print(reverse_flag)
    
        encoded_seq = self.bed_encoder(seq, reverse_flag)
        
        print(key, len(encoded_seq))
        print(encoded_seq)
        
        self.overlap_detector(['test_reads.bam', 'test_reads.bam'])
        self.transitioner()

def main():
    """Access class and call driver member function to begin programatic cascade."""

    class_access = Encoder()
    class_access.driver()


if __name__ == '__main__':
    main()
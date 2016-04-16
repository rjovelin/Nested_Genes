from elegans_nested_genes import *


def convert_elegans_genome(genome_fasta):
    '''
    (file) -> dict
    Return a dictionnary with chromosome as key and sequence as value
    '''
    
    # open file for reading    
    infile = open(genome_fasta, 'r')
    # read content of fasta file
    content = infile.read()
    # close file after reading
    infile.close()
    # create a dict
    genome = {}
    # replace end of line
    content = content.replace('\n', '').upper()
    # split on >
    content = content.split('>')
    # remove empty strings
    content.remove('')    
    # populate dict
    for seq in content:
        if seq[0] == 'M':
            genome['MtDNA'] = seq.replace('MtDNA', '')
        elif seq.count('I') == 1 and 'V' not in seq:
            genome['I'] = seq.replace('I', '')
        elif seq.count('I') == 2:
            genome['II'] = seq.replace('I', '')
        elif seq.count('I') == 3:
            genome['III'] = seq.replace('I', '')
        elif seq[0] == 'X':
            genome['X'] = seq[1:]
        elif seq[0] == 'V':
            genome['V'] = seq[1:]
        elif 'IV' in seq:
            genome['IV'] = seq[2:]
    
    return genome


def take_reverse_complement(dna):
    '''
    (str) -> (str)
    Return the reverse complementary sequence of string dna

    >>> take_reverse_complement('atcg')
    'cgat'
    '''
    # create a dict with nucleotides: complements pairs 
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'Y': 'R', 'R': 'Y', 'W': 'W', 'S': 'S',
                   'K': 'M', 'M': 'K', 'D': 'H', 'V': 'B',
                   'H': 'D', 'B': 'V', 'N': 'N'}

    # convert dna in upper 
    dna2 = dna.upper()
    # create complement
    reverse_comp_dna = ''
    # loop over reversed dna, take the complement and add to reverse_comp_dna    
    for i in reversed(dna2):
        nucleotide = complements[i]
        reverse_comp_dna += nucleotide
        
    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
        
    return reverse_comp_dna


def three_prime_UTR_coordinates(gff_file):
    '''
    (file) -> dict
    Returns a dictionnary with the 3' UTR coordinates of each transcript
    '''

    # create a dictionnary to stote the coordinates {transcript_name : [chromo, start, end, orienation]}
    annotated_three_prime = {}

    # open file for reading
    infile = open(gff_file, 'r')
    for line in infile:
        # check that coordinates correspond to 3'UTR
        if 'WormBase' in line and 'three_prime_UTR' in line:
            line = line.rstrip()        
            if line != '':
                line = line.split()
                # get chromo
                chromo = line[0]
                # get start and end positions
                start, end  = int(line[3]), int(line[4])
                # get orientation
                sense = line[6]
                # get transcript name
                transcript = line[8][line[8].index(':')+1:]
                annotated_three_prime[transcript] = [chromo, start, end, sense]
                
    infile.close()
    return annotated_three_prime


def get_three_prime_UTR_sequences(gff_file, genome_seq):
    '''
    (file, file) -> dict
    Returns a dictionnary with the gene name as key and the 3'UTR sequence
    as value. Use a single transcript per gene if multiple 3'UTRs are annotated
    per gene
    '''

    # convert genome to dict
    genome = convert_elegans_genome(genome_seq)
        
    # get the 3'UTR coordinates
    UTRs = three_prime_UTR_coordinates(gff_file)
    
    # get the transcript: gene pairs
    transcripts = transcripts_to_gene(gff_file)
    
    # create dict to store transcript: UTR sequence pairs
    UTR_sequences = {}
    
    # create a set of genes for which 3'UTR sequence is recorded 
    recorded_genes = set()
    
    # loop over UTRs:
    for transcript in UTRs:
        # get the coordinates
        chromo = UTRs[transcript][0]
        # use 0-based index
        start = UTRs[transcript][1] -1 
        end = UTRs[transcript][2]
        orientation = UTRs[transcript][-1]
        # grab the corresponding sequence
        UTR_seq = genome[chromo][start : end]
        # if UTR on minus strand, take reverse complement
        if orientation == '-':
            UTR_seq = take_reverse_complement(UTR_seq)
        # get the gene name
        gene = transcripts[transcript]
        # check if the parent gene's UTR is already recorded in dict
        if gene not in recorded_genes:
            # populate the UTR dict with transcript: UTR_seq
            UTR_sequences[gene] = UTR_seq
            # add the parent gene to set of recorded genes
            recorded_genes.add(gene)
    
    return UTR_sequences
        

def UTR_table_for_targetscan(gff_file, genome_seq, outputfile):
    '''
    (file, file, file) -> file
    Generate the input file for TargetScan with the 3'UTR sequences of
    Celegans transcripts. Use a single transcript per gene. 
    '''
        
    # get the UTR sequences
    UTR_sequences = get_three_prime_UTR_sequences(gff_file, genome_seq)    
    
    # open outfile for writing (transcript_ID species_ID UTR_seq)
    newfile = open(outputfile, 'w')
    for gene in UTR_sequences:
        # convert T to U
        newfile.write(gene + '\t' + '6239' + '\t' + UTR_sequences[gene].replace('T', 'U') + '\n')
    
    newfile.close()  


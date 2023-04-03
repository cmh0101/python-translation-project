#! /usr/bin/env python3

import sys
import re

def vet_nucleotide_sequence(sequence):
    """
    Return None if `sequence` is a valid RNA or DNA sequence, else raise exception. 

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if sequence is valid, otherwise raise an
        exception.

    Examples
    --------
    >>> vet_nucleotide_sequence('ACGTACGT') == None
    True

    >>> vet_nucleotide_sequence('not a valid sequence')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'not a valid sequence'

    Don't allow mixing of DNA and RNA!
    >>> vet_nucleotide_sequence('AUTGC')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: 'AUTGC'

    Don't allow whitespace (or other characters) before, within, or after!
    >>> vet_nucleotide_sequence(' ACGT ACGT ')
    Traceback (most recent call last):
        ...
    Exception: Invalid sequence: ' ACGT ACGT '

    But, an empty string should be deemed valid
    >>> vet_nucleotide_sequence('') == None
    True
    """


    rna_pattern_str = r'^[^T]*$' 
    dna_pattern_str = r'^[^U]*$'

    rna_pattern = re.compile(rna_pattern_str)
    dna_pattern = re.compile(dna_pattern_str)

    if rna_pattern.match(sequence):
        return
    if dna_pattern.match(sequence):
        return
    else:
        raise Exception("Invalid sequence: {0!r}".format(sequence))



def vet_codon(codon):
    """
    Return None if `codon` is a valid RNA codon, else raise an exception. 

    Parameters
    ----------
    codon : str
        A string representing a codon (upper or lower-case)

    Returns
    -------
    None
        Return nothing (None) if codon is valid, otherwise raise an
        exception.

    Examples
    --------
    Valid codon
    >>> vet_codon('AUG') == None
    True

    lower-case is also vaild 
    >>> vet_codon('aug') == None
    True

    DNA is not valid
    >>> vet_codon('ATG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'ATG'

    A codon must be exactly 3 RNA bases long
    >>> vet_codon('AUGG')
    Traceback (most recent call last):
        ...
    Exception: Invalid codon: 'AUGG'
    """

    codon_pattern_str = r'([aug]{3,3}$|[AUG]{3,3}$)'
    codon_pattern = re.compile(codon_pattern_str)

    if codon_pattern.match(codon):
        return
    else:
        raise Exception("Invalid codon: {0!r}".format(codon))


def find_first_orf(sequence,
        start_codons = ['AUG'],
        stop_codons = ['UAA', 'UAG', 'UGA']):
    """
    Return the first open-reading frame in the DNA or RNA `sequence`.

    An open-reading frame (ORF) is the part of an RNA sequence that is
    translated into a peptide. It must begin with a start codon, followed by
    zero or more codons (triplets of nucleotides), and end with a stop codon.
    If there are no ORFs in the sequence, an empty string is returned.

    Parameters
    ----------
    sequence : str
        A string representing a DNA or RNA sequence (upper or lower-case)
    start_codons : list of strings
        All possible start codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.
    stop_codons : list of strings
        All possible stop codons. Each codon must be a string of 3 RNA bases,
        upper or lower-case.

    Returns
    -------
    str
        An uppercase string of the first ORF found in the `sequence` that
        starts with any one of the `start_codons` and ends with any one of the
        `stop_codons`. If no ORF is found an empty string is returned.

    Examples
    --------
    When the whole RNA sequence is an ORF:
    >>> find_first_orf('AUGGUAUAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When the whole DNA sequence is an ORF:
    >>> find_first_orf('ATGGTATAA', ['AUG'], ['UAA'])
    'AUGGUAUAA'

    When there is no ORF:
    >>> find_first_orf('CUGGUAUAA', ['AUG'], ['UAA'])
    ''

    When there is are bases before and after ORF:
    >>> find_first_orf('CCAUGGUAUAACC', ['AUG'], ['UAA'])
    'AUGGUAUAA'
    """
    # Make sure the sequence is valid
    vet_nucleotide_sequence(sequence)

    # Make sure the codons are valid
    for codon in start_codons:
        vet_codon(codon)
    for codon in stop_codons:
        vet_codon(codon)

    # Get copies of everything in uppercase
    seq = sequence.upper()
    starts = [c.upper() for c in start_codons]
    stops = [c.upper() for c in stop_codons]
    # Make sure seq is RNA
    seq = seq.replace('T', 'U')


    orf_pattern_str = r"(?:" + "|".join(start_codons) + ")" + r"(?:[AUGC]{3})*?" + r"(?:" + "|".join(stop_codons) + ")"
 
    # Create the regular expression object
    orf_pattern = re.compile(orf_pattern_str)
    # Search the sequence
    match_object = orf_pattern.search(seq)
    if match_object:
        return match_object.group()
    return ''


def parse_sequence_from_path(path):
    # Try to open the path to read from it, and handle exceptions if they
    # arrise
    try:
        file_stream = open(path, 'r')
    except FileNotFoundError as e:
        sys.stderr.write("Sorry, couldn't find path {}".format(path))
        raise e
    except IsADirectoryError as e:
        sys.stderr.write("Sorry, path {} appears to be a directory".format(
                path))
        raise e
    except:
        sys.stderr.write("Sorry, something went wrong when trying to open {}".format(
                path))
        raise
    # If we've reached here, the file is open and ready to read
    sequence = ''
    # A for loop to visit each line in the file
    for line in file_stream:
        # Strip whitespace from the line and concatenate it to the end of the
        # sequence
        sequence += line.strip()
    return sequence


def main():
    import argparse

    # Create a command-line parser object
    parser = argparse.ArgumentParser()

    default_start_codons = ['AUG']
    default_stop_codons = ['UAA', 'UAG', 'UGA']

    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A start codon. This option can be used multiple times '
                    'if there are multiple start codons. '
                    'Default: {0}.'.format(" ".join(default_start_codons))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('A stop codon. This option can be used multiple times '
                    'if there are multiple stop codons. '
                    'Default: {0}.'.format(" ".join(default_stop_codons))))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    # Check to see if start/stop codons were provided by the caller. If not,
    # use the defaults.
    if not args.start_codon:
        args.start_codon = default_start_codons
    if not args.stop_codon:
        args.stop_codon = default_stop_codons

    orf = find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)
    sys.stdout.write('{}\n'.format(orf))


if __name__ == '__main__':
    main()
def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """
    rna_sequence = rna_sequence.upper()
    amino = []
    base = list(rna_sequence)
    while True:
        if len(base) > 2:
            codon = ''.join(base[0:3])
            del base[0:3]
        else:
            break
        aa = genetic_code[codon]
        if aa == '*':
            break
        amino.append(aa)
    return ''.join(amino)

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.
        
    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).
    
    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.
        
    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.
            
    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).
    
    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    
    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    rna_sequence = rna_sequence.upper()
    num_base = len(rna_sequence)
    last_first_base_index = num_base - 3
    
    polypeptide_list = []
    for i in range(last_first_base_index + 1):
        i_end = i + 3
        next_three = rna_sequence[i:i_end]
        if next_three == 'AUG':
            polypeptide = translate_sequence(rna_sequence[i:], genetic_code)
            polypeptide_list.append(polypeptide)
    return polypeptide_list
        
def get_reverse(sequence):
    """Reverse orientation of `sequence`.
        
    Returns a string with `sequence` in the reverse order.
        
    If `sequence` is empty, an empty string is returned.
    
    Examples
    --------
    >>> get_reverse('AUGC')
  get_reverse = sequence[::-1]
    return get_reverse.upper()

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.
    
    Returns a string with the complementary sequence of `sequence`.
    
    If `sequence` is empty, an empty string is returned.
    
    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    
    sequence = sequence.upper()
    
    complement = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    
    compsequence = ''
       
    for i in sequence:
        compsequence += complement[i]
    return compsequence
    
def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.
        
    Returns a string that is the reversed and complemented sequence
    of `sequence`.
            
    If `sequence` is empty, an empty string is returned.
    
    Examples
    --------
    >>> reverse_and_complement('AUGC')   
    'GCAU'
    """
        
    compsequence = get_complement(sequence)
    rev_compsequence = get_reverse(compsequence)
    return rev_compsequence
    
def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.
    
    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.
    
    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).
    
    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').
    
    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
        
    revcomp_rna_sequence = reverse_and_complement(rna_sequence)
    revcomp_rna_array = get_all_translations(revcomp_rna_sequence, genetic_code)
    rna_array = get_all_translations(rna_sequence, genetic_code)
    longest_array = rna_array
        
    maxl = -1
    longest = ''  
    index = -1
    
    
    for i in range(len(rna_array)):
        test = len(rna_array[i])
        if test > maxl:
           maxl = test
           index = i
           longest = rna_array[index]
    
    for j in range(len(revcomp_rna_array)):
        revcomp_test = len(revcomp_rna_array[j])
        if revcomp_test > maxl:
           maxl = revcomp_test
           index = j
           longest_array = revcomp_rna_array
           longest = longest_array[index]
    
    return longest
    
if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': $
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a s$
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(  
            rna_seq,
            longest_peptide) 
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
    
    


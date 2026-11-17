# rosalind_solutions2
problems to be solved:

DNA - given a DNA string of a certain length, goal is to return a count of how many of each nucleotide there are in the order A, C, G, T; in the format #A #C #G #T.

RNA - given a DNA string of a certain length, goal is to return a RNA string which is a transcribed twin of the DNA string. This means to ignore 5'/3' stuff and essentially just swap T's with U's.

REVC - given a DNA string of a certain length, goal is to return the reverse complement of the DNA string. If the given string is arranged from 5' to 3', the return string should be arranged 5' to 3'.

GC - given a long, single string made up of many DNA sequences, separated by FASTA ids, goal is to return the FASTA id of the sequence with the largest GC content, along with its corresponding GC content, expressed as a float number between 0 to 100.

HAMM - given two DNA strings of equal length, goal is to return the Hamming distance between the strings.

PROT - given a string of RNA, goal is to return the string of amino acids to which this would correspond if the RNA were translated into protein.

SUBS - given two strings, where the first is longer than the second and the second is meant to be a substring of the first, goal is to return every starting position where the substring can be found the longer string. The assumption is that positions are 1-based instead of 0-based. Python is natively 0-based.

PRTM - given a amino acid string, goal is to return the atomic mass of the peptide, rounded to the thousandth.

REVP - given a DNA string of a certain length in FASTA format, goal is to return the starting position and length of every reverse palindrome within the string which is between 4 to 12 characters long. The pairs may return in any order, but should be separated by newlines, while the reported position and length are separated by a single space.

TRAN - given two DNA strings of equal length in FASTA format, goal is to return the transition/transversion ratio of string 1 to string 2.

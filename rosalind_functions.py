#DNA
def base_count(x):
    dict_nuc_count = {}
    for i in x:
        if i in dict_nuc_count:
            dict_nuc_count[i] += 1
        else:
            dict_nuc_count[i] = 1
    ret_ans = str(dict_nuc_count['A']) + ' ' + str(dict_nuc_count['C']) + ' ' + str(dict_nuc_count['G']) + ' ' + str(dict_nuc_count['T'])
    return ret_ans

#RNA
def DNA_txn(x):
    RNA = x.replace('T', 'U')
    return RNA

#REVC
def rev_compl(x):
    rev_compl = x.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    return rev_compl

#GC
def gc_content(x):
    ct = 0
    strhold = ''
    idhold = ''
    full_dict = {}
    while ct < len(x):
        if x[ct] == '>':
            idhold = x[ct:ct+14]
            strhold = ''
            full_dict[idhold] = strhold
        elif x[ct] != '>':
            strhold = strhold + x[ct]
            full_dict[idhold] = strhold
        ct += 1
    for i in full_dict:
        full_dict[i] = full_dict[i].strip('>Rosalind_0987654321\n').replace('\n','')
    dict_gc = {}
    for i in full_dict:
        dict_gc[i] = round(float(full_dict[i].count('G') + full_dict[i].count('C'))/len(full_dict[i]), 6)
    return str(max(dictgc, key = dictgc.get)) + '\n' + str(max(dictgc.values())*100)

#HAMM
def hamm_dist(x,y):
    x_pos = range(len(x))
    y_pos = range(len(y))
    hd = 0
    for i in x_pos:
        if x[i] != y[i]:
            hd += 1
    return hd

#PROT
def transl8n(x):
    codons = [x[i:i+3] for i in range(0, len(x), 3)]
    prot_seq = []
    for t in codons:
        if 'UUU' in t:
            prot_seq.append('F')
        if 'UUC' in t:
            prot_seq.append('F')
        if 'UUA' in t:
            prot_seq.append('L')
        if 'UUG' in t:
            prot_seq.append('L')
        if 'UCU' in t:
            prot_seq.append('S')
        if 'UCC' in t:
            prot_seq.append('S')
        if 'UCA' in t:
            prot_seq.append('S')
        if 'UCG' in t:
            prot_seq.append('S')
        if 'UAU' in t:
            prot_seq.append('Y')
        if 'UAC' in t:
            prot_seq.append('Y')
        if 'UAA' in t:
            prot_seq.append('')
        if 'UAG' in t:
            prot_seq.append('')
        if 'UGA' in t:
            prot_seq.append('')
        if 'UGU' in t:
            prot_seq.append('C')
        if 'UGC' in t:
            prot_seq.append('C')
        if 'UGG' in t:
            prot_seq.append('W')
        if 'CUU' in t:
            prot_seq.append('L')
        if 'CUC' in t:
            prot_seq.append('L')
        if 'CUA' in t:
            prot_seq.append('L')
        if 'CUG' in t:
            prot_seq.append('L')
        if 'CCU' in t:
            prot_seq.append('P')
        if 'CCC' in t:
            prot_seq.append('P')
        if 'CCA' in t:
            prot_seq.append('P')
        if 'CCG' in t:
            prot_seq.append('P')
        if 'CAU' in t:
            prot_seq.append('H')
        if 'CAC' in t:
            prot_seq.append('H')
        if 'CAA' in t:
            prot_seq.append('Q')
        if 'CAG' in t:
            prot_seq.append('Q')
        if 'CGU' in t:
            prot_seq.append('R')
        if 'CGC' in t:
            prot_seq.append('R')
        if 'CGA' in t:
            prot_seq.append('R')
        if 'CGG' in t:
            prot_seq.append('R')
        if 'AUU' in t:
            prot_seq.append('I')
        if 'AUC' in t:
            prot_seq.append('I')
        if 'AUA' in t:
            prot_seq.append('I')
        if 'AUG' in t:
            prot_seq.append('M')
        if 'ACU' in t:
            prot_seq.append('T')
        if 'ACC' in t:
            prot_seq.append('T')
        if 'ACA' in t:
            prot_seq.append('T')
        if 'ACG' in t:
            prot_seq.append('T')
        if 'AAU' in t:
            prot_seq.append('N')
        if 'AAC' in t:
            prot_seq.append('N')
        if 'AAA' in t:
            prot_seq.append('K')
        if 'AAG' in t:
            prot_seq.append('K')
        if 'AGU' in t:
            prot_seq.append('S')
        if 'AGC' in t:
            prot_seq.append('S')
        if 'AGA' in t:
            prot_seq.append('R')
        if 'AGG' in t:
            prot_seq.append('R')
        if 'GUU' in t:
            prot_seq.append('V')
        if 'GUC' in t:
            prot_seq.append('V')
        if 'GUA' in t:
            prot_seq.append('V')
        if 'GUG' in t:
            prot_seq.append('V')
        if 'GCU' in t:
            prot_seq.append('A')
        if 'GCC' in t:
            prot_seq.append('A')
        if 'GCA' in t:
            prot_seq.append('A')
        if 'GCG' in t:
            prot_seq.append('A')
        if 'GAU' in t:
            prot_seq.append('D')
        if 'GAC' in t:
            prot_seq.append('D')
        if 'GAA' in t:
            prot_seq.append('E')
        if 'GAG' in t:
            prot_seq.append('E')
        if 'GGU' in t:
            prot_seq.append('G')
        if 'GGC' in t:
            prot_seq.append('G')
        if 'GGA' in t:
            prot_seq.append('G')
        if 'GGG' in t:
            prot_seq.append('G')
    return ''.join(prot_seq)

#SUBS
def substr_pos(x, y):
    x_pos = range(len(x))
    list_pos = []
    for i in x_pos:
        if x[i:i+len(y)] == y:
            list_pos.append(str(i+1))
    return ' '.join(list_pos)

#PRTM
def pep_mass(x):
    mass = 0
    aa_dict = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
               'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
               'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
               'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
               'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333}
    for i in x:
        if i in aa_dict:
            mass = mass + aa_dict[i]
    return round(mass, 3)

#REVP
def rev_pal_pl(x):
    y = x.strip('>Rosalind_1234567890\n').replace('\n','')
    import string
    def rev_compl(z):
        return z.translate(string.maketrans('ATCG', 'TAGC'))[::-1]
    def rev_pal_chk(z):
        return z == rev_compl(z)
    pos_len_list = []
    for i in range(4, 13):
        for o  in range(len(y) - i + 1):
            substr = y[o:o+i]
            if rev_pal_chk(substr):
                pos_len_list.append((o+1, i))
    ans = ''
    for i in range(len(pos_len_list)):
        ans = ans + str(pos_len_list[i][0]) + ' ' + str(pos_len_list[i][1])
        ans = ans + '\n'
    return ans

#TRAN
def tran_tran_rat(x):
    ct = 0
    strhold = ''
    idhold = ''
    seq_dict = {}
    while ct < len(x):
        if x[ct] == '>':
            idhold = x[ct:ct+14]
            strhold = ''
            seq_dict[idhold] = strhold
        elif x[ct] != '>':
                strhold = strhold + x[ct]
                seq_dict[idhold] = strhold
        ct += 1
    seq_list = []
    for i in seq_dict:
            seq_dict[i] = seq_dict[i].strip('>Rosalind_0987654321\n').replace('\n','')
            seq_list.append(seq_dict[i])
    zero_pos = range(len(seq_list[0]))
    one_pos = range(len(seq_list[1]))
    transt = 0
    transv = 0
    for i in zero_pos:
        if seq_list[0][i] == seq_list[1][i]:
            continue
        elif seq_list[0][i] == 'A' and seq_list[1][i] == 'G' or seq_list[0][i] == 'G' and seq_list[1][i] == 'A':
            transt = transt + 1
        elif seq_list[0][i] == 'C' and seq_list[1][i] == 'T' or seq_list[0][i] == 'T' and seq_list[1][i] == 'C':
            transt = transt + 1
        elif seq_list[0][i] == 'A' and seq_list[1][i] == 'C' or 'T' or seq_list[0][i] == 'G' and seq_list[1][i] == 'C' or 'T':
            transv = transv + 1
        elif seq_list[0][i] == 'C' and seq_list[1][i] == 'A' or 'G' or seq_list[0][i] == 'T' and seq_list[1][i] == 'A' or 'G':
            transv = transv + 1
    return float(transt)/transv
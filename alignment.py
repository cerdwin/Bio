#!/usr/bin/python
import argparse
import csv
import pprint
#from bio import SeqIO
import string
from collections import defaultdict
 
import numpy as np
from Bio import SeqIO
 
 
def csv_dict_list(argument):
    reader = csv.DictReader(open(str(args.e), 'r'))
    dict_list = []
    for line in reader:
        dict_list.append(line)
    return dict_list
 
 
def find_me(mapper, a, b):
    for key, value in mapper.items():
        if key == (a, b):
            return value
 
 
def find_max(x, y, my_table, seqencea, sequenceb, my_grid, penalty):
    find = find_me(my_grid, seqencea[x-1], sequenceb[y-1])
    A = int(my_table[x-1][y])-penalty
    B = int(my_table[x-1][y-1])+find
    C = int(my_table[x][y-1])-int(penalty)
 
 
    if A >= B and A >= C:
        return A, 1
    if B >= A and B >= C:
        return B, 2
    if C >= A and C >= B:
        return C, 3
    #if A >= B and A >= C:
      #  return A, 1
   # elif C >= B and A >= C:
     #   return A, 1
    #elif B >= A and B >= C:
     #   return B, 2
    #else:
      #  return C, 3
 
 
def global_trace(x, y, track, sequenceb, seqencea):
    if track[x][y] == 1:
        return x-1, y, seqencea[x - 2], sequenceb[y-1]
    elif track[x][y] == 2:
        return x-1, y-1, seqencea[x - 2], sequenceb[y-2]
    elif track[x][y] == 3:
        return x, y-1, seqencea[x - 1], sequenceb[y-2]
    else:
        return -1, -1, '-1', '-1'
 
 
def my_global(penalisation, seqencea, sequenceb, matrix):
 
    first = []
    second = []
    for i in range(len(seqencea)):
        first.append(seqencea[i])
 
    for i in range(len(sequenceb)):
        second.append(sequenceb[i])
 
    my_table = [[0 for i in range(len(seqencea)+1)] for j in range(len(sequenceb) + 1)]
    trackback = [[0 for i in range(len(seqencea)+1)] for j in range(len(sequenceb) + 1)]
 
    for i in range(1, len(seqencea)+1):
        my_table[0][i] = my_table[0][i-1]-penalisation
        trackback[0][i] = 3
 
    for i in range(1, len(sequenceb)+1):
        my_table[i][0] = my_table[i-1][0]-penalisation
        trackback[i][0] = 1
 
    my_table[0][0] = 0
    trackback[0][0] = -1
    a = 0
    b = 0
 
    for x in range(1, len(sequenceb) + 1):
        for y in range(1, len(seqencea)+1):
            table_result, trackback[x][y] = find_max(x, y, my_table, sequenceb, seqencea, matrix, penalisation)
            my_table[x][y] = str(table_result)
 
    first_sequence = []
    second_sequence = []
    besta = []
    bestb = []
 
    first_sequence.append(seqencea[y - 1])
    second_sequence.append(sequenceb[x - 1])
 
    besta.append(seqencea[-1])
    bestb.append(sequenceb[-1])
 
    str_a = ''
    str_b = ''
    orientation = trackback[x][y]
    result = int(my_table[x][y])
    finish = int(my_table[x][y])
    while True:
        original_x = x
        original_y = y
        x, y, str_a, str_b = global_trace(x, y, trackback, seqencea, sequenceb)
        if x == 0 and y == 0:
            break
        elif x == 0 and y != 0:
            if original_y == y:
                besta[-1] = '-'
            besta.append(seqencea[y - 1])
            bestb.append('-')
        elif y == 0 and x != 0:
            if original_x == x:
                bestb[-1] = '-'
            bestb.append(sequenceb[x-1])
            besta.append('-')
        elif original_x != x and original_y != y:
            besta.append(seqencea[y-1])
            bestb.append(sequenceb[x-1])
        elif original_x == x and original_y != y:
            besta.append(seqencea[y-1])
            bestb[-1] = '-'
            bestb.append(sequenceb[x - 1])
        elif original_x!= x and original_y == y:
            besta[-1] = '-'
            besta.append(seqencea[y-1])
            bestb.append(sequenceb[x-1])
 
        #result = result + int(my_table[x][y])
 
        if x == 0 and y == 0:
            break
 
    for i in range(1, 6):
        for y in range(1, 10):
            print(trackback[i][y], ', ')
        print('\n')
 
    # print(np.matrix(my_table))
    #print(trackback)
    besta.reverse()
    bestb.reverse()
    print(''.join(besta))
    print('and')
    print(''.join(bestb))
    print(finish)
 
 
 
 
def find_local_max(x, y, my_table, seqencea, sequenceb, my_grid, penalty):
 
    find = find_me(my_grid, seqencea[x-1], sequenceb[y-1])
    A = int(my_table[x-1][y])-penalty
    B = int(my_table[x-1][y-1])+find
    C = int(my_table[x][y-1])-int(penalty)
    D = 0
    if A >= B and A >= C and A >= D:
        return A, 1
    elif B >= A and B >= C and B >= 0:
        return B, 2
    elif D >= A and D >= B and D >= C:
        return D, -1
    else:
        return C, 3
 
 
def my_local(penalisation, seqencea, sequenceb, matrix):
 
    first = []
    second = []
    # i place strings into arrays
 
    for i in range(len(seqencea)):
        first.append(seqencea[i])
 
    for i in range(len(sequenceb)):
        second.append(sequenceb[i])
 
    # I create a table where Im going to fill alignment scores
    my_table = [[0 for i in range(len(seqencea)+1)] for j in range(len(sequenceb) + 1)]
    trackback = [[0 for i in range(len(seqencea)+1)] for j in range(len(sequenceb) + 1)]
 
    # first row and column are zeros
    my_table[0][0] = 0
    for i in range(1, len(seqencea)+1):
        my_table[0][i] = my_table[0][i-1]
    for i in range(1, len(sequenceb)+1):
        my_table[i][0] = my_table[i-1][0]
 
 
    minimum = -999999999999999999
    a = 0
    b = 0
    for x in range(1, len(sequenceb) + 1):
        for y in range(1, len(seqencea)+1):
            table_result, trackback[x][y] = find_local_max(x, y, my_table, sequenceb, seqencea, matrix, penalisation)
            my_table[x][y] = str(table_result)
            if int(my_table[x][y]) > minimum:
                minimum = int(my_table[x][y])
                a = x
                b = y
 
    ###############
    # pprint.pprint(np.matrix(my_table))
    first_sequence = []
    second_sequence = []
    besta = []
    bestb = []
    first_sequence.append(seqencea[b-1])
    second_sequence.append(sequenceb[a-1])
    str_a = ''
    str_b = ''
    x = a
    y = b
    result = int(my_table[x][y])
    besta.append(seqencea[y-1])
    bestb.append(sequenceb[x-1])
    while True:
        original_x = x
        original_y = y
        x, y, str_a, str_b = global_trace(x, y, trackback, seqencea, sequenceb)
        if my_table[x][y] == str(0):
            break
        elif x == 0 and y != 0:
            if original_y == y:
                besta[-1] = '-'
            besta.append(seqencea[y - 1])
            bestb.append('-')
        elif y == 0 and x != 0:
            if original_x == x:
                bestb[-1] = '-'
            bestb.append(sequenceb[x-1])
            besta.append('-')
 
        elif original_x != x and original_y != y:
            besta.append(seqencea[y-1])
            bestb.append(sequenceb[x-1])
        elif original_x == x:
            besta.append(seqencea[y-1])
            bestb[-1] = '-'
            bestb.append(sequenceb[x - 1])
        else:
            #
            besta[-1] = '-'
            besta.append(seqencea[y-1])
            bestb.append(sequenceb[x-1])
 
        # result = result + int(my_table[x][y])
        if my_table[x][y] == str(0):
            break
 
    # pprint.pprint(np.matrix(my_table))
    # pprint.pprint(trackback)
    besta.reverse()
    bestb.reverse()
    print(''.join(besta))
    print(''.join(bestb))
    print(result)
 
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", help="global", action='store_true')
    parser.add_argument("-l", help="local", action='store_true')
    parser.add_argument("-s1", help="the path to first sequence")
    parser.add_argument("-s2", help="path to second sequence")
    parser.add_argument("-e", help="path to the score matrix in CSV format follows")
    parser.add_argument("-p", help="gap penalization")
    parser.add_argument('-pe', help="gap_extension_penalisation", action='store_true')
 
    args = parser.parse_args()
    # extracting ranking table
    table2 = []
    my_matrix = defaultdict(int)
 
    with open(str(args.e), newline='') as csv_source:
        matrixReader = csv.reader(csv_source, delimiter=',', quotechar='"')
        names = next(matrixReader)
        for i in range(1, len(names)):
            row = next(matrixReader)
            for j in range(1, len(names)):
                my_matrix[(row[0], names[j])] = int(row[j])
 
    csv_dictionary = csv_dict_list(args.e)
 
    ##### PARSING FROM FASTA #####
    for s1 in SeqIO.parse(str(args.s1), "fasta"):
        seq1 = s1.seq
    for s2 in SeqIO.parse(str(args.s2), "fasta"):
        seq2 = s2.seq
    print('shut the heck up')
    print(seq1)
    print('insufferable pain')
    print(seq2)
    print('done')
 
    ######################################
    seq1st_examplea = 'HEAGAWGHEE'#"CGTCCAAGTG"# #PEKREMTECHLSDEEIRKLNRDLRILIATNGTLTRILNVLANDEIVVEIVKQQIQDAAPEMDGCDHSSIGRVLRRDIVLKGRRSGIPFVAAESFIAIDLLPPEIVASLLETHRPIGEVMAASCIETFKEEAKVWAGESPAWLELDRRRNLPPKVVGRQYRVIAEGRPVIIITEYFLRSVFEDNSREEPIRHQRSVGTSARSGRSICT"
    seq2st_exampleb = 'PAWHEAE'#"TACGAA"#  #LSREEIRKLDRDLRILVATNGTLTRVLNVVANEEIVVDIINQQLLDVAPKIPELENLKIGRILQRDILLKGQKSGILFVAAESLIVIDLLPTAITTYLTKTHHPIGEIMAASRIETYKEDAQVWIGDLPCWLADYGYWDLPKRAVGRRYRIIAGGQPVIITTEYFLRSVFQDTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL"
    #seq2 = 'MTNRTLSREEIRKLDRDLRILVATNGTLTRVLNVVANEEIVVDIINQQLLDVAPKIPELENLKIGRILQRDILLKGQKSGILFVAAESLIVIDLLPTAITTYLTKTHHPIGEIMAASRIETYKEDAQVWIGDLPCWLADYGYWDLPKRAVGRRYRIIAGGQPVIITTEYFLRSVFQDTPREELDRCQYSNDIDTRSGDRFVLHGRVFKNL'
    #seq1 = 'MLAVLPEKREMTECHLSDEEIRKLNRDLRILIATNGTLTRILNVLANDEIVVEIVKQQIQDAAPEMDGCDHSSIGRVLRRDIVLKGRRSGIPFVAAESFIAIDLLPPEIVASLLETHRPIGEVMAASCIETFKEEAKVWAGESPAWLELDRRRNLPPKVVGRQYRVIAEGRPVIIITEYFLRSVFEDNSREEPIRHQRSVGTSARSGRSICT'
    #seq1 = 'HEAGAWGHEE'
    #seq2 = 'PAWHEAE'
    #my_local(8, seq1, seq2, my_matrix)
    #print(seq1)
    #print('sesesd')
    #print(seq2)
    #print('leggo')
    if args.l:
        my_local(int(args.p), seq1, seq2, my_matrix)
    elif args.g:
        my_global(int(args.p), seq1, seq2, my_matrix)

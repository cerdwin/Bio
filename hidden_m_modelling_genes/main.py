#!/usr/bin/python

import sys
import math

def parse_test():
    '''
    Takes in input files and returns a list of segments along with coding positions and total size of genes to test for
    :return:
    '''
    size = 0
    global_coordinates = []
    test_file = open('test.seqs', 'r')
    for segment in test_file:
        tmp = [i.strip(',[]') for i in segment.split()]
        current_coordinates = {}
        current_coordinates['from'] = current_coordinates['to'] = -1
        current_coordinates['from'] += int(tmp[0])
        current_coordinates['to'] += int(tmp[1])
        current_coordinates['positions'] = []
        coordinate_list = []
        for i in tmp[2:]:
            coordinate_list.append(int(i)-1)
            if len(coordinate_list) ==2:
                size+=coordinate_list[1]-coordinate_list[0]
                current_coordinates['positions'].append(tuple(coordinate_list))
                coordinate_list = []
        global_coordinates.append(current_coordinates)
    return global_coordinates, size
def parse_input(): # checked
    '''
    Hopefully self-explanatory
    :return: genome: sequence and a dictionary with positions of genes
    '''
    genome_file = open("sequenceOfEcoliStrainM54.txt", "r")
    train_file = open("train.seqs", "r")
    genome_split = genome_file.read().split()
    genome = ''
    for i in genome_split:
        genome += i
    global_coordinates = []
    # I shall have a list of dictionaries which will specify the start and end of each coding parts
    for segment in train_file:
        tmp = segment.split()
        tmp = [i.strip(',[]') for i in tmp]
        current_coordinates = {}
        current_coordinates['from'] = int(tmp[0]) - 1
        current_coordinates['to'] = int(tmp[1]) - 1
        coordinate_list = []
        temp = []
        for i in range(2, 2 + len(tmp[2:])):
            if i % 2 == 1:
                temp.append(int(tmp[i]) - 1)
                coordinate_list.append(tuple(temp))
                temp = []
            else:
                temp.append(int(tmp[i]) - 1)
        current_coordinates['positions'] = coordinate_list
        global_coordinates.append(current_coordinates)
    return genome, global_coordinates
def starting(start, global_coordinates):
    '''
    Works out the probabilities of the initial nucleaotide
    :param start: dict storing probs.
    :param global_coordinates: dict. storing segment positions
    :return: Dictionary with relative probabilities of each base
    '''
 # 1. Start
    for read in global_coordinates:
        start[gene[read['from']]] += 1
    # COUNTING PROBABILITIES
    suma = sum(start.values())
    for k, v in start.items():
        start[k] = math.log(v/suma)
    return start

def noncoding_parser(read, gene):
    '''
    Hopefully also self-explanatory
    :param read:
    :param gene:
    :return:
    '''
    # 2.1.1. extracting the very last bit after last read
    if read['to'] <= read['positions'][-1][1]:
        ending = []
    else:
        ending = [gene[read['positions'][-1][1] + 1:read['to']]]
    # 2.1.2. extracting the bit before any coding part
    res = [gene[read['from']: read['positions'][0][0]]] # extracting whatever comes before the coding part
    # 2.1.3. extracting the bits in between any coding parts
    for i in range(len(read['positions'])):
        if i == len(read['positions'])-1:
            break
        res.append(gene[read['positions'][i][1] + 1:read['positions'][i + 1][0]])  # parsing out all the bits between coding reads
    if len(ending[0])>0: #
        res.append(ending[0]) #
        ending = True#
    else:#
        ending = False#
    return [res, ending] # ending is a single letter

def noncoding(nonc, global_coordinates):
    '''
    Works out likelihoods for each base in the non-gene parts
    :param nonc: A dictionary with default values
    :param global_coordinates: A dictionary of positions of reads
    :return: A dictionary with relative base occurrences
    '''
    for read in global_coordinates:
        result = noncoding_parser(read, gene)
        for segment in result[0]:
            for j in range(len(segment)):
                if j == len(segment)-1:
                    break
                nonc[segment[j]][segment[j+1]] += 1
        if result[1]:
            nonc[result[0][-1][len(result[0][-1])- 1]]['transitions_to_end'] += 1
        # 2.2. Now, we're dealing with transitions between noncoding parts and start codon
        for coords in read['positions']:
            nonc[gene[coords[0]-1]]['transitions_to_start_codon'][gene[coords[0]]] +=1
    for n, val in nonc.items():
                midsum = 0
                for k,v in nonc[n].items():
                    if k in ['a', 't', 'c', 'g']:
                        midsum+=v
                for tn,value in nonc[n].items():
                    if tn in ['a', 't', 'c', 'g']:
                        nonc[n][tn] = math.log(nonc[n][tn]/midsum)
                        nonc[n]['transitions_to_start_codon'][tn] = math.log(nonc[n]['transitions_to_start_codon'][tn]/midsum)
                nonc[n]['transitions_to_end'] = nonc[n]['transitions_to_end']/midsum
    return nonc
##############################################
def define_genes(begin, track):
    '''
    Extracts positions of genes through backtrack
    :param begin:
    :param track:
    :return:
    '''
    fin = []
    begins = []
    adjusted_track = []
    for i in range(len(track)):
        for _ in range(3):
            adjusted_track.append(track[i])
            if track[i]>=7 or track[i]<=4:
                if 3 > adjusted_track[-1] > 1:
                    begins.append(len(adjusted_track)-1)
                break
    for initials in begins:
        for i in range(len(adjusted_track)):
            if adjusted_track[i] >5 and i >=initials and adjusted_track[i]<7:
                fin.append(i)
                break
    return [[begin + begins[i], begin + fin[i]] for i in range(len(begins))] #
##############################################
def internal_parser(read, gene):
    res = []
    ending = []
    following_bases = []
    for i in range(len(read['positions'])):
        res.append(gene[read['positions'][i][0] + 3:read['positions'][i][1] -2])  # MAKE SURE TO CHECK INDEXES HERE!!!
        if len(res[-1])%3 != 0:
            if debug_mode:
                print('we have a length of', len(res[-1]))
        ending.append(res[-1][-3:])
        following_bases.append(gene[read['positions'][i][1]-2:read['positions'][i][1]+1])
    return res, ending, following_bases  # ending here is a list of letters, be careful

def stop_parser(read, gene):
    res = []
    for i in range(len(read['positions'])):
        tmp = gene[read['positions'][i][1] - 2:read['positions'][i][1] + 1]
        if read['positions'][i][1]+1 < len(gene):
            tmp = tmp+gene[read['positions'][i][1]+1]
        else:
            if debug_mode:
                print("there is no gene left after this stop codon")
        res.append(tmp)
    return res
#######################################################

def helper_nonc(current_nucleotide, previous_nucleotide, preceeding_state, previous_codon, aggregate):
    ret = []
    res = -999999999999
    after_gene = None
    if preceeding_state == 'nongene_part':
        res = aggregate['nongene_part'][previous_nucleotide][current_nucleotide]
    elif preceeding_state == 'last_codon':
        if previous_codon in ['taa', 'tag', 'tga']:
            after_gene = True
            res = aggregate['last_codon'][previous_codon]['changes_to_noncoding'][current_nucleotide]
    ret.append(res)
    ret.append(after_gene)
    return ret

def helperinternal(each_nucleotide, codonx, previous_codon, preceeding_state, aggregate):
    res = None
    if each_nucleotide >=4 and preceeding_state == 'internals' and previous_codon not in ['tga', 'tag', 'taa']:
            if codonx is None:
                if debug_mode:
                    print('we got an empty codon')
            elif len(previous_codon) == len(codonx) == 3 :
                res = aggregate['internals'][previous_codon][codonx]
    if preceeding_state == 's3':
        if codonx is None:
            if debug_mode:
                print('we got an empty codon')
        elif len(previous_codon)==len(codonx) == 3:
            if previous_codon[-2] == 't' and previous_codon[-1] == 'g':
                res = aggregate['first_codon']['3_rd']['transition_codons'][codonx]
    if not res:
        return -999999999999
    return res
#######################################################

def dp_function(aggregate, codonx, previous_codon, T1, T2, current_sequence):
    positions = ['start', 'nongene_part', 's1', 's2', 's3', 'internals',
         'last_codon', 'end']
    d = [i for i in range(len(positions))]
    for each_nucleotide in range(len(current_sequence)):
        for av in range(len(positions)):
            if each_nucleotide == 0:
                T1[av][each_nucleotide] = 0
                continue
            T1[av][each_nucleotide] = -999999999999
            t2_ret = -len(positions)
            current_nucleotide = current_sequence[each_nucleotide]
            previous_nucleotide = current_sequence[each_nucleotide-1]
            if each_nucleotide > 2:
                previous_codon = current_sequence[each_nucleotide-3]+current_sequence[each_nucleotide-2]+current_sequence[each_nucleotide-1]
            if each_nucleotide <= len(current_sequence)-4:
                codonx = current_sequence[each_nucleotide]+current_sequence[each_nucleotide+1]+current_sequence[each_nucleotide+2]
            for p_av in range(8):
                offset = -2
                if positions[av] == 'end' or 'start' == positions[av]:
                    break
                elif positions[p_av] != 'start' and positions[av] == 'nongene_part': # chybi tu doresit previous codon
                    res = helper_nonc(current_nucleotide, previous_nucleotide, positions[p_av], previous_codon, aggregate)
                    if res[1] and res[0]+T1[p_av][each_nucleotide-3] > T1[av][each_nucleotide]:
                        T1[av][each_nucleotide] = res[0]+T1[p_av][each_nucleotide-3]
                        t2_ret = p_av
                    elif not res[1] and res[0]+T1[p_av][each_nucleotide - 1]>T1[av][each_nucleotide]:
                        T1[av][each_nucleotide] = res[0]+T1[p_av][each_nucleotide - 1]
                        t2_ret = p_av
                elif positions[av] == 's1':
                    if each_nucleotide+2 <len(current_sequence): # i.e. we cannot be at a start of a gene, if we dont have the entire start codon
                        if current_sequence[each_nucleotide+1] != 't' or current_sequence[each_nucleotide+2] != 'g':
                            break
                    if positions[p_av] == 'nongene_part':
                        tmp = T1[p_av][each_nucleotide-1]+aggregate['nongene_part'][previous_nucleotide]['transitions_to_start_codon'][current_nucleotide]
                        if tmp >T1[av][each_nucleotide]:
                            T1[av][each_nucleotide] = tmp
                            t2_ret = p_av
                elif positions[av] == 's2':
                    if each_nucleotide+1 < len(current_sequence):
                        if current_sequence[each_nucleotide+1] != 'g':
                            break
                    if current_nucleotide == 't' and positions[p_av] == 's1' and T1[p_av][each_nucleotide-1]>T1[av][each_nucleotide]:
                        T1[av][each_nucleotide] = T1[p_av][each_nucleotide-1]
                        t2_ret = p_av
                elif positions[av] == 's3':
                    if previous_nucleotide != 't':
                        break
                    if positions[p_av] == 's2' and current_nucleotide == 'g' and T1[p_av][each_nucleotide-1]>T1[av][each_nucleotide]:
                        T1[av][each_nucleotide] = T1[p_av][each_nucleotide-1]
                        t2_ret = p_av
                elif positions[av] == 'internals':
                    if positions[p_av] == 's3':
                        offset = 0
                    if 2<each_nucleotide:
                        tmp = T1[p_av][each_nucleotide+offset-1]+helperinternal( each_nucleotide, codonx, previous_codon, positions[p_av], aggregate)
                        if tmp>T1[av][each_nucleotide]:
                            T1[av][each_nucleotide] = tmp
                            t2_ret = p_av
                elif positions[av] == 'last_codon':
                    if codonx in ['taa', 'tag', 'tga'] and positions[p_av] == 'internals' and each_nucleotide >=3:
                        tmp = T1[p_av][each_nucleotide-3]+aggregate['internals'][previous_codon][codonx]
                        if tmp>T1[av][each_nucleotide]:
                            T1[av][each_nucleotide] = tmp
                            t2_ret = p_av
            T2[av][each_nucleotide] = max(-1, t2_ret)
    return T2

def startcodon(starts, global_coordinates):
    internals = {}
    starts['3_rd']['transition_codons']  = {}
    for first in ['t', 'g', 'c', 'a']:
        for second in ['t', 'g', 'c', 'a']:
            for third in ['t', 'g', 'c', 'a']:
                new_codon = first+second+third
                internals[new_codon] = {}
    for k,v in internals.items():
        for first in ['t', 'g', 'c', 'a']:
            for second in ['t', 'g', 'c', 'a']:
                for third in ['t', 'g', 'c', 'a']:
                    new_codon = first + second + third
                    temp_v = v
                    temp_v[new_codon] = 1
                    internals[k] = temp_v
        temp_v = v
        temp_v['transitions_to_end_codon'] = {'taa': 1, 'tag': 1, 'tga': 1}
        internals[k] = temp_v
        starts['3_rd']['transition_codons'][k] = 1

    for read in global_coordinates:
        for coords in read['positions']:
            starts['3_rd']['transition_codons'][gene[coords[0]+3:coords[0]+6]] +=1
            starts['3_rd']['transition_count']+=1
    suma = sum(starts['3_rd']['transition_codons'].values())
    for k, v in starts['3_rd']['transition_codons'].items():
        starts['3_rd']['transition_codons'][k] = math.log(v / suma)

    return starts, internals

def stopcodon(stopcodon, global_coordinates):
    for read in global_coordinates:
        stop_bases = stop_parser(read, gene)
        for i in range(len(stop_bases)):
            actual_stop_base = stop_bases[i][:3]
            stopcodon[actual_stop_base]['value'] += 1
            stopcodon['value'] += 1
            stopcodon[actual_stop_base]['changes_to_noncoding'][stop_bases[i][-1]] += 1
    # Stop codons
    for i in stopcodon:
        if i in ['taa', 'tag', 'tga']:
            new_sum = sum(stopcodon[i]['changes_to_noncoding'].values())
            for x in stopcodon[i]['changes_to_noncoding']:
                if x in ['a', 't', 'g', 'c']:
                    stopcodon[i]['changes_to_noncoding'][x] = math.log(stopcodon[i]['changes_to_noncoding'][x]/new_sum)
    return stopcodon
def internalcodon(internals, global_coordinates):
    permutation_list = []
    for first in ['a', 'c', 't', 'g']:
        for second in ['a', 'c', 't', 'g']:
            for third in ['a', 'c', 't', 'g']:
                tmp = first
                tmp = tmp+second
                tmp = tmp+third
                if tmp not in permutation_list:
                    permutation_list.append(tmp)
    for read in global_coordinates:
        internal_bases, ending_letters, following_bases = internal_parser(read, gene)
        for i in range(len(internal_bases)):
            for x in range(0, -5+len(internal_bases[i])):
                if x%3 == 0:
                    internals[internal_bases[i][x:x+3]][internal_bases[i][x+3:x+6]] += 1
        for z in range(len(ending_letters)):
            if following_bases[z] in ['tga', 'tag', 'taa']:
                internals[ending_letters[z]]['transitions_to_end_codon'][following_bases[z]] += 1
                internals[ending_letters[z]][following_bases[z]] += 1
    # Internal codons
    for i in internals:
        midsum = 0
        if midsum == 0:
            for k, v in internals[i].items():
                if len(k) == 3:
                    midsum +=v
        for tc, tv in internals[i].items():
            if tc in permutation_list:
                if internals[i][tc] == 0:
                    print('this case should not happen')
                    return # should not happen
                internals[i][tc] = math.log(internals[i][tc] / midsum)  ## question, should we also account for transitions to end codons
            elif tc == 'transitions_to_end_codon':  # transition to stop codon
                transition_sum = sum(internals[i][tc].values())
                for tr in ['taa', 'tag', 'tga']:
                    internals[i][tc][tr] = math.log(
                        internals[i][tc][tr] / transition_sum)
    return internals

def compute(aggregate, test_segments, gene):
    '''
    Most of the Viterbi algorithm, as per pseudocode
    :param aggregate:
    :param test_segments:
    :param gene:
    :return:
    '''
    results = []
    for segment in test_segments:
        current_sequence = gene[segment['from']:segment['to']]  # this is the sample we are going to be testing
        #######################################################################
        T1 = [[0 for _ in current_sequence] for _ in ['start', 'nongene_part', 's1', 's2', 's3',
            'internal_codons', 'end_codon', 'end']]
        switch = 1
        cesta = []
        T2 = [[-1 for _ in current_sequence] for _ in T1]
        for i in range(8):
            if i == 1:
                T2[i][0] = 0
                T1[i][0] = aggregate['start'][current_sequence[0]]
            else:
                T1[i][0] = -999999999999
                T2[i][0] = -1
        if debug_mode:
            print("ag keys", aggregate.keys())
        midresult = dp_function(aggregate, None, None, T1, T2, current_sequence)
        vzdalenost = 0
        for i in range(len(current_sequence) - switch, 0, -1):
            tmp = [midresult[switch][i - vzdalenost]]
            tmp.extend(b for b in cesta)
            cesta = tmp
            switch = midresult[switch][i - vzdalenost]
            if switch == 0:
                break
            if 4 < switch < 7:
                vzdalenost += 2
        results.append(define_genes(segment['from'], cesta))
    return results

if __name__ == "__main__":
    write_up = True
    debug_mode = True
    start = {} # start base
    nonc = {} # noncoding segment
    starts = {} # start codon
    for i in ['a', 't', 'g', 'c']:
        start[i] = 1
    for i in start.keys():
        nonc[i] = {}
        for x in start.keys():
            nonc[i][x] = 0
        nonc[i]['transitions_to_end'] = 1
        nonc[i]['transitions_to_start_codon'] = {'a': 1, 't': 1, 'g': 1, 'c': 1}
    # 1. parse out gene and coordinates of individual reads
    gene, global_coordinates = parse_input() # gene is training sequence, global coordinates are corresponding
    # Work out probabilities of start section
    start = starting(start, global_coordinates)
     # 2. Noncoding
    nonc = noncoding(nonc, global_coordinates)
    # 3. Start codon
    starts['1_st']={}
    starts['2_nd'] = {}
    for k,v in starts.items():
        starts[k]['prob_to_next'] = 1
    starts['3_rd'] = {'transition_count': 0, 'transition_codons': {'a': 1, 't': 1, 'g': 1, 'c': 1}}
    starts, internals = startcodon(starts, global_coordinates)
    # 4. Internal codons ###################
    internals = internalcodon(internals, global_coordinates)
    # 5. STOP codons ############
    stop_codon = {'value':0}
    for combo in ['taa', 'tag', 'tga']:
        stop_codon[combo] = {'value':0}
    for k,v in stop_codon.items():
        if len(k) == 3:
            stop_codon[k]['changes_to_noncoding'] = {'a': 1, 't': 1, 'g': 1, 'c': 1}
    stop_codon = stopcodon(stop_codon, global_coordinates)
    # Retrieving test cases
    test_segments, size = parse_test()
    aggregate = {'last_codon':stop_codon, 'first_codon':starts, 'internals':internals, 'nongene_part':nonc,'start':start }
    results = compute(aggregate, test_segments, gene)
    if write_up:
            predictions_file = open("predictions.txt", "w")
            one = 0
            three = 0
            two = {}
            four = 0 # overlap

            for c in range(len(results)):
                predictions_file.write(str(test_segments[c]['from']+1)+" "+str(test_segments[c]['to']+1))
                for i in results[c]:
                    if i == []:
                        break
                    one+=1
                    tmp = gene[i[0]:i[1]+3]
                    three+=len(tmp)
                    broken = False
                    for ix in test_segments:
                        for p in ix['positions']:
                            p = (p[0]+1, p[1]+1)
                            print('p is:', p)
                            new_min = max(p[0], i[0])
                            new_max = min(p[1], i[1]+2)
                            if new_min<=new_max:
                                print('od', new_min, 'do', new_max, 'adding', 1+new_max-new_min)
                                two[i[0]] = 1
                                four+=new_max-new_min
                            if p[0]<= i[0] and p[1]>=i[0]:
                                broken = True
                                break
                            if p[0]<=i[1]+2 and p[1]>=i[1]+2:
                                broken = True
                                break
                        if broken:
                            break
                    predictions_file.write(' '+'['+str(i[0])+', '+str(i[1]+2)+']')
                predictions_file.write('\n')
            ### Comparing accuraccy
            ## two, number of genes which have nonempty overlap with one of the truth genes
            predictions_file = open("accuracy.txt", "w")
            predictions_file.write(str(one)+'\n')
            predictions_file.write(str(sum(two.values())) + '\n')
            predictions_file.write(str(three) + '\n')
            predictions_file.write(str(four) + '\n')
            predictions_file.write(str(four/(size + three-four)) + '\n')














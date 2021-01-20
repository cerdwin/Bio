#!/usr/bin/python
 
import sys
import math
 
 
def classify(sequence, cpg_matrix, null_matrix):
    column_dict = {'A':0, 'T':1, 'G':2, 'C':3, '\n':4}
    row_dict = {'A':0, 'T':1, 'G':2, 'C':3, 'start':4}
    iterator = 0;
    cpg_score = 0
    null_score = 0
    for i in range(len(sequence)):
        if i == 0:
            cpg_score += math.log(cpg_matrix[row_dict['start']][column_dict[sequence[i]]])
            null_score += math.log(null_matrix[row_dict['start']][column_dict[sequence[i]]])
        elif i != len(sequence)-1:
            cpg_score += math.log(cpg_matrix[row_dict[sequence[i]]][column_dict[sequence[i+1]]])
            null_score += math.log(null_matrix[row_dict[sequence[i]]][column_dict[sequence[i + 1]]])
    # print('cpg:', cpg_score, 'null:', null_score)
    if cpg_score > null_score:
        return 1
    return 0
 
def find_item(where, item):
    my_counter = 0
    for i in range(len(where)):
        if where[i] == item:
            my_counter += 1
    return my_counter
 
def transition_counter(v, i, my_array):
    my_counter = 0
    for s in range(len(my_array)):
        if (my_array[s] == v) and ((s+1)< len(my_array)) and (my_array[s+1] == i):
            my_counter +=1
    return  my_counter
 
cpg_train_file = open("cpg_train.txt","r")
null_train_file = open("null_train.txt","r")
 
cpg_array = []
null_array = []
cpg_size = 0 # count of items in the cpg file
null_size = 0 # count of items in the null file
# get cpg file into an array
for line in cpg_train_file:
    for letter in line:
        cpg_array.append(letter)
 
# get null file into an array
for line in null_train_file:
    for letter in line:
        null_array.append(letter)
cpg_size = len(cpg_array)
null_size = len(null_array)
# column legend: [*, A, T, G, C, end], row: [*, A, T, G, C, start]
cpg_matrix = [[0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0]]
null_matrix = [[0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0]]
column_legend = ['A', 'T', 'G', 'C', '\n']
row_legend = ['A', 'T', 'G', 'C', 'start']
iterator = 0
for i in column_legend:
    cpg_matrix[4][iterator] = (find_item(cpg_array, i)+1)/(len(cpg_array)+5)
    null_matrix[4][iterator] =(find_item(null_array, i)+1)/(len(null_array)+5)
    iterator +=1
# print(cpg_array)
for i in range(len(column_legend)):
    for v in range(len(row_legend)):
        if i != 4:
            # print('col_legend:', column_legend[v], 'row_legend:', row_legend[i], '\n')
            cpg_matrix[i][v] = (transition_counter(row_legend[i], column_legend[v],cpg_array)+1)/(len(cpg_array)+5)
            null_matrix[i][v] = (transition_counter(row_legend[i], column_legend[v], null_array)+1)/(len(null_array)+5)
 
 
# print(cpg_matrix)
 
# maybe we should build the classifier here
test_array = []
results = []
test_sequences_file = open("seqs_test.txt","r")
for line in test_sequences_file:
    for letter in line:
        test_array.append(letter)
    results.append(classify(test_array, cpg_matrix, null_matrix))
    test_array = []
# print(results)
# probably we should test the sequences to get predictions ...
 
# ... and store the predictions here
predictions_file = open("predictions.txt","w")
for i in results:
    predictions_file.write(str(i))
    predictions_file.write('\n')
# ... and finally, we should compare the predictions to the ground truth
correct = 0
wrong = 0
 
null_train_file = open("classes_test.txt","r")
iterator = 0
trues = 0
falses = 0
TP = 0
TN = 0
FP = 0
FN = 0
P = 0
N=0
for line in null_train_file:
    if line[0] == '1':
        P+=1
    if line[0] == '0':
        N+=1
    if line[0] == '1' and results[iterator] == 1:
        correct += 1
        trues +=1
        TP +=1
    elif line[0] == '0' and results[iterator] == 0:
        correct+=1
        falses+=1
        TN+=1
    elif line[0] == '0'and results[iterator] == 1:
        FP +=1
        wrong+=1
    else:
        FN +=1
        wrong+=1
    iterator +=1
 
# print('correct:', correct, 'wrong:', wrong)
accuracy1 = correct/(correct+wrong)
# 1. Number of correct predictions: TN+TP
correct = TN+TP
# 2. Number of wrong predictions: FN+FP
wrong = FN+FP
# 3. accuracy: (TP+TN)/(P+N)
accuracy = (TP+TN)/(P+N)
# 4. precision: TP/(TP+FP)
precision = TP/(TP+FP)
# 5. recall: TP/P
recall = TP/P
 
# ... and store the accuracy, precision, and recall
accuracy_file = open("accuracy.txt","w")
accuracy_file.write(str(correct))
accuracy_file.write('\n')
 
accuracy_file.write(str(wrong))
accuracy_file.write('\n')
 
accuracy_file.write(str(accuracy))
accuracy_file.write('\n')
 
accuracy_file.write(str(precision))
accuracy_file.write('\n')
 
accuracy_file.write(str(recall))
accuracy_file.write('\n')

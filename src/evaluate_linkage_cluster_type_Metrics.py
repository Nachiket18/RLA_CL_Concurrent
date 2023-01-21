import math

import numpy as np
import pandas as pd
import csv
import random
from math import comb


# Reads a cluster file in Mamun's format
# Returns a vector of clusters where clusters are vector of SSN in the cluster
def get_rlaCL_clusters(path):
    total_cluster_vec = []
    cluster_vec = []
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        for row in csv_reader:
            if len(row) == 1:
                cluster_vec = []
            if len(row) > 1:
                cur_ssn = row[0]
                cluster_vec.append(cur_ssn)
            if len(row) < 1:
                if len(cluster_vec) != 0:
                    cluster_vec.sort()
                    total_cluster_vec.append(cluster_vec)
                    cluster_vec = []
    # ruru
    # for i in total_cluster_vec:
    #     print(i)
    return total_cluster_vec

def get_my_clusters(path):
    total_cluster_vec = []
    cluster_vec = []
    with open(path, "r", encoding="utf8") as file:
        csv_reader = csv.reader(file, delimiter=",")
        for row in csv_reader:
            cluster_vec = []
            if len(row) > 0:
                for uid in row:
                    if len(uid) > 1:
                        cluster_vec.append(uid)
            cluster_vec.sort()
            total_cluster_vec.append(cluster_vec)

    #ruru
    # for i in total_cluster_vec:
    #     print(i)
    return total_cluster_vec

# Takes input vector of clusters
# Outputs (i) a dictionary where SSN are keys and a tuple nis the value. Tuple contains number of records of key
# found together in clusters
# Outputs (ii) a dictionary where SSN are keys and value is the total count of the SSN found in all clusters combined
# Outputs (iii) a vector of dictionary where each dictionary contains which ssn in contained
# in the corresponding cluster how many times
def get_ssn_info(cluster_vec):
    ssn_group_info = {}
    ssn_counts = {}
    ssn_groups_in_single_cluster = []
    for x in cluster_vec:
        cur_cluster_ssn_group = {}
        for y in x:
            if y in ssn_counts:
                ssn_counts[y] = ssn_counts[y] + 1
            else:
                ssn_counts[y] = 1
            if y in cur_cluster_ssn_group:
                cur_cluster_ssn_group[y] = cur_cluster_ssn_group[y] + 1
            else:
                cur_cluster_ssn_group[y] = 1
        for c in cur_cluster_ssn_group:
            if c in ssn_group_info:
                ssn_group_info[c].append(cur_cluster_ssn_group[c])
            else:
                ssn_group_info[c] = []
                ssn_group_info[c].append(cur_cluster_ssn_group[c])

        ssn_groups_in_single_cluster.append(cur_cluster_ssn_group)
    # for x in ssn_group_info:
    #     print(x)
    #     print(ssn_counts[x])
    #     print(ssn_group_info[x])
    return ssn_counts, ssn_group_info, ssn_groups_in_single_cluster


# Inputs: vector of dictionary {ssn:count} corresponding to a cluster, dictionary containing total number of each SSN
# Output: array[4] containing count of clusters of type i in i-1 th index
def get_cluster_types(ssn_groups_in_single_cluster, ssn_counts):
    t_counts = [0, 0, 0, 0]
    for my_cluster_dict in ssn_groups_in_single_cluster:
        is_type_three = False
        if len(my_cluster_dict) == 1:
            key = next(iter(my_cluster_dict))
            if my_cluster_dict[key] == ssn_counts[key]:
                t_counts[0] = t_counts[0] + 1  # Type 1: Cluster contains ALL records of only one SSN
            else:
                t_counts[1] = t_counts[1] + 1  # Type 2: Cluster contains SOME records of only one SSN
        else:
            for x in my_cluster_dict:
                if my_cluster_dict[x] == ssn_counts[x]:
                    is_type_three = True  # Type 3: Cluster contains ALL records of an SSN and some of other
                    break
            if is_type_three:
                t_counts[2] = t_counts[2] + 1
            else:
                t_counts[3] = t_counts[3] + 1  # Type 4: Cluster contains some records of multiple SSNs
    return t_counts


# Input: dictionary of vectors of numbers of same ssn records found together in clusters
# Output: Total True Positive Links
def get_linkage_tp(ssn_groups):
    tp_count = 0
    for x in ssn_groups:
        for y in ssn_groups[x]:
            tp_count = tp_count + math.comb(y, 2)
    return tp_count


# Input: vector of dictionary {ssn:count} corresponding to a cluster
# Output: Total False Positive Links
def get_linkage_fp(cluster_ssn_groups):
    fp_count = 0
    for x in cluster_ssn_groups:
        c_records = 0
        for y in x:
            c_records = c_records + x[y]
        fp_in_x = 0
        for y in x:
            fp_in_x = fp_in_x + (x[y] * (c_records - x[y]))
        fp_in_x = fp_in_x / 2
        # print(f"Cluster SSNs: {x}, total records in this cluster: {c_records}, fp links: {fp_in_x}")
        fp_count = fp_count + fp_in_x
    return int(fp_count)


# Input: dictionary {ssn:Totalcount}, dictionary (ssn: vector[c1,c2,...,cK]) corresponding to numbers of records of key
#       in clusters 1,2,...k
# Output: Total True Negative Links
def get_linkage_fn(ssn_dict, ssn_group_dict):
    fn_count = 0
    for x in ssn_group_dict:
        cur_ssn_fn = math.comb(ssn_dict[x], 2)  # Max possible links (upper bound)
        for y in ssn_group_dict[x]:
            cur_ssn_fn = cur_ssn_fn - math.comb(y, 2)  # Cancelling True positives
        fn_count = fn_count + cur_ssn_fn
    return fn_count


# Returns total number of records
def get_total_records(ssn_dict):
    r_count = 0
    for x in ssn_dict:
        r_count = r_count + ssn_dict[x]
    return r_count


# Input: vector of dictionary {ssn:count} corresponding to a cluster
# Output: Dict{size:count} of cluster sizes and number of them
def get_cluster_sizes(clusters):
    c_size_dict = {}
    for x in clusters:
        sum = 0
        for y in x:
            sum = sum + x[y]
        if sum in c_size_dict:
            c_size_dict[sum] = c_size_dict[sum] + 1
        else:
            c_size_dict[sum] = 1
    return c_size_dict

### main ###

# file_path = "/Users/joyanta/Downloads/output_edit_ds2.1.txtOutSingle"
file_path = "/home/nachiket/RLA_CL_parallel/output_parallel_rla"

### Calculate Cluster Accuracy

cluster_members_vec = get_my_clusters(file_path)
ssn_tot_dict, ssn_group_dict, per_cluster_ssn_group_dict = get_ssn_info(cluster_members_vec)
total_clusters = len(cluster_members_vec)
types_counts = get_cluster_types(per_cluster_ssn_group_dict, ssn_tot_dict)

print(f"Total Clusters: {total_clusters}")
for i in range(len(types_counts)):
    print(f"Type {i + 1} Count: {types_counts[i]} Percentage: {(types_counts[i] / total_clusters) * 100}")

total_records = get_total_records(ssn_tot_dict)
print(f"Total records: {total_records} of {len(ssn_tot_dict)} persons")
linkage_TP = get_linkage_tp(ssn_group_dict)
print(f"True Positive Links: {linkage_TP}")

linkage_FP = get_linkage_fp(per_cluster_ssn_group_dict)
print(f"False Positive Links: {linkage_FP}")

linkage_FN = get_linkage_fn(ssn_tot_dict, ssn_group_dict)
print(f"False Negative Links: {linkage_FN}")

linkage_TN = math.comb(total_records, 2) - linkage_TP - linkage_FP - linkage_FN
print(f"True Negative Links: {linkage_TN}")

# P = tp / (tp+fp)
precision = float(linkage_TP) / float(linkage_TP + linkage_FP)
# R = tp / (tp+fn)
recall = float(linkage_TP) / float(linkage_TP + linkage_FN)
# A = (tp + tn) / (tp + tn + fp + fn)
accuracy = float(linkage_TP + linkage_TN) / float(linkage_TP + linkage_TN + linkage_FN + linkage_FP)
# f1 = 2tp / (2tp + fp + fn)
f1_score = float(2 * linkage_TP) / float(2 * linkage_TP + linkage_FN + linkage_FP)

print(f"Linkage Precision: {precision}")
print(f"Linkage Recall: {recall}")
print(f"Linkage f1 score: {f1_score}")
print(f"Linkage Accuracy: {accuracy}")
cluster_sizes_dict = get_cluster_sizes(per_cluster_ssn_group_dict)
print(cluster_sizes_dict)
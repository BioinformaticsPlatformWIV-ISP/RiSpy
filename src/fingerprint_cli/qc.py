import pandas as pd
import yaml
from Bio import SeqIO
import csv
from collections import Counter


def calculate_mean_and_median_depth(input: str, output: str, config: str) -> None:
    """
    Calculating the mean and median depth of coverage from a depth file
    """
    total_depth = 0
    total_count = 0
    depth_counts = Counter()
    chunk_size = 100  # increase for speed if enough memory
    for chunk in pd.read_table(input, header=0, names=['chr', 'pos', 'depth'], chunksize=chunk_size):
        total_depth += chunk['depth'].sum()
        total_count += chunk['depth'].count()
        # count occurrences of each depth value
        depth_counts.update(chunk['depth'].astype(int).tolist())
    mean = total_depth / total_count
    # Compute median from histogram
    cum_count = 0
    median_idx = (total_count + 1) // 2
    for depth in sorted(depth_counts):
        cum_count += depth_counts[depth]
        if cum_count >= median_idx:
            median = depth
            break
    with open(output, 'w') as handle:
        handle.write(f"Mean depth: {mean}\tMedian depth: {median}\n")

    results = {
        "Mean Depth of Coverage": str(round(mean, 2)),
        "Median Depth of Coverage": str(median)
    }
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)

def calculate_ref_length(ref: str) -> int:
    """
    Calculating the length of the reference genome
    """
    reference =  SeqIO.parse(ref, "fasta")
    ref_len = 0
    for sequence in reference:
        ref_len += len(sequence.seq)
    return ref_len

def calculate_breadth_of_coverage(input: str, output: str, cutoff: int, ref_length: int, config: str) -> None:
    i = 0
    with open(input, 'r') as handle:
        csv_reader = csv.reader(handle, delimiter='\t')
        #skip first row (=header)
        next(csv_reader)
        for row in csv_reader:
            #Only count site if depth above a certain cutoff
            if int(row[2]) > cutoff:
                i += 1
    breadth = (float(i)/ref_length)*100
    with open(output, 'w') as handle:
        handle.write(str(breadth))
        handle.write('\n')
    results = {"Breadth of Coverage":str(round(breadth, 2))}
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)


def calculate_fold80(cov_hist: str, depth: str, count_file: str, ref_length: int, output: str, config: str) -> None:
    total_depth = 0
    total_count = 0
    chunk_size = 100  # increase for speed if enough memory
    for chunk in pd.read_table(depth, header=0, names=['chr', 'pos', 'depth'], chunksize=chunk_size):
        total_depth += chunk['depth'].sum()
        total_count += chunk['depth'].count()
    mean_coverage = total_depth / total_count
    percentile_20 = ref_length * 0.2
    cov = pd.read_csv(cov_hist,sep="\t")
    cov['#Coverage'] = cov['#Coverage'].astype(int)
    cov['Number of genomic locations'] = cov['Number of genomic locations'].astype(int)
    #Set 'Coverage' as the index
    cov.set_index('#Coverage',inplace=True)
    #Create a complete range of coverage values
    full_range = pd.RangeIndex(start=0,stop=cov.index.max() + 1,step=1)
    #Reindex the DataFrame to include all coverage values, filling missing entries with 0
    cov = cov.reindex(full_range,fill_value=0)
    cov.reset_index(inplace=True)
    cov.rename(columns={'index': 'Coverage'},inplace=True)
    # get the list of coverage values. As we do not change the order, the index is the coverage
    # e.g. [7, 100, 500, 100] -> 7 bases 0X, 100 bases 1X, 500 bases 2X, 100 bases 3X
    cov_values = list(cov['Number of genomic locations'])
    # calculate cumsum for the histogram:
    # e.g. [7, 100, 500, 100] -> [7, 107, 607, 707]: 7 bases 0X, 107 bases max. 1X, 607 bases max. 2X, 707 bases max. 3X
    histo_cumsum = [sum(cov_values[:i + 1]) for i, value in enumerate(cov_values)]
    #add the input percentile_20 to the cumsum list and sort\
    tmp = histo_cumsum + [int(percentile_20)]
    tmp.sort()
    # Find the coverage of the 20th percentile
    # e.g. 20% 707 = 141 -> [7, 107, 141, 607, 707]
    find_index = tmp.index(int(percentile_20))
    val_before = histo_cumsum[find_index - 1]
    val_after = histo_cumsum[find_index]
    #in our example, find_index will be 2, val_before: 107, val_after: 607
    cov_20 = (((val_after - percentile_20) / (val_after - val_before)) + find_index - 1)
    #ideally, as close as possible to 1
    fold80 = mean_coverage / cov_20
    with open(output, 'w') as handle:
        handle.write(str(fold80))
        handle.write('\n')
    with open(count_file,"r") as c:
        count = int(c.read().strip())
    results = {"Number of Primary Mapped Reads":str(count), "Fold80":str(round(fold80, 2))}
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)

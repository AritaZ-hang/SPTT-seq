import argparse
import pysam
import pandas as pd
from tqdm import tqdm
import pickle

def load_bcdict(dict_path):
    if dict_path is not None:
        dict = pickle.load(open(dict_path, "rb"))
        return dict
    else:
        return None

def mapping_list_preprocessing(mapping_dict):
    if mapping_dict is not None:
        all_mapped_bead = list(mapping_dict.keys())
    else:
        all_mapped_bead = None
    return all_mapped_bead

def args_parse():
    parser = argparse.ArgumentParser(description="Summarise UTI tag. Beads barcodes is the 1-12nt in R1, UMI tag is either 13-22nt in R1 or 1-10nt in R2 based on the batch of experiment.")
    parser.add_argument("-in_fq_1", type = str, default = None, help = "The filename of input fastq1.")
    parser.add_argument("-in_fq_2", type = str, default=None, help = "The filename of input fastq2.")
    parser.add_argument("-workdir", type = str, default=None, help = "The working directory.")
    parser.add_argument("-mapping_dic", type = str, default=None, help = "The raw barcodes-degen barcodes mapping dic.")
    parser.add_argument("-umi_loc", type = str, default=None, help = "Whether the umi is located in R1 or R2. possible inputs are: r1, r2.")
    return parser.parse_args()

def writeDictToCSV(dict_to_summarize):
    target_csv = pd.DataFrame.from_dict(dict_to_summarize, orient = "index")
    target_csv = target_csv.fillna(0)
    return target_csv

def uti_summarise(in_fq_1: str, 
                  in_fq_2: str, 
                  workdir: str, 
                  mapping_dic:str, 
                  umi_loc:str):
    
    degen_flag = False
    all_mapped_bead = mapping_list_preprocessing(mapping_dic)
    
    # Whether to degen the bead barcodes
    mapping_dic = load_bcdict(mapping_dic)
    if mapping_dic is not None:
        degen_flag = True
    
    beads_umi_dic = dict()
    with pysam.FastxFile(in_fq_1) as fastq_file1, pysam.FastxFile(in_fq_2) as fastq_file2:
        for read1, read2 in tqdm(zip(fastq_file1, fastq_file2)):
            fastq_read1 = {"name": read1.name, 
                           "sequence": read1.sequence, 
                           "comment": read1.comment if read1.comment else "", 
                           "quality": read1.quality if read1.quality else ""}
            fastq_read2 = {"name": read2.name, 
                           "sequence": read2.sequence, 
                           "comment": read2.comment if read2.comment else "", 
                           "quality": read2.quality if read2.quality else ""}
            barcodes = fastq_read1["sequence"][0:12]
            if umi_loc == "r1":
                umi = fastq_read1["sequence"][12:22]
            elif umi_loc == "r2":
                umi = fastq_read2["sequence"][0:10]

            if degen_flag == True:
                if barcodes in all_mapped_bead:
                    barcodes = mapping_dic[barcodes]
            
            if barcodes not in beads_umi_dic.keys():
                beads_umi_dic[barcodes] = {}
            if umi not in beads_umi_dic[barcodes].keys():
                beads_umi_dic[barcodes][umi] = 1
            else:
                beads_umi_dic[barcodes][umi] += 1

    with open(workdir + "/beads_umi_dic.pkl", "wb") as f:
        pickle.dump(beads_umi_dic, f)

    utis = {}
    for group in beads_umi_dic.keys():
        if group not in utis.keys():
            utis[group] = len(beads_umi_dic[group].keys())
    
    with open(workdir + "/utis.pkl", "wb") as f:
        pickle.dump(utis, f)

    utis_df = writeDictToCSV(utis)
    utis_df.columns = ["reads"]
    utis_df.to_csv(workdir + "./utis.csv")
            
    return

def main():
    args = args_parse()
    uti_summarise(in_fq_1 = args.in_fq_1, 
                  in_fq_2 = args.in_fq_2, 
                  workdir = args.workdir, 
                  mapping_dic = args.mapping_dic, 
                  umi_loc = args.umi_loc)
    return()

if __name__ == "__main__":
    main()
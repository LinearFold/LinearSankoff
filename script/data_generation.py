from asyncio import selector_events
from collections import defaultdict
from genericpath import isdir
from turtle import st
from unittest import defaultTestLoader
from Bio import SeqIO
import subprocess
import random
import os
import sys

random.seed(1)

def cal_seq_iden(seq1, seq2):
    assert len(seq1) == len(seq2)

    seq1 = seq1.upper().replace('T', "U")
    seq2 = seq2.upper().replace('T', "U")

    total = 0
    num_same = 0
    for nuc1, nuc2 in zip(seq1, seq2):
        if nuc1 in "-N" and nuc1 == nuc2:
            continue
        total += 1
        if nuc1 == nuc2: 
            num_same += 1

    return num_same / total


def sample(data_dir, file, fam, outdir, samplesize):
    data = {}
    alndata = {}
    for item in SeqIO.parse(file, "fasta"):
        # get seq from *.seq
        name = item.name
        alnseq = str(item.seq)

        dotfile = os.path.join(data_dir, name + ".dotbracket")
        if not os.path.exists(dotfile): continue
        with open(dotfile, "r") as f:
            lines = f.readlines()
        struc = lines[-1].strip()
        unstruc = struc.replace("(", "").replace(")", "").replace("<", "").replace(">", "")
        # print(unstruc)
        ratio = len(unstruc)/len(struc)
        if ratio > 0.5:
            # print(name, ratio, struc)
            continue
        
        seqfile = os.path.join(data_dir, name + ".seq")
        with open(seqfile, "r") as f:
            lines = f.readlines()
        seq = lines[-1].strip()[:-1].replace("T", "U")

        if len(alnseq.replace("-", "")) != len(seq): continue

        if "N" in seq: continue
        if "SRP" in fam:
            if len(seq) < 200: continue
        elif "16S" in fam:
            if len(seq) < 1400: continue

        if "tmRNA" in fam:
            # tmRNA at least 4 branch
            # get ct file from *.seq
            ctfile = os.path.join(data_dir, name + ".ct")
            # ct to dot-bracket format
            dbfile = os.path.join(data_dir, "%s.dotbracket" % name)

            if not os.path.exists(ctfile): continue
            if not os.path.exists(dbfile): # should exist
                print(ctfile)
                subprocess.run(["/mnt/home/sizhen/benchmark/modified/RNAstructure/exe/ct2dot", ctfile, "1",
                                dbfile])

            with open(dbfile, "r") as dbf:
                struc = dbf.readlines()[-1].strip()
                # print(struc)
                stack = []
                num = 0
                for i, n in enumerate(struc):
                    if n in ".()": continue
                    if n == "<":
                        stack.append(i)
                        continue
                    if n == ">":
                        stack.pop()
                        if len(stack) == 0:
                            num += 1

                # print(num)
                if num != 4: 
                    continue

        data[name] = seq

        newalnseq = ""
        pos = 0
        for i in range(len(alnseq)):
            if alnseq[i] == "-":
                newalnseq += "-"
                continue
            if pos >= len(seq):
                print(name, alnseq)
                print(seq)
            newalnseq += seq[pos]
            pos += 1
        alndata[name] = newalnseq

        unknown_chars = seq.replace("-", "").replace("A", "").replace("U", "").replace("G", "").replace("C", "")
        if len(unknown_chars) > 0:
            print(fam, name, unknown_chars)

    # calculate pairwise seq identity and order by it
    print("pairwise sequence calculation...")
    pairwise_seqs = defaultdict(set)
    all_names = list(data.keys())
    num_seqs = len(all_names)
    print("number of seqs: ", num_seqs)
    for i in range(num_seqs):
        name1 = all_names[i]
        for j in range(i+1, num_seqs):
            name2 = all_names[j]
            seq_idntty = cal_seq_iden(alndata[name1], alndata[name2])
            if seq_idntty == 1.0:
                continue

            if seq_idntty >= 0.9:
                pairwise_seqs[0.9].add((name1, name2))

            elif seq_idntty >= 0.8:
                pairwise_seqs[0.8].add((name1, name2))

            elif seq_idntty >= 0.7:
                pairwise_seqs[0.7].add((name1, name2))

            elif seq_idntty >= 0.6:
                pairwise_seqs[0.6].add((name1, name2))

            elif seq_idntty >= 0.5:
                pairwise_seqs[0.5].add((name1, name2))

            elif seq_idntty >= 0.4:
                pairwise_seqs[0.4].add((name1, name2))

            elif seq_idntty >= 0.3:
                pairwise_seqs[0.3].add((name1, name2))

            else:
                pairwise_seqs[0].add((name1, name2))
                # print(name1, name2, seq_idntty)
                # print(alndata[name1])
                # print(alndata[name2])

    count = 1
    for key, items in pairwise_seqs.items():
        print(key, len(items))
        namepairs = random.choices(list(items), k=samplesize)
        for (name1, name2) in namepairs:
            seq_idntty = cal_seq_iden(alndata[name1], alndata[name2])
            with open(os.path.join(outdir, "%s_%f.fasta" % (fam, seq_idntty)), "w") as f:
                f.write(">%s %d\n" % (name1, len(data[name1])))
                f.write("%s\n" % data[name1]) # 
                f.write(">%s %d\n" % (name2, len(data[name2])))
                f.write("%s\n" % data[name2]) # 
            count += 1

    # # return;
    # for i in range(1, 21):
    #     names = random.choices(list(data.keys()), k=2)
    #     # homologous seqs
    #     with open(os.path.join(outdir, "%s_group%d.fasta" % (fam, i)), "w") as f:
    #         for name in names:
    #             # print(data[name].replace("-", "").replace("A", "").replace("U", "").replace("G", "").replace("C", ""))
    #             if len(data[name].replace("-", "").replace("A", "").replace("U", "").replace("G", "").replace("C", "")) > 0:
    #                 print(outdir, name,
    #                     data[name].replace("-", "").replace("A", "").replace("U", "").replace("G", "").replace("C", ""))
    #             f.write(">%s %d\n" % (name, len(data[name].replace("-", ""))))
    #             f.write("%s\n" % data[name].replace("-", "")) # 


def train_data_generation():
    # homologs seqs
    # four family, each family sample 10 pairs
    # 5S ribosomal RNA (Eubacteria subfamily), group I intron (IC1 subfamily), tmRNA, and tRNA families
    data_dir = "/scratch/lisiz/dataset/RNAStrAlign"
    out_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/train_data_new"

    fam_dirs = os.listdir(data_dir)
    for fam_dir in fam_dirs:
        dir_path = os.path.join(data_dir, fam_dir)
        if not os.path.isdir(dir_path):
            continue

        if "group" not in fam_dir: continue

        if "tRNA" in fam_dir or "tmRNA" in fam_dir or "5S" in fam_dir or "group" in fam_dir:
            out_dir_fam = os.path.join(out_dir, fam_dir)
            if not os.path.exists(out_dir_fam):
                os.mkdir(out_dir_fam)

            fam_name =  fam_dir.replace("_database", "")
            if "tRNA" in fam_dir or "tmRNA" in fam_dir: 
                fasta_file = os.path.join(dir_path, "%s.fasta" % fam_name)
                sample(os.path.join(data_dir, fam_dir), fasta_file, fam_name, out_dir_fam, 2)
            else:
                subfam_dirs = os.listdir(dir_path)
                for subfam_dir in subfam_dirs:
                    if subfam_dir == "CVS": continue
                    if subfam_dir.startswith("."): continue

                    out_dir_sumfam = os.path.join(out_dir_fam, subfam_dir)
                    if not os.path.exists(out_dir_sumfam):
                        os.mkdir(out_dir_sumfam)
                    
                    # subfam_name = subfam_dir
                    fasta_file = os.path.join(dir_path, subfam_dir, "%s.fasta" % subfam_dir)
                    # print(fasta_file)
                    sample(os.path.join(data_dir, fam_dir, subfam_dir), fasta_file, "%s_%s" % (fam_name, subfam_dir), out_dir_sumfam, 1)

        
def cal_avg_seq_iden():
    aln_file = ""

    aln_data = {}
    for item in SeqIO.parse(aln_file, "fasta"):
        name = item.id
        seq = str(item.seq)
        aln_data[name] = seq
    
    sum_seq_iden = 0 
    seq_id_data = {}
    for name1 in aln_data:
        for name2 in aln_data:
            if name1 != name2: continue
            if (name1, name2) in seq_id_data or (name2, name1) in seq_id_data: continue
            seq_iden = cal_seq_iden(aln_data[name1], aln_data[name2])
            seq_id_data[(name1, name2)] = seq_iden
            sum_seq_iden += seq_iden
            

    print(aln_file)
    print(sum_seq_iden / len(aln_data))

def seq_identity(aln_file):
    aln_data = {}
    for item in SeqIO.parse(aln_file, "fasta"):
        name = item.id
        seq = str(item.seq)
        aln_data[name] = seq
    
    count = 0
    sum_seq_identity = 0
    for name1 in aln_data:
        for name2 in aln_data:
            if name1 == name2: continue
            seq1 = aln_data[name1]
            seq2 = aln_data[name2]

            aln_len = 0
            match = 0
            for nuc1, nuc2 in zip(seq1, seq2):
                if nuc1 == "-" and nuc2 == "-": continue
                if nuc1 == "N" or nuc2 == "N": continue
                aln_len += 1
                if nuc1 == nuc2:
                    match +=1
            sum_seq_identity += match / aln_len
            count += 1

    print(sum_seq_identity / count)


def seq_identityALL():
    data_dir = "/Users/sizhenli/Desktop/bioinfo/dataset/RNAStrAlign-master-aedbaf9e95ebad34d15752fe1b750f5ff3960bca"

    fam_dirs = os.listdir(data_dir)
    for fam_dir in fam_dirs:
        dir_path = os.path.join(data_dir, fam_dir)
        if not os.path.isdir(dir_path):
            continue

        fam_name =  fam_dir.replace("_database", "")
        if "tRNA" in fam_dir or "tmRNA" in fam_dir: 
            fasta_file = os.path.join(dir_path, "%s.fasta" % fam_name)
            print(fasta_file)
            # sample(fasta_file, fam_name, out_dir)
            seq_identity(fasta_file)
        elif "5S" in fam_dir or "group" in fam_dir:
            subfam_dirs = os.listdir(dir_path)
            for subfam_dir in subfam_dirs:
                if subfam_dir == "CVS": continue
                if subfam_dir.startswith("."): continue
                
                if "5S" in fam_dir:
                    if "Bacteria" not in subfam_dir: 
                        continue
                if "group" in fam_dir:
                    if "IC1" not in subfam_dir: 
                        continue
                
                subfam_name = subfam_dir
                fasta_file = os.path.join(dir_path, subfam_dir, "%s.fasta" % subfam_dir)
                print(fasta_file)
                # sample(fasta_file, "%s_%s" % (fam_name, subfam_dir), out_dir)
                seq_identity(fasta_file)
        else:
            pass


def test_data_generation():
    # homologs seqs
    # four family, each family sample 10 pairs
    # 5S ribosomal RNA (Eubacteria subfamily), group I intron (IC1 subfamily), tmRNA, and tRNA families
    data_dir = "/scratch/lisiz/dataset/RNAStrAlign"
    out_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/test_data_new"

    fam_dirs = os.listdir(data_dir)
    for fam_dir in fam_dirs:
        dir_path = os.path.join(data_dir, fam_dir)
        if not os.path.isdir(dir_path):
            continue

        if "SRP" in fam_dir or "RNaseP" in fam_dir or "telomerase" in fam_dir or "16S" in fam_dir:
            out_dir_fam = os.path.join(out_dir, fam_dir)
            if not os.path.exists(out_dir_fam):
                os.mkdir(out_dir_fam)

            fam_name =  fam_dir.replace("_database", "")
            if "telomerase" in fam_dir: 
                fasta_file = os.path.join(dir_path, "%s.fasta" % fam_name)
                sample(os.path.join(data_dir, fam_dir), fasta_file, fam_name, out_dir_fam, 2)
            else:
                subfam_dirs = os.listdir(dir_path)
                for subfam_dir in subfam_dirs:
                    if subfam_dir == "CVS": continue
                    if subfam_dir.startswith("."): continue

                    out_dir_sumfam = os.path.join(out_dir_fam, subfam_dir)
                    if not os.path.exists(out_dir_sumfam):
                        os.mkdir(out_dir_sumfam)
                    
                    # subfam_name = subfam_dir
                    fasta_file = os.path.join(dir_path, subfam_dir, "%s.fasta" % subfam_dir)
                    # print(fasta_file)
                    sample(os.path.join(data_dir, fam_dir, subfam_dir), fasta_file, "%s_%s" % (fam_name, subfam_dir), out_dir_sumfam, 2)


def analyze_dataset():
    sample_dir = sys.argv[1]
    ground_truth_dir = sys.argv[2]

    files = os.listdir(sample_dir)
    for filename in files:
        filepath = os.path.join(sample_dir, filename)
        with open(filepath, "r") as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                name = line.strip()[1:].split()[0]
                # print(name)

                struc_file = os.path.join(ground_truth_dir, name + ".dotbracket")
                # print(struc_file)

                # if "tmRNA" in sample_dir:
                # # of pseudoknot = 4
                with open(struc_file, "r") as f2:
                    struc = f2.readlines()[2].strip()
                
                stack = []
                page = 0
                struc_ratio = 0
                for i, item in enumerate(struc):
                    if item != ".":
                        struc_ratio += 1
                    if item == "<":
                        stack.append(i)
                    elif item == ">":
                        stack.pop()
                        if len(stack) == 0:
                            page += 1

                max_length_of_loop = 0 
                start = 0
                while struc[start] == ".":
                    start += 1
                end = len(struc) - 1
                while struc[end] == ".":
                    end -= 1
                newstruc = struc[start: end+1]
                loop_len = 0
                for i, item in enumerate(struc):
                    if item == ".":
                        loop_len += 1
                    else:
                        if loop_len > 0:
                            max_length_of_loop = max(loop_len, max_length_of_loop)
                        loop_len = 0
                        

                assert len(stack) == 0
                print(filename, name, page, struc_ratio/len(struc), max_length_of_loop)
                # print(struc)


                # max unpaired regions 

def deviation():
    ret_dir = sys.argv[1]
    files = os.listdir(ret_dir)

    for filename in files:
        filepath = os.path.join(ret_dir, filename)
        with open(filepath, "r") as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            line = line.strip()
            if not line: continue
            if line.startswith("viterbi path score"):
                aln_best_score = float(line.split()[-1])
                # print("aln_best_score", aln_best_score)
                continue
            if line.startswith("backward score"):
                seq1_mfe = lines[i+1].strip().split()[-1][1:-1]
                seq2_mfe = lines[i+2].strip().split()[-1][1:-1]
                # print("seq mfe", seq1_mfe, seq2_mfe)
                continue
            if line.startswith("inside"):
                # inside: -4572.73 2650 2370 -191.855
                items = line.split()
                print(filename, aln_best_score, seq1_mfe, seq2_mfe, int(items[2]), int(items[3]), float(items[4]), float(items[1]))
                break


def debug2():
    debug_file = sys.argv[1]
    with open(debug_file, "r") as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line: continue
            # print(line)

            # if "P2P: 2 28 42 27 41 27 43 26 42" in line:
            # if "M2=M+P: 27 43 26 42 10 43 10 42" in line:
            # if "M2 to M: 10 43 10 42" in line:
            # if "M + U 10 43 10 42" in line:
            if "M + U 10 44 10 43" in line:
                count = 0
                while count < 30:
                    print(lines[i+count].strip())
                    count += 1
                break


def multilign_data_generation():
    config_dir = "/mnt/home/sizhen/benchmark/RNAstructure/exe/configs"
    out_dir = "/mnt/home/sizhen/dataset/ltf_data/test/test2_fasta"
    files = os.listdir(config_dir)

    for filename in files:
        filepath = os.path.join(config_dir, filename)
        data = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if not line: continue
                if "SequenceNumber" in line: continue
                if line.startswith("Seq"):
                    # print(line)
                    seq_path = line.split(" = ")[-1]
                    # print(seq_path)
                    name = os.path.basename(seq_path).split("-", 1)[-1]
                    with open(seq_path, "r") as f1:
                        seq = f1.readlines()[2].strip()[:-1]
                        # print(name, seq)
                    data.append((name, seq))

        out_file = os.path.join(out_dir, filename.split(".")[0] + ".fasta")
        print(out_file)
        with open(out_file, "w") as f:
            for (name, seq) in data:
                f.write(">%s\n" % name)
                f.write("%s\n" % seq) 


def replace_rnasep_seqs():
    data_dir = "/mnt/home/sizhen/linearsankoff/data/test_data/RNaseP"
    new_dir = "/mnt/home/sizhen/linearsankoff/data/test_data/RNaseP_new"
    files = os.listdir(data_dir)

    seq_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"

    for filename in files:
        filepath = os.path.join(data_dir, filename)
        newfile = os.path.join(new_dir, filename)

        with open(newfile, "w") as wf:
            with open(filepath, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if line.startswith(">"):
                        name = line[1:].split()[0]
                        print(name)

                        seqfile = os.path.join(seq_dir, name + ".seq")
                        print(seqfile)

                        with open(seqfile, "r") as f:
                            seq = f.readlines()[-1].strip()[:-1]
                        print(seq)

                        # write
                        wf.write(line + "\n")
                        wf.write(seq + "\n")
                    else:
                        print(line.strip())
                    

# debug_data_generation()
# train_data_generation()
# debug_data()
# seq_identityALL()

# debug()
# debug2()

test_data_generation()
# analyze_dataset()

# deviation()

# multilign_data_generation()


# replace_rnasep_seqs()
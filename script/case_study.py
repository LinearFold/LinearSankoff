from html.entities import name2codepoint
import os
from random import sample
import subprocess
from Bio import SeqIO
import sys

def get_num_hairpins(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
        seq = lines[1].strip()
        struc = lines[2].strip()

        # detect number of hairpins 
        pre = ""
        num_hairpin = 0
        for s in struc:
            if s == ".":
                continue
            if s == "(":
                pre = "("
                continue
            if s == ")":
                if pre == "(": # hairpin 
                    num_hairpin += 1
                pre = ")"

    return num_hairpin

# branch insertion tRNA
def branch_insertion(data_dir):
    # data_dir = "/scratch/lisiz/dataset/RNAStrAlign/tRNA_database"
    # data_dir = "/scratch/lisiz/dataset/RNAStrAlign/SRP_database/archael"
    
    files = os.listdir(data_dir)
    dotbrackets_files = [(data_dir, file) for file in files if file.endswith(".dotbracket")]

    # data_dir = "/scratch/lisiz/dataset/RNAStrAlign/SRP_database/long_bacterial" # archael, protozoan
    # files = os.listdir(data_dir)
    # dotbrackets_files.extend([(data_dir, file) for file in files if file.endswith(".dotbracket")])

    branches = {}
    name2seq = {}
    name2struc = {}
    for data_dir, filename in dotbrackets_files:
        filepath = os.path.join(data_dir, filename)
        with open(filepath, "r") as f:
            lines = f.readlines()
            seq = lines[1].strip()
            name2seq[filename] = seq
            struc = lines[2].strip()
            name2struc[filename] = struc

            # detect number of hairpins 
            pre = ""
            num_hairpin = 0
            for s in struc:
                if s == ".":
                    continue
                if s == "(":
                    pre = "("
                    continue
                if s == ")":
                    if pre == "(": # hairpin 
                        num_hairpin += 1
                    pre = ")"

            if num_hairpin not in branches:
                branches[num_hairpin] = set()
            branches[num_hairpin].add(filename)

    for num_hairpin in branches:
        print(num_hairpin, len(branches[num_hairpin]))
        if num_hairpin == 9:
            for filename in branches[num_hairpin]:
                print(filename, name2seq[filename], name2struc[filename])


# branch insertion 
def branch_insertion2(data_dir):
    files = os.listdir(data_dir)
    dotbrackets_files = [(data_dir, file) for file in files if file.endswith(".dotbracket")]

    for data_dir, filename in dotbrackets_files:
        filepath = os.path.join(data_dir, filename)

        branches = []
        with open(filepath, "r") as f:
            lines = f.readlines()
            # seq = lines[1].strip()
            struc = lines[2].strip()

            # detect number of hairpins 
            pairs = {}
            stack = []
            pre = ""
            close_pairs = {}
            pre_pair = ()
            for i, s in enumerate(struc):
                if s == ".":
                    continue
                if s == "(":
                    stack.append(i)
                    if pre == ")":
                        left, right = pre_pair
                        print(pre_pair)
                        close_pairs[left] = right
                    pre = "("
                    continue
                if s == ")":
                    left = stack.pop()
                    if len(stack) == 0:
                        close_pairs[left] = i
                    pre_pair = (left, i)
                    pre = ")"
                    pairs[left] = i

        print(close_pairs, struc)
        for left, right in close_pairs.items():
            i = left + 1
            num_branch = 1
            while i < right:
                # print(left, right, i)
                if i in pairs:
                    i = pairs[i]
                    num_branch += 1
                else:
                    i += 1
            # if num_branch > 1:
            #     print(left, right, num_branch)


def ct2dot():
    data_dir = "/scratch/lisiz/dataset/RNAStrAlign/"
    fams = os.listdir(data_dir)
    for fam in fams:
        if fam == "ReadMe.txt":
            continue
        if "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
            files = os.listdir(os.path.join(data_dir, fam))
            for filename in files:
                if "CVS" in filename: continue
                if not filename.endswith(".ct"): continue 
                filepath = os.path.join(data_dir, fam, filename)
                outpath = os.path.join(data_dir, fam, filename.replace(".ct", ".dotbracket"))
                if os.path.exists(outpath): continue
                subprocess.run(["/nfs/stak/users/lisiz/RNAstructure/exe/ct2dot", filepath, "1", outpath])
        else:
            subfams = os.listdir(os.path.join(data_dir, fam))
            for subfam in subfams:
                files = os.listdir(os.path.join(data_dir, fam, subfam))
                for filename in files:
                    if "CVS" in filename: continue
                    if not filename.endswith(".ct"): continue 
                    filepath = os.path.join(data_dir, fam, subfam, filename)
                    outpath = os.path.join(data_dir, fam, subfam, filename.replace(".ct", ".dotbracket"))
                    if os.path.exists(outpath): continue
                    subprocess.run(["/nfs/stak/users/lisiz/RNAstructure/exe/ct2dot", filepath, "1", outpath])

# branch_insertion
def dynalignii_result_analysis():
    ret_dir = "/scratch/lisiz/lso_ret/dynalignii_ret/four_four"
    subdirs = os.listdir(ret_dir)
    for subdir in subdirs:
        files = os.listdir(os.path.join(ret_dir, subdir))
        db_files = [os.path.join(ret_dir, subdir, file) for file in files if file.endswith(".dotbracket")]
        num1 = get_num_hairpins(db_files[0])
        num2 = get_num_hairpins(db_files[1])
        # print(num1, num2)
        # print(get_num_hairpins(db_files[0]), get_num_hairpins(db_files[1]))
        if num1 == 3 and num2 == 3:
            pass
        else:
            print(get_num_hairpins(db_files[0]), get_num_hairpins(db_files[1]))

def get_clean_align(aln1, aln2):
    seq1_pos = 0
    seq2_pos = 0
    mapping = {}
    for a1, a2 in zip(aln1, aln2):
        if a1 == "-" and a1 == a2:
            continue
        
        if a1 != "-":
            seq1_pos += 1
        if a2 != "-":
            seq2_pos += 1
            mapping[seq1_pos] = seq2_pos
        else:
            mapping[seq1_pos] = 0

    return mapping

def align_inter_loop(name2aln, name2struc):
    for name1 in name2struc:
        for name2 in name2struc:
            if name1 == name2:
                continue
            mapping = get_clean_align(name2aln[name1], name2aln[name2])

            struc1 = name2struc[name1]
            struc2 = name2struc[name2]
            size2 = len(struc2)
            stack = []
            prel = 0
            prer = 0
            for i, s in enumerate(struc1, 1):
                if s == ".":
                    continue
                if s == "(":
                    stack.append(i)
                    continue
                if s == ")":
                    left = stack.pop()
                    right = i
                    if prel:
                        # print(left, prel, prer, right)
                        if left+1 < prel or prer < right-1: # not stacking
                            loop1 = len(struc1[left: prel-1].replace(".", ""))
                            loop2 = len(struc1[prer: right-1].replace(".", ""))
                            if loop1 == 0 and loop2 == 0: # internal loop
                                #left, prel, right, prer
                                if mapping[prel] and mapping[left] and mapping[right] and mapping[prer]:
                                    left2 = mapping[left]
                                    prel2 = mapping[prel]
                                    prer2 = mapping[prer]
                                    right2 = mapping[right]

                                    if struc2[left2-1] == "(" and struc2[right2-1] == ")":
                                        pass
                                    elif struc2[prel2-1] == "(" and struc2[prer2-1] == ")":
                                        pass
                                    else:  
                                        print(left, prel, prer, right, struc1[left-1: prel], struc1[prer-1: right])
                                        print(left2, prel2, prer2, right2, struc2[left2-1: prel2], struc2[prer2-1: right2],struc2[left2-1-5: prel2+5], struc2[prer2-1-5: right2+5])
                                        print("")
                                    
                                else:
                                    # print(name2aln[name1], name2aln[name2], struc1, struc2)
                                    # break
                                    pass

                    prel = left
                    prer = right


def length_diff_bw_align_inter_loop():
    # train data
    data_dir = "/scratch/lisiz/dataset/RNAStrAlign/"
    fams = os.listdir(data_dir)
    for fam in fams:
        if fam == "ReadMe.txt":
            continue

        if "tRNA" in fam or "tmRNA" in fam:
            # read alignment
            alnfile = os.path.join(data_dir, fam, "%s.fasta" % fam.replace("_database", ""))
            name2aln = {}
            for record in SeqIO.parse(alnfile, "fasta"):
                name2aln[record.id] = str(record.seq)
        
            # read structures
            name2struc = {}
            files = os.listdir(os.path.join(data_dir, fam))
            for filename in files:
                if "CVS" in filename: continue
                if not filename.endswith(".dotbracket"): 
                    continue 
                filepath = os.path.join(data_dir, fam, filename)
                with open(filepath, "r") as f:
                    lines = f.readlines()
                    name2struc[filename.replace(".dotbracket", "")] = lines[2].strip()
               
            # print(fam, len(name2aln), len(name2struc))
            align_inter_loop(name2aln, name2struc)


        elif "5S" in fam or "intron" in fam:
            if "5S" in fam:
                subfam = "Bacteria"
            else:
                subfam = "IC1"

            # read alignment
            alnfile = os.path.join(data_dir, fam, subfam, "%s.fasta" % subfam)
            name2aln = {}
            for record in SeqIO.parse(alnfile, "fasta"):
                name2aln[record.id] = str(record.seq)
             
            # read structures
            name2struc = {}
            files = os.listdir(os.path.join(data_dir, fam, subfam))
            for filename in files:
                if "CVS" in filename: continue
                if not filename.endswith(".dotbracket"): 
                    continue 
                filepath = os.path.join(data_dir, fam, subfam, filename)
                with open(filepath, "r") as f:
                    lines = f.readlines()
                    name2struc[filename.replace(".dotbracket", "")] = lines[2].strip()

            # print(fam, len(name2aln), len(name2struc))
            align_inter_loop(name2aln, name2struc)

        else:
            pass


def layout_in_alignment(ret_dir):
    fam_dirs = os.listdir(ret_dir)
    fam_data = {}
    for fam in fam_dirs:
        out_dir = os.path.join(ret_dir, fam)
        files = os.listdir(out_dir)
        for filename in files:
            filepath = os.path.join(out_dir, filename)
            # print(filepath)
            with open(filepath, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if "sum of seq length:" in line:
                        seq1_len = int(line.strip().split()[-2])
                        seq2_len = int(line.strip().split()[-1])
                        continue
                    if "average range:" in line:
                        average = float(line.strip().split()[-1])
                        # print(filename, seq1_len * average, seq1_len * seq2_len, average/seq2_len)
                        if fam not in fam_data:
                            fam_data[fam] = []
                        fam_data[fam].append(average * seq1_len / (seq1_len + seq2_len) / 100)
                        break

    for fam, data in fam_data.items():
        print(fam, sum(data)/len(data), max(data))


# branch_insertion2(sys.argv[1])

# dynalignii_result_analysis()
# ct2dot()

# length_diff_bw_align_inter_loop()

layout_in_alignment(sys.argv[1])
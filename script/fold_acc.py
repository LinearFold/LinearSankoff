from __future__ import division

import os
import sys
import subprocess
from xmlrpc.client import FastParser

from collections import defaultdict
from hopcroftkarp import HopcroftKarp

ret_dir = sys.argv[1]

# ground_truth = "/scratch2/mtdata/lisiz/dataset/RNAStrAlign"
ground_truth = "/scratch/lisiz/dataset/RNAStrAlign"


def pairs(s):
    stack1, stack2, stack3, stack4 = [], [], [], []
    ret = set()

    for i, x in enumerate(s, 1):
        if x == '(':
            stack1.append(i)
        elif x == ')':
            ret.add((stack1.pop(), i))
        elif x == '[':
            stack2.append(i)
        elif x == ']':
            ret.add((stack2.pop(), i))
        elif x == '{':
            stack3.append(i)
        elif x == '}':
            ret.add((stack3.pop(), i))
        elif x == '<':
            stack4.append(i)
        elif x == '>':
            ret.add((stack4.pop(), i))

    return ret


def eval(gold_pairs, test_pairs):
    slip_matched, non_slip_matched = 0, 0
    gold, test = len(gold_pairs), len(test_pairs)

    graph = defaultdict(set)

    for (i, j) in test_pairs:
        # non slip
        if (i, j) in gold_pairs:
            # print i,j
            non_slip_matched += 1

        for (x, y) in [(i, j), (i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]:
            if (x, y) in gold_pairs:
                graph[(i, j)].add((str(x), str(y)))

    slip_matched = len(HopcroftKarp(graph).maximum_matching()) / 2

    ###########

    return slip_matched / (test + 1e-6), slip_matched / (gold + 1e-6), non_slip_matched / (
                test + 1e-6), non_slip_matched / (gold + 1e-6)

def acc(name, gold, test):
    P_slip, R_slip, P_non_slip, R_non_slip = eval(pairs(gold), pairs(test))

    if (P_slip + R_slip) > 0: 
        f1 = 2 * P_slip * R_slip / (P_slip + R_slip)
    else: 
        f1 = 0

    return P_slip, R_slip, f1, P_non_slip, R_non_slip 

def process_file(fam, file_path, ref_dir):
    # print(file)
    with open(file_path, "r") as f:
        lines = f.readlines()[:-1]
    # if "time" not in lines[-1]: continue
    # lines = lines[:-1]
    if len(lines) == 0: return 
    if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.":
        lines = lines[:-1]
    if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": return
    # print(file)

    seq1_name = lines[0].strip().split()[0][1:]
    seq1_struc = lines[-4].strip()
    seq2_name = lines[2].strip().split()[0][1:]
    seq2_struc = lines[-3].strip()

    ## remove singleton base pairs
    new1_struc = list(seq1_struc)
    new2_struc = list(seq2_struc)
    # remove singleton base pairs: struc1 
    stack = []
    pairs = set()
    unpaired = set()
    for i, s in enumerate(seq1_struc):
        if s == "(":
            stack.append(i)
        elif s == ")":
            left = stack.pop()
            pairs.add((left, i))
        else:
            assert s == ".", seq1_struc
            unpaired.add(i)
    for (left, right) in pairs:
        if (left-1) in unpaired and (left+1) in unpaired:
            if (right-1) in unpaired and (right+1) in unpaired:
                new1_struc[left] = "."
                new1_struc[right] = "."
    # remove singleton base pairs: struc2
    stack = []
    pairs = set()
    unpaired = set()
    for i, s in enumerate(seq2_struc):
        if s == "(":
            stack.append(i)
        elif s == ")":
            left = stack.pop()
            pairs.add((left, i))
        else:
            assert s == "."
            unpaired.add(i)
    for (left, right) in pairs:
        if (left-1) in unpaired and (left+1) in unpaired:
            if (right-1) in unpaired and (right+1) in unpaired:
                new2_struc[left] = "."
                new2_struc[right] = "."
    # remove singleton base pairs: convert list to string 
    new1_struc = "".join(new1_struc)
    new2_struc = "".join(new2_struc)

    # print(seq1_name)
    # print(seq2_struc)
    # print(seq2_name)
    # print(seq2_struc)

    seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
    seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

    # print(seq1_ground_truth)
    assert os.path.exists(seq1_ground_truth), seq1_ground_truth
    assert os.path.exists(seq2_ground_truth), seq2_ground_truth

    seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
    seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

    if not os.path.exists(seq1_ground_truth_db):
        subprocess.run(["/mnt/home/sizhen/benchmark/modified/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
    if not os.path.exists(seq2_ground_truth_db):
        subprocess.run(["/mnt/home/sizhen/benchmark/modified/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

    # accuracy
    with open(seq1_ground_truth_db, "r") as f:
        seq1_ref = f.readlines()[2].strip()
    # subprocess.run(["python", "new_eval_He_Zhang.py", fam, seq1_ref, new1_struc])
    P_slip1, R_slip1, f11, P_non_slip1, R_non_slip1 = acc(fam, seq1_ref, new1_struc)
    
    with open(seq2_ground_truth_db, "r") as f:
        seq2_ref = f.readlines()[2].strip()
    P_slip2, R_slip2, f12, P_non_slip2, R_non_slip2 = acc(fam, seq2_ref, new2_struc)

    return P_slip1, R_slip1, f11, P_slip2, R_slip2, f12


def process_file2(seq_name, ct_path, ref_dir):
    struc_path = ct_path.replace(".ct", ".dotbracket")
    # print(os.path.exists(struc_path))
    if not os.path.exists(struc_path):
        subprocess.run(["./ct2dot", ct_path, "1", struc_path])
    with open(struc_path, "r") as f:
        seq_struc = f.readlines()[2].strip()
    
    ground_truth = os.path.join(ref_dir, "%s.dotbracket" % seq_name)
    assert os.path.exists(ground_truth), ground_truth
    with open(ground_truth, "r") as f:
        seq_ref = f.readlines()[2].strip()

    print(seq_ref)
    print(seq_struc)
    
    # accuracy
    # subprocess.run(["python", "new_eval_He_Zhang.py", fam, seq1_ref, new1_struc])
    P_slip1, R_slip1, f11, P_non_slip1, R_non_slip1 = acc(seq_name, seq_ref, seq_struc)
    return  P_slip1, R_slip1, f11, P_non_slip1, R_non_slip1


# folding accuracy
def ours():
    result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7_newnew/%s" % ret_dir
    print(result_dir)
    files = os.listdir(result_dir)

    ppv_acc = defaultdict(list)
    sen_acc = defaultdict(list)

    fams = os.listdir(result_dir)
    for fam in fams:
        # print(fam)
        fam_dir = os.path.join(result_dir, fam)

        # if "groupIintron" in fam: continue

        if "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database")
        elif "5S" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_rRNA_database", "Bacteria")
        elif "16S" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_rRNA_database", "Alphaproteobacteria")
        elif "RNaseP" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database", "b_bacterial")
        elif "SRP" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database", "protozoan")
        elif "groupIintron" in fam:
            ref_dir = os.path.join(ground_truth, "group_I_intron_database", "IC1")
        else:
            pass

        files = os.listdir(fam_dir)
        for filename in files:
            # print(filename)
            filepath = os.path.join(fam_dir, filename)
            P_slip1, R_slip1, f11, P_slip2, R_slip2, f12 = process_file(filename, filepath, ref_dir)
            print(filename, P_slip1, R_slip1, f11*100)
            print(filename, P_slip2, R_slip2, f12*100)
            
            # seq_idntty = int(float(filename.split("_")[-1].replace(".fasta", "")) * 10)
            # print("seq_idntty", filename, seq_idntty)

            ppv_acc[fam].extend([P_slip1, P_slip2])
            sen_acc[fam].extend([R_slip1, R_slip2])

        # else:
        #     subfams = os.listdir(fam_dir)
        #     for subfam in subfams:
        #         ref_dir = os.path.join(ground_truth, fam, subfam)

        #         subfam_dir = os.path.join(fam_dir, subfam)
        #         files = os.listdir(subfam_dir)
        #         for filename in files:
        #             filepath = os.path.join(subfam_dir, filename)
        #             P_slip1, R_slip1, f11, P_slip2, R_slip2, f12 = process_file(filename, filepath, ref_dir)

        #             # seq_idntty = int(float(filename.split("_")[-1].replace(".fasta", "")) * 10)
        #             # print("seq_idntty", filename, seq_idntty)

        #             ppv_acc[fam].extend([P_slip1, P_slip2])
        #             sen_acc[fam].extend([R_slip1, R_slip2])
    
    for fam in ppv_acc.keys():
        ppv = sum(ppv_acc[fam])/len(ppv_acc[fam])
        sen = sum(sen_acc[fam])/len(sen_acc[fam])
        print(fam, len(ppv_acc[fam]), ppv, sen, 2*ppv*sen/(ppv+sen))
    
def dynalign():
    result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7_newnew/%s" % ret_dir
    print(result_dir)
    files = os.listdir(result_dir)

    ppv_acc = defaultdict(list)
    sen_acc = defaultdict(list)
    ppv_acc_noslip = defaultdict(list)
    sen_acc_noslip = defaultdict(list)

    fams = os.listdir(result_dir)
    for fam in fams:
        if "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database")
        elif "5S" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_rRNA_database", "Bacteria")
        elif "16S" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_rRNA_database", "Alphaproteobacteria")
        elif "RNaseP" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database", "b_bacterial")
        elif "SRP" in fam:
            ref_dir = os.path.join(ground_truth, fam + "_database", "protozoan")
        elif "groupIintron" in fam:
            ref_dir = os.path.join(ground_truth, "group_I_intron_database", "IC1")
        else:
            pass
        

        fam_dir = os.path.join(result_dir, fam)

        # if "tRNA" in fam or "tmRNA" in fam:
        # ref_dir = os.path.join(ground_truth, fam)

        outdirs = os.listdir(fam_dir)
        for outdir in outdirs:
            outdirpath = os.path.join(fam_dir, outdir)

            files = os.listdir(outdirpath)
            ct_files = [filename for filename in files if filename.endswith(".ct")]
            # print(fam, ct_files)

            for ct_file in ct_files:
                ctpath = os.path.join(outdirpath, ct_file)
                P_slip, R_slip, f1, P_noslip, R_noslip = process_file2(ct_file.replace(".ct", ""), ctpath, ref_dir)

                # seq_idntty = int(float(outdir.split("_")[-1].replace(".fasta", "")) * 10)

                ppv_acc[fam].append(P_slip)
                sen_acc[fam].append(R_slip)
                ppv_acc_noslip[fam].append(P_noslip)
                sen_acc_noslip[fam].append(R_noslip)
                print(ct_file, P_slip, R_slip, P_noslip, R_noslip)

        # else:
        #     subfams = os.listdir(fam_dir)
        #     for subfam in subfams:
        #         ref_dir = os.path.join(ground_truth, fam, subfam)

        #         subfam_dir = os.path.join(fam_dir, subfam)
        #         outdirs = os.listdir(subfam_dir)
        #         for outdir in outdirs:
        #             outdirpath = os.path.join(subfam_dir, outdir)

        #             files = os.listdir(outdirpath)
        #             ct_files = [filename for filename in files if filename.endswith(".ct")]
                    
        #             for ct_file in ct_files:
        #                 ctpath = os.path.join(outdirpath, ct_file)
        #                 P_slip, R_slip, f1, _, _ = process_file2(ct_file.replace(".ct", ""), ctpath, ref_dir)

        #                 seq_idntty = int(float(outdir.split("_")[-1].replace(".fasta", "")) * 10)

        #                 ppv_acc[fam].append(P_slip)
        #                 sen_acc[fam].append(R_slip)

    
    for fam in ppv_acc.keys():
        ppv = sum(ppv_acc[fam])/len(ppv_acc[fam])
        sen = sum(sen_acc[fam])/len(sen_acc[fam])
        ppv_noslip = sum(ppv_acc_noslip[fam])/len(ppv_acc_noslip[fam])
        sen_noslip = sum(sen_acc_noslip[fam])/len(sen_acc_noslip[fam])
        print(fam, len(ppv_acc[fam]), ppv, sen, 2*ppv*sen/(ppv+sen), ppv_noslip, sen_noslip, 2*ppv_noslip*sen_noslip/(ppv_noslip+sen_noslip))

def ltf():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/%s" % ret_dir
    ltf_ret_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/ltf_ret/test_data"
    files = os.listdir(our_result_dir)

    for file in files:
        #if "group_I_intron" not in file: continue
        #if int(file.replace(".fasta", "")[-2:]) >= 20: continue
        file_path = os.path.join(our_result_dir, file)
        # print(file_path)
        if "tRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tRNA")
            # continue
        elif "5S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "5S")
            # continue
        elif "tmRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tmRNA")
            # continue
        elif "group_I_intron" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1"
            dyn_dir = os.path.join(ltf_ret_dir, "groupIintron")
            # continue
        elif "telomerase" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database"
            dyn_dir = os.path.join(ltf_ret_dir, "telomerase")
            # continue
        elif "SRP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan"
            dyn_dir = os.path.join(ltf_ret_dir, "SRP")
            # continue
        elif "RNaseP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"
            dyn_dir = os.path.join(ltf_ret_dir, "RNaseP")
            # continue

        elif "16S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "16S")
            # continue
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        seq1_len = lines[0].strip().split()[1]
        seq2_len = lines[2].strip().split()[1]
        # print(seq1_name, seq2_name)

        out_dir = os.path.join(dyn_dir, file)
        seq1_file = os.path.join(out_dir, "1_" + seq1_name + "_%s.ct" % seq1_len)
        seq2_file = os.path.join(out_dir, "2_" + seq2_name + "_%s.ct" % seq2_len)
        # print(seq1_file)
        # print(seq2_file)

        if not  os.path.exists(seq1_file): continue
        if not  os.path.exists(seq2_file): continue

        assert os.path.exists(seq1_file), seq1_file
        assert os.path.exists(seq2_file), seq2_file

        seq1_file_db = os.path.join(out_dir, file + "_1.dotbracket")
        seq2_file_db = os.path.join(out_dir, file + "_2.dotbracket")

        if not os.path.exists(seq1_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_file, "1", seq1_file_db])
        if not os.path.exists(seq2_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_file, "1", seq2_file_db])

        # ref
        seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
        seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

        assert os.path.exists(seq1_ground_truth), seq1_ground_truth
        assert os.path.exists(seq2_ground_truth), seq2_ground_truth

        seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
        seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

        if not os.path.exists(seq1_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
        if not os.path.exists(seq2_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

        # accuracy
        with open(seq1_ground_truth_db, "r") as f:
            seq1_ref = f.readlines()[2].strip()
        # print(seq1_struc)
        # print(seq1_ref)
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq1_ref, open(seq1_file_db, "r").readlines()[2].strip()])
        with open(seq2_ground_truth_db, "r") as f:
            seq2_ref = f.readlines()[2].strip()
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq2_ref, open(seq2_file_db, "r").readlines()[2].strip()])



def fold():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/%s" % ret_dir
    ltf_ret_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/fold_ret/test_data"
    files = os.listdir(our_result_dir)

    for file in files:
        #if "telomerase" not in file: continue
        #if int(file.replace(".fasta", "")[-2:]) >= 20: continue
        file_path = os.path.join(our_result_dir, file)
        # print(file_path)
        if "tRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tRNA")
            # continue
        elif "5S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "5S")
            # continue
        elif "tmRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tmRNA")
            # continue
        elif "group_I_intron" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1"
            dyn_dir = os.path.join(ltf_ret_dir, "groupIintron")
            # continue
        elif "telomerase" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database"
            dyn_dir = os.path.join(ltf_ret_dir, "telomerase")
            # continue
        elif "SRP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan"
            dyn_dir = os.path.join(ltf_ret_dir, "SRP")
            # continue
        elif "RNaseP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"
            dyn_dir = os.path.join(ltf_ret_dir, "RNaseP")
            # continue

        elif "16S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "16S")
            # continue
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        seq1_len = lines[0].strip().split()[1]
        seq2_len = lines[2].strip().split()[1]
        # print(seq1_name, seq2_name)

        out_dir = os.path.join(dyn_dir, file)
        seq1_file = os.path.join(out_dir + "_1.ct")
        seq2_file = os.path.join(out_dir + "_2.ct")
        # print(seq1_file)
        # print(seq2_file)

        if not  os.path.exists(seq1_file): continue
        if not  os.path.exists(seq2_file): continue

        assert os.path.exists(seq1_file), seq1_file
        assert os.path.exists(seq2_file), seq2_file

        seq1_file_db = os.path.join(out_dir + "_1.dotbracket")
        seq2_file_db = os.path.join(out_dir + "_2.dotbracket")

        if not os.path.exists(seq1_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_file, "1", seq1_file_db])
        if not os.path.exists(seq2_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_file, "1", seq2_file_db])

        # ref
        seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
        seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

        assert os.path.exists(seq1_ground_truth), seq1_ground_truth
        assert os.path.exists(seq2_ground_truth), seq2_ground_truth

        seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
        seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

        if not os.path.exists(seq1_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
        if not os.path.exists(seq2_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

        # accuracy
        with open(seq1_ground_truth_db, "r") as f:
            seq1_ref = f.readlines()[2].strip()
        # print(seq1_struc)
        # print(seq1_ref)
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq1_ref, open(seq1_file_db, "r").readlines()[2].strip()])
        with open(seq2_ground_truth_db, "r") as f:
            seq2_ref = f.readlines()[2].strip()
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq2_ref, open(seq2_file_db, "r").readlines()[2].strip()])



def ltf():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/%s" % ret_dir
    ltf_ret_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/ltf_ret/test_data"
    files = os.listdir(our_result_dir)

    for file in files:
        #if "group_I_intron" not in file: continue
        #if int(file.replace(".fasta", "")[-2:]) >= 20: continue
        file_path = os.path.join(our_result_dir, file)
        # print(file_path)
        if "tRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tRNA")
            # continue
        elif "5S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "5S")
            # continue
        elif "tmRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tmRNA")
            # continue
        elif "group_I_intron" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1"
            dyn_dir = os.path.join(ltf_ret_dir, "groupIintron")
            # continue
        elif "telomerase" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database"
            dyn_dir = os.path.join(ltf_ret_dir, "telomerase")
            # continue
        elif "SRP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan"
            dyn_dir = os.path.join(ltf_ret_dir, "SRP")
            # continue
        elif "RNaseP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"
            dyn_dir = os.path.join(ltf_ret_dir, "RNaseP")
            # continue

        elif "16S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "16S")
            # continue
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        seq1_len = lines[0].strip().split()[1]
        seq2_len = lines[2].strip().split()[1]
        # print(seq1_name, seq2_name)

        out_dir = os.path.join(dyn_dir, file)
        seq1_file = os.path.join(out_dir, "1_" + seq1_name + "_%s.ct" % seq1_len)
        seq2_file = os.path.join(out_dir, "2_" + seq2_name + "_%s.ct" % seq2_len)
        # print(seq1_file)
        # print(seq2_file)

        if not  os.path.exists(seq1_file): continue
        if not  os.path.exists(seq2_file): continue

        assert os.path.exists(seq1_file), seq1_file
        assert os.path.exists(seq2_file), seq2_file

        seq1_file_db = os.path.join(out_dir, file + "_1.dotbracket")
        seq2_file_db = os.path.join(out_dir, file + "_2.dotbracket")

        if not os.path.exists(seq1_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_file, "1", seq1_file_db])
        if not os.path.exists(seq2_file_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_file, "1", seq2_file_db])

        # ref
        seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
        seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

        assert os.path.exists(seq1_ground_truth), seq1_ground_truth
        assert os.path.exists(seq2_ground_truth), seq2_ground_truth

        seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
        seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

        if not os.path.exists(seq1_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
        if not os.path.exists(seq2_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

        # accuracy
        with open(seq1_ground_truth_db, "r") as f:
            seq1_ref = f.readlines()[2].strip()
        # print(seq1_struc)
        # print(seq1_ref)
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq1_ref, open(seq1_file_db, "r").readlines()[2].strip()])
        with open(seq2_ground_truth_db, "r") as f:
            seq2_ref = f.readlines()[2].strip()
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq2_ref, open(seq2_file_db, "r").readlines()[2].strip()])



def lf():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/%s" % ret_dir
    ltf_ret_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/lf_ret_infb/train_data"
    files = os.listdir(our_result_dir)

    for file in files:
        if "tmRNA" not in file: continue
        #if int(file.replace(".fasta", "")[-2:]) >= 20: continue
        file_path = os.path.join(our_result_dir, file)
        # print(file_path)
        if "tRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tRNA")
            # continue
        elif "5S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "5S")
            # continue
        elif "tmRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database"
            dyn_dir = os.path.join(ltf_ret_dir, "tmRNA")
            # continue
        elif "group_I_intron" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1"
            dyn_dir = os.path.join(ltf_ret_dir, "groupIintron")
            # continue
        elif "telomerase" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database"
            dyn_dir = os.path.join(ltf_ret_dir, "telomerase")
            # continue
        elif "SRP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan"
            dyn_dir = os.path.join(ltf_ret_dir, "SRP")
            # continue
        elif "RNaseP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"
            dyn_dir = os.path.join(ltf_ret_dir, "RNaseP")
            # continue

        elif "16S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria"
            dyn_dir = os.path.join(ltf_ret_dir, "16S")
            # continue
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        
        outfile = os.path.join(dyn_dir, file+".out")
        with open(outfile, "r") as f:
            lines = f.readlines()
        lines = [line.strip().split()[0] for line in lines if "(" in line]
        # print(lines)

        # ref
        seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
        seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

        assert os.path.exists(seq1_ground_truth), seq1_ground_truth
        assert os.path.exists(seq2_ground_truth), seq2_ground_truth

        seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
        seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

        if not os.path.exists(seq1_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
        if not os.path.exists(seq2_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

        # accuracy
        with open(seq1_ground_truth_db, "r") as f:
            seq1_ref = f.readlines()[2].strip()
        # print(seq1_struc)
        # print(seq1_ref)
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq1_ref, lines[0]])
        with open(seq2_ground_truth_db, "r") as f:
            seq2_ref = f.readlines()[2].strip()
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq2_ref, lines[1]])


# folding accuracy
def locarna():
    result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7_newnew/%s" % ret_dir
    print(result_dir)
    ltf_ret_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/locarna_p_ret/train_set"
    files = os.listdir(result_dir)

    for file in files:

        # if "16S" not in file: continue
        file_path = os.path.join(result_dir, file)
        # print(file, file_path)
        if "tRNA" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database"
            # continue
        elif "5S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria"
            # continue
        elif "tmRNA" in file:
            # print(file)
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database"
            # continue

        elif "group_I_intron" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1"
            # continue
        elif "telomerase" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database"
            # continue
        elif "SRP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan"
            # continue
        elif "RNaseP" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial"
            # continue

        elif "16S" in file:
            ref_dir = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria"
            # continue

        else:
            pass
        # print(file)
        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if len(lines) == 0: continue
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.":
            lines = lines[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        # print(file)

        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]

        # pred file
        pred_file = os.path.join(ltf_ret_dir, file)
        # get consensus structure
        align_seqs = {}
        consensus_struc = ""
        with open(os.path.join(ret_dir, ".out"), "r") as f:
            lines = f.readlines()
        flag = 0
        flag1 = 0
        for line in lines:
            line = line.strip()
            if not line: continue
            if line.startswith("Perform progressive alignment"):
                flag = 1
                continue
            if not flag: continue
            if line.startswith("reliability"): break
            if line.startswith("alifold"):
                flag1 = 1
                _, seq = line.split()
                consensus_struc += seq
                continue
            if flag1:
                seq = line
                if "(" in line:
                    seq = line.split()[0]
                consensus_struc += seq
                continue
            # print(line.split())
            name, seq = line.split()
            if name in align_seqs:
                align_seqs[name] += seq
            else:
                align_seqs[name] = seq

        ## remove singleton base pairs
        new1_struc = list(consensus_struc)
        new2_struc = list(consensus_struc)
        # remove singleton base pairs: struc1 
        stack = []
        pairs = set()
        unpaired = set()
        for i, s in enumerate(consensus_struc):
            if s == "(":
                stack.append(i)
            elif s == ")":
                left = stack.pop()
                pairs.add((left, i))
            else:
                assert s == "."
                unpaired.add(i)
        for (left, right) in pairs:
            if (left-1) in unpaired and (left+1) in unpaired:
                if (right-1) in unpaired and (right+1) in unpaired:
                    new1_struc[left] = "."
                    new1_struc[right] = "."
        # remove singleton base pairs: struc2
        stack = []
        pairs = set()
        unpaired = set()
        for i, s in enumerate(consensus_struc):
            if s == "(":
                stack.append(i)
            elif s == ")":
                left = stack.pop()
                pairs.add((left, i))
            else:
                assert s == "."
                unpaired.add(i)
        for (left, right) in pairs:
            if (left-1) in unpaired and (left+1) in unpaired:
                if (right-1) in unpaired and (right+1) in unpaired:
                    new2_struc[left] = "."
                    new2_struc[right] = "."
        # remove singleton base pairs: convert list to string 
        new1_struc = "".join(new1_struc)
        new2_struc = "".join(new2_struc)

        # print(seq1_name)
        # print(seq2_struc)
        # print(seq2_name)
        # print(seq2_struc)

        seq1_ground_truth = os.path.join(ref_dir, "%s.ct" % seq1_name)
        seq2_ground_truth = os.path.join(ref_dir, "%s.ct" % seq2_name)

        # print(seq1_ground_truth)
        assert os.path.exists(seq1_ground_truth), seq1_ground_truth
        assert os.path.exists(seq2_ground_truth), seq2_ground_truth

        seq1_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq1_name)
        seq2_ground_truth_db = os.path.join(ref_dir, "%s.dotbracket" % seq2_name)

        if not os.path.exists(seq1_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/modified/RNAstructure/exe/ct2dot", seq1_ground_truth, "1", seq1_ground_truth_db])
        if not os.path.exists(seq2_ground_truth_db):
            subprocess.run(["/mnt/home/sizhen/benchmark/modified/RNAstructure/exe/ct2dot", seq2_ground_truth, "1", seq2_ground_truth_db])

        # accuracy
        with open(seq1_ground_truth_db, "r") as f:
            seq1_ref = f.readlines()[2].strip()
        # print(seq1_struc)
        # print(seq1_ref)
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq1_ref, new1_struc])
        with open(seq2_ground_truth_db, "r") as f:
            seq2_ref = f.readlines()[2].strip()
        subprocess.run(["python", "new_eval_He_Zhang.py", file, seq2_ref, new2_struc])
        
# ours()
dynalign()

# ltf()
# fold()
# lf()

# locarna()

import os
import sys
import subprocess
from Bio import SeqIO
from collections import defaultdict

ret_dir = sys.argv[1]
# ground_truth = "/scratch2/mtdata/lisiz/dataset/RNAStrAlign"
ground_truth = "/scratch/lisiz/dataset/RNAStrAlign"

def aln_acc(ref1, ref2, pred1, pred2):
    assert len(ref1) == len(ref2)
    assert len(pred1) == len(pred2)

    # print(ref1)
    # print(ref2)
    # print(pred1)
    # print(pred2)

    ref_pairs = set()
    seq1_pos = 0
    seq2_pos = 0
    ref_aln_len = 0
    seq_idntty = 0
    
    # newref1 = ""
    # newref2 = ""
    for nuc1, nuc2 in zip(ref1, ref2):
        if nuc1 == "-" and nuc2 == "-":
            continue
        ref_aln_len += 1

        # if nuc1 == "N" and nuc2 == "N": continue
        
        # newref1 += nuc1
        # newref2 += nuc2

        if nuc1 == "-":
            ref_pairs.add((-1, seq2_pos))
            seq2_pos += 1
        elif nuc2 == "-":
            ref_pairs.add((seq1_pos, -1))
            seq1_pos += 1
        else:
            ref_pairs.add((seq1_pos, seq2_pos))
            seq1_pos += 1
            seq2_pos += 1

        if nuc1 == nuc2:
            seq_idntty += 1

    # print(newref1)
    # print(newref2)

    seq1_pos = 0
    seq2_pos = 0
    pred_aln_len = 0
    num_same = 0
    for nuc1, nuc2 in zip(pred1, pred2):
        if nuc1 == "-" and nuc2 == "-":
            assert False
        pred_aln_len += 1

        # if nuc1 == "N" and nuc2 == "N": continue

        if nuc1 == "-":
            if (-1, seq2_pos) in ref_pairs: 
                num_same += 1
            seq2_pos += 1
        elif nuc2 == "-":
            if (seq1_pos, -1) in ref_pairs:
                num_same += 1
            seq1_pos += 1
        else:
            if (seq1_pos, seq2_pos) in ref_pairs:
                num_same += 1
            # elif (seq1_pos+1, seq2_pos) in ref_pairs:
            #     num_same += 1
            # elif (seq1_pos, seq2_pos+1) in ref_pairs:
            #     num_same += 1
            # elif (seq1_pos-1, seq2_pos) in ref_pairs:
            #     num_same += 1
            # elif (seq1_pos, seq2_pos-1) in ref_pairs:
            #     num_same += 1
            seq1_pos += 1
            seq2_pos += 1

    seq_idntty /= ref_aln_len
    p = num_same/pred_aln_len
    r = num_same/ref_aln_len

    # print(seq_idntty, p, r, 2 * p * r / (p + r))
    if p==0 and r==0:
        return seq_idntty, p, r, 0
    else:
        return seq_idntty, p, r, 2 * p * r / (p + r)
    
def process_file(aln_file, ref_aln):
    with open(aln_file, "r") as f:
        lines = f.readlines()
        seq1_aln = lines[-3].strip()
        seq2_aln = lines[-2].strip()

        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
    
    # print(seq1_aln, seq2_aln)

    # ref
    seq1_ref_aln = ""
    seq2_ref_aln = ""
    for item in SeqIO.parse(ref_aln, "fasta"):
        name = item.id
        if name == seq1_name:
            seq1_ref_aln = str(item.seq)
        if name == seq2_name:
            seq2_ref_aln = str(item.seq)
        if seq1_ref_aln and seq2_ref_aln:
            break
    
    assert seq1_ref_aln
    assert seq2_ref_aln

    seq_idntty, p, r, f1 = aln_acc(seq1_ref_aln, seq2_ref_aln, seq1_aln, seq2_aln)

    return seq_idntty, p, r, f1


def process_file2(names, aln_file, ref_aln):
    with open(aln_file, "r") as f:
        lines = f.readlines()
        seq1_aln = lines[1].strip()
        seq2_aln = lines[2].strip()

    seq1_name, seq2_name = names
    # ref
    seq1_ref_aln = ""
    seq2_ref_aln = ""
    for item in SeqIO.parse(ref_aln, "fasta"):
        name = item.id
        if name == seq1_name:
            seq1_ref_aln = str(item.seq)
        if name == seq2_name:
            seq2_ref_aln = str(item.seq)
        if seq1_ref_aln and seq2_ref_aln:
            break
    
    assert seq1_ref_aln
    assert seq2_ref_aln

    if len(seq1_aln.replace("-", "")) ==  len(seq1_ref_aln.replace("-", "")):
        seq_idntty, p, r, f1 = aln_acc(seq1_ref_aln, seq2_ref_aln, seq1_aln, seq2_aln)
    else:
        seq_idntty, p, r, f1 = aln_acc(seq1_ref_aln, seq2_ref_aln, seq2_aln, seq1_aln)

    return seq_idntty, p, r, f1

# folding accuracy
def ours():
    result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7_newnew/%s" % ret_dir
    print(result_dir)
    files = os.listdir(result_dir)

    ppv_acc = defaultdict(list)
    sen_acc = defaultdict(list)

    fams = os.listdir(result_dir)
    for fam in fams:
        fam_dir = os.path.join(result_dir, fam)

        # if "16S" in fam: continue
        
        if 
        elif "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
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

        refpathname = os.path.basename(ref_dir).replace("_database", "") + ".fasta"
        refpath = os.path.join(ref_dir, refpathname)
        # print(refpath)

        # if "tRNA" in fam or "tmRNA" in fam:
        files = os.listdir(fam_dir)
        for filename in files:
            filepath = os.path.join(fam_dir, filename)
            # refpath = os.path.join(ref_dir, "%s.fasta" % (fam.replace("_database", "")))

            # print(filepath, refpath)
            seq_idntty, p, r, f1 = process_file(filepath, refpath)
            
            # seq_idntty = int(float(filename.split("_")[-1].replace(".fasta", "")) * 10)
            # print("seq_idntty", filename, seq_idntty)

            ppv_acc[fam].append(p)
            sen_acc[fam].append(r)

        # else:
        #     subfams = os.listdir(fam_dir)
        #     for subfam in subfams:
        #         ref_dir = os.path.join(ground_truth, fam, subfam)

        #         subfam_dir = os.path.join(fam_dir, subfam)
        #         files = os.listdir(subfam_dir)
        #         for filename in files:
        #             filepath = os.path.join(subfam_dir, filename)
        #             refpath = os.path.join(ref_dir, "%s.fasta" % subfam)

        #             seq_idntty, p, r, f1 = process_file(filepath, refpath)

        #             seq_idntty = int(float(filename.split("_")[-1].replace(".fasta", "")) * 10)
        #             # print("seq_idntty", filename, seq_idntty)

        #             ppv_acc[fam].append(p)
        #             sen_acc[fam].append(r)
    
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

    fams = os.listdir(result_dir)
    for fam in fams:
        if "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
            refpath = os.path.join(ground_truth, fam + "_database", "%s.fasta" % fam)
        elif "5S" in fam:
            refpath = os.path.join(ground_truth, fam + "_rRNA_database", "Bacteria", "Bacteria.fasta")
        elif "16S" in fam:
            refpath = os.path.join(ground_truth, fam + "_rRNA_database", "Alphaproteobacteria", "Alphaproteobacteria.fasta")
        elif "RNaseP" in fam:
            refpath = os.path.join(ground_truth, fam + "_database", "b_bacterial", "b_bacterial.fasta")
        elif "SRP" in fam:
            refpath = os.path.join(ground_truth, fam + "_database", "protozoan", "protozoan.fasta")
        elif "groupIintron" in fam:
            refpath = os.path.join(ground_truth, "group_I_intron_database", "IC1", "IC1.fasta")
        else:
            pass

        fam_dir = os.path.join(result_dir, fam)

        # if "tRNA" in fam or "tmRNA" in fam:
        # ref_dir = os.path.join(ground_truth, fam)

        outdirs = os.listdir(fam_dir)
        for outdir in outdirs:
            outdirpath = os.path.join(fam_dir, outdir)

            files = os.listdir(outdirpath)
            names = [filename.replace(".ct", "") for filename in files if filename.endswith(".ct")]
            # print(fam, ct_files)

            alnpath = os.path.join(outdirpath, "aln.out")
            # refpath = os.path.join(ref_dir, "%s.fasta" % (fam.replace("_database", "")))
                
            seq_idntty, p, r, f1 = process_file2(names, alnpath, refpath)

            # seq_idntty = int(float(outdir.split("_")[-1].replace(".fasta", "")) * 10)

            ppv_acc[fam].append(p)
            sen_acc[fam].append(r)

        # else:
        #     subfams = os.listdir(fam_dir)
        #     for subfam in subfams:
        #         ref_dir = os.path.join(ground_truth, fam, subfam)

        #         subfam_dir = os.path.join(fam_dir, subfam)
        #         outdirs = os.listdir(subfam_dir)
        #         for outdir in outdirs:
        #             outdirpath = os.path.join(subfam_dir, outdir)

        #             files = os.listdir(outdirpath)
        #             names = [filename.replace(".ct", "") for filename in files if filename.endswith(".ct")]
                    
        #             alnpath = os.path.join(outdirpath, "aln.out")
        #             refpath = os.path.join(ref_dir, "%s.fasta" % subfam)

        #             seq_idntty, p, r, f1 = process_file2(names, alnpath, refpath)

        #             seq_idntty = int(float(outdir.split("_")[-1].replace(".fasta", "")) * 10)

        #             ppv_acc[fam].append(p)
        #             sen_acc[fam].append(r)
    
    for fam in ppv_acc.keys():
        ppv = sum(ppv_acc[fam])/len(ppv_acc[fam])
        sen = sum(sen_acc[fam])/len(sen_acc[fam])
        print(fam, len(ppv_acc[fam]), ppv, sen, 2*ppv*sen/(ppv+sen))

def mafft():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir # "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/%s" % ret_dir
    files = os.listdir(our_result_dir)
    ltf_dir = "/mnt/storage/idl-0/bio/sizhen/mafft_ret/"

    for file in files:
        file_path = os.path.join(our_result_dir, file)
        # if "16S_rRNA_Alphaproteobacteria_group1.fasta" in file: continue
        if "tRNA" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database/tRNA.fasta"
            dyn_dir = os.path.join(ltf_dir, "tRNA")
        elif "5S" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria/Bacteria.fasta"
            dyn_dir = os.path.join(ltf_dir, "5S")
        elif "tmRNA" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database/tmRNA.fasta"
            dyn_dir = os.path.join(ltf_dir, "tmRNA")
        elif "group_I_intron" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1/IC1.fasta"
            dyn_dir = os.path.join(ltf_dir, "groupIintron")

        elif "SRP" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan/protozoan.fasta"
            dyn_dir = os.path.join(ltf_dir, "SRP")
        elif "telomerase" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database/telomerase.fasta"
            dyn_dir = os.path.join(ltf_dir, "telomerase")
        elif "RNaseP" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial/b_bacterial.fasta"
            dyn_dir = os.path.join(ltf_dir, "RNaseP")
        elif "16S" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria/Alphaproteobacteria.fasta"
            dyn_dir = os.path.join(ltf_dir, "16S")
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        # print(file)

        # out_dir = os.path.join(dyn_dir, file)
        aln_file = os.path.join(dyn_dir, file + ".ali")

        i_seq = 0
        for item in SeqIO.parse(aln_file, "fasta"):
            name = item.id
            seq = str(item.seq)
            if i_seq == 0:
                seq1_aln = seq.upper()
                i_seq += 1
            else:
                seq2_aln = seq.upper()

        # ref
        seq1_ref_aln = ""
        seq2_ref_aln = ""
        for item in SeqIO.parse(ref_aln, "fasta"):
            name = item.id
            if name == seq1_name:
                seq1_ref_aln = str(item.seq)
            if name == seq2_name:
                seq2_ref_aln = str(item.seq)
            if seq1_ref_aln and seq2_ref_aln:
                break
        
        assert seq1_ref_aln
        assert seq2_ref_aln

        # acc
        seq_idntty, p, r, f1 = aln_acc(seq1_ref_aln, seq2_ref_aln, seq1_aln, seq2_aln)
        print(file, seq_idntty, p, r, f1)
    
def ltf():
    ret_dir = sys.argv[1]
    our_result_dir = ret_dir #  "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/array_version7/test_data_w30_b500_newnewret/"
    files = os.listdir(our_result_dir)
    ltf_dir = "/mnt/storage/idl-0/bio/sizhen/linearsankoff_rets/ltf_ret/test_data"

    for file in files:
        file_path = os.path.join(our_result_dir, file)
        if "tRNA" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/tRNA_database/tRNA.fasta"
            dyn_dir = os.path.join(ltf_dir, "tRNA")
        elif "5S" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/5S_rRNA_database/Bacteria/Bacteria.fasta"
            dyn_dir = os.path.join(ltf_dir, "5S")
        elif "tmRNA" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/tmRNA_database/tmRNA.fasta"
            dyn_dir = os.path.join(ltf_dir, "tmRNA")
        elif "group_I_intron" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/group_I_intron_database/IC1/IC1.fasta"
            dyn_dir = os.path.join(ltf_dir, "groupIintron")

        elif "SRP" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/SRP_database/protozoan/protozoan.fasta"
            dyn_dir = os.path.join(ltf_dir, "SRP")
        elif "telomerase" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/telomerase_database/telomerase.fasta"
            dyn_dir = os.path.join(ltf_dir, "telomerase")
        elif "RNaseP" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/RNaseP_database/b_bacterial/b_bacterial.fasta"
            dyn_dir = os.path.join(ltf_dir, "RNaseP")
        elif "16S" in file:
            ref_aln = "/mnt/home/sizhen/dataset/RNAStrAlign/16S_rRNA_database/Alphaproteobacteria/Alphaproteobacteria.fasta"
            dyn_dir = os.path.join(ltf_dir, "16S")
        else:
            pass

        with open(file_path, "r") as f:
            lines = f.readlines()[:-1]
        if lines[-4][0] not in ")(." or lines[-3][0] not in ")(.": continue
        seq1_name = lines[0].strip().split()[0][1:]
        seq2_name = lines[2].strip().split()[0][1:]
        # print(file)
        
        out_dir = os.path.join(dyn_dir, file)
        # aln_file = os.path.join(out_dir, file + ".aln")
        aln_file = os.path.join(out_dir, "output.aln")

        # print(aln_file)
        # with open(aln_file, "r") as f:
        num = 0
        for item in SeqIO.parse(aln_file, "fasta"):
            # lines = f.readlines()
            # seq1_aln = lines[1].strip()
            # seq2_aln = lines[2].strip()
            if num == 0: 
                seq1_aln = str(item.seq)
                num += 1
            else:
                seq2_aln = str(item.seq)

        # ref
        seq1_ref_aln = ""
        seq2_ref_aln = ""
        for item in SeqIO.parse(ref_aln, "fasta"):
            name = item.id
            if name == seq1_name:
                seq1_ref_aln = str(item.seq)
            if name == seq2_name:
                seq2_ref_aln = str(item.seq)
            if seq1_ref_aln and seq2_ref_aln:
                break
        
        assert seq1_ref_aln
        assert seq2_ref_aln

        # acc
        seq_idntty, p, r, f1 = aln_acc(seq1_ref_aln, seq2_ref_aln, seq1_aln, seq2_aln)
        print(file, seq_idntty, p, r, f1)



ours()
# dynalign()
# mafft()
# ltf()

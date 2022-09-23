import sys
import os
import subprocess

fam = sys.argv[1]
version = int(sys.argv[2])

data_root = "/scratch/lisiz/dataset/RNAStrAlign"

data_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/train_data"

if version == 1:
    template_config = "/nfs/stak/users/lisiz/RNAstructure/exe/dynalign.config"
    out_dir = "/scratch/lisiz/lso_ret/dynalign_ret"
else:
    template_config = "/nfs/stak/users/lisiz/RNAstructure/exe/dynalignii.config"
    out_dir = "/scratch/lisiz/lso_ret/dynalignii_ret"

# if fam == "5S":
#     fam = "5S_rRNA_database"
# elif fam == "groupIintron":
#     fam = "group_I_intron_database"
# else:
#     fam = fam + "_database"


# if "tRNA" in fam or "tmRNA" in fam:
fam_dir = os.path.join(data_dir, fam)
fam_dir_out = os.path.join(out_dir, fam)
if not os.path.exists(fam_dir_out):
    os.makedirs(fam_dir_out)

if "tRNA" in fam or "tmRNA" in fam or "telomerase" in fam:
    fam_data_root = os.path.join(data_root, fam + "_database")
elif "5S" in fam:
    fam_data_root = os.path.join(data_root, fam + "_rRNA_database", "Bacteria")
elif "16S" in fam:
    fam_data_root = os.path.join(data_root, fam + "_rRNA_database", "Alphaproteobacteria")
elif "RNaseP" in fam:
    fam_data_root = os.path.join(data_root, fam + "_database", "b_bacterial")
elif "SRP" in fam:
    fam_data_root = os.path.join(data_root, fam + "_database", "protozoan")
elif "groupIintron" in fam:
    fam_data_root = os.path.join(data_root, "group_I_intron_database", "IC1")
else:
    pass
        
# fam_data_root = os.path.join(data_root, fam)

files = os.listdir(fam_dir)
for filename in files:
    filepath = os.path.join(fam_dir, filename)
    
    finalout = os.path.join(fam_dir_out, filename)
    if not os.path.exists(finalout):
        os.makedirs(finalout)

    with open(filepath, "r") as f:
        lines = f.readlines()
        name1 = lines[0].strip().split()[0][1:]
        name2 = lines[2].strip().split()[0][1:]
    # print(name1, name2)

    seq1_file = os.path.join(fam_data_root, name1+".seq")
    seq2_file = os.path.join(fam_data_root, name2+".seq")

    ct1_file = os.path.join(finalout, name1+".ct")
    ct2_file = os.path.join(finalout, name2+".ct")
    aout = os.path.join(finalout, "aln.out")

    # build config 
    if version == 1:
        newconfig = "configs/%s.dynalign.config" % filename
    else:
        newconfig = "configs/%s.dynalignii.config" % filename
    # if not os.path.exists(newconfig):
    
    with open(newconfig, 'w') as wf:
        with open(template_config, "r") as rf:
            lines = rf.readlines()
            for line in lines:
                if line.startswith("#"):
                    wf.write(line)
                    continue
                if "inseq1" in line:
                    wf.write("inseq1 = %s\n" % seq1_file)
                elif "inseq2" in line:
                    wf.write("inseq2 = %s\n" % seq2_file)
                elif "outct" in line and "outct2" not in line:
                    wf.write("outct = %s\n" % ct1_file)
                elif "outct2" in line:
                    wf.write("outct2 = %s\n" % ct2_file)
                elif "aout" in line:
                    wf.write("aout = %s\n" % aout)
                else:
                    wf.write(line)

    if version == 1:
        print("nohup ./dynalign %s &" % newconfig)
    else:
        print("nohup ./dynalign_ii %s &" % newconfig)

# else:
#     fam_dir = os.path.join(data_dir, fam)
#     fam_dir_out = os.path.join(out_dir, fam)
#     if not os.path.exists(fam_dir_out):
#         os.makedirs(fam_dir_out)

#     fam_data_root = os.path.join(data_root, fam)

#     sub_fams = os.listdir(fam_dir)
#     for sub_fam in sub_fams:
#         sub_fam_dir = os.path.join(fam_dir, sub_fam)
#         sum_fam_dir_out = os.path.join(fam_dir_out, sub_fam)
#         if not os.path.exists(sum_fam_dir_out):
#             os.makedirs(sum_fam_dir_out)
        
#         subfam_data_root = os.path.join(fam_data_root, sub_fam)

#         files = os.listdir(sub_fam_dir)
#         for filename in files:
#             filepath = os.path.join(sub_fam_dir, filename)
            
#             finalout = os.path.join(sum_fam_dir_out, filename)
#             if not os.path.exists(finalout):
#                 os.makedirs(finalout)

#             with open(filepath, "r") as f:
#                 lines = f.readlines()
#                 name1 = lines[0].strip().split()[0][1:]
#                 name2 = lines[2].strip().split()[0][1:]
#             # print(name1, name2)

#             seq1_file = os.path.join(subfam_data_root, name1+".seq")
#             seq2_file = os.path.join(subfam_data_root, name2+".seq")

#             ct1_file = os.path.join(finalout, name1+".ct")
#             ct2_file = os.path.join(finalout, name2+".ct")
#             aout = os.path.join(finalout, "aln.out")

#             # build config 
#             newconfig = "configs/%s.config" % filename
#             # if not os.path.exists(newconfig):
#             with open("configs/%s.config" % filename, 'w') as wf:
#                 with open("dynalign.config", "r") as rf:
#                     lines = rf.readlines()
#                     for line in lines:
#                         if line.startswith("#"):
#                             wf.write(line)
#                             continue
#                         if "inseq1" in line:
#                             wf.write("inseq1 = %s\n" % seq1_file)
#                         elif "inseq2" in line:
#                             wf.write("inseq2 = %s\n" % seq2_file)
#                         elif "outct" in line and "outct2" not in line:
#                             wf.write("outct = %s\n" % ct1_file)
#                         elif "outct2" in line:
#                             wf.write("outct2 = %s\n" % ct2_file)
#                         elif "aout" in line:
#                             wf.write("aout = %s\n" % aout)
#                         else:
#                             wf.write(line)

#             if version == 1:
#                 print("nohup ./dynalign %s &" % newconfig)
#             else:
#                 print("nohup ./dynalign_ii %s &" % newconfig)

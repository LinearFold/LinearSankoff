import sys
import os
import subprocess

fam = sys.argv[1]
w = sys.argv[2] # 40
b = sys.argv[3] # 100
LFb = sys.argv[4] # -1
energyDiff = sys.argv[5] # 0.3

if "four" in fam or "five" in fam:
    data_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/SRP_data"
    out_dir = "/scratch2/mtdata/lisiz/lso_ret/SRP_case_study/w%s_b%s_LFb%s_energyDiff%s_ret" % (w, b, LFb, energyDiff)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

elif "time_mem_data" == fam:
    data_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/"
    out_dir = "/scratch2/mtdata/lisiz/lso_ret/time_mem_data/w%s_b%s_LFb%s_energyDiff%s_ret" % (w, b, LFb, energyDiff)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

else:
    data_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/train_data"
    out_dir = "/scratch2/mtdata/lisiz/lso_ret/branch_insertion/train_w%s_b%s_LFb%s_energyDiff%s_ret" % (w, b, LFb, energyDiff)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if fam in ["16S", "SRP", "RNaseP", "telomerase"]:
        data_dir = "/nfs/stak/users/lisiz/LinearSankoff/dataset/test_data"
        out_dir = "/scratch2/mtdata/lisiz/lso_ret/branch_insertion/test_w%s_b%s_LFb%s_energyDiff%s_ret" % (w, b, LFb, energyDiff)

fam_dir = os.path.join(data_dir, fam)
fam_dir_out = os.path.join(out_dir, fam)
if not os.path.exists(fam_dir_out):
    os.makedirs(fam_dir_out)

files = os.listdir(fam_dir)
for filename in files:
    filepath = os.path.join(fam_dir, filename)
    outfile = os.path.join(fam_dir_out, filename)

    print("nohup cat %s | /nfs/stak/users/lisiz/LinearSankoff/LinearSankoff -w %s -b %s --LFb %s --energyDiff %s --astar > %s &" % (filepath, w, b, LFb, energyDiff, outfile))



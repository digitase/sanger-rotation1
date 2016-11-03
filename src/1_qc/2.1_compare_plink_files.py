
import sys

f1_file = sys.argv[1]
f2_file = sys.argv[2]

if sys.argv[3] == "fam":
    comm_col = 0
elif sys.argv[3] == "bim" :
    comm_col = 1
else:
    raise ValueError("Invalid filetype:", sys.argv[3])

with open(f1_file + "." + sys.argv[3]) as f1: 
    with open(f2_file + "." + sys.argv[3]) as f2: 
        f1_set = set((l.strip().split()[comm_col] for l in f1.readlines()))
        f2_set = set((l.strip().split()[comm_col] for l in f2.readlines()))
        print('Compare {} files. Intersect: {}, Only f1: {}, Only f2: {}'.format(sys.argv[3], len(f1_set.intersection(f2_set)), len(f1_set.difference(f2_set)), len(f2_set.difference(f1_set))))


import Needleman
import UPGMA
# read score
f = open(r"socore.txt")
raw_matrix = [line.split() for line in f]
f.close()
raw_dicts = [dict() for x in range(len(raw_matrix[0]))]
for i in range(len(raw_matrix[0])):
    raw_dicts[i] = dict(zip(raw_matrix[0], map(int, raw_matrix[i + 1])))
s_matrix = dict()
for i in range(len(raw_matrix[0])):
    s_matrix[raw_matrix[0][i]] = raw_dicts[i]
print(s_matrix)

# read sequences to sequences
sequences = []
fseq = open(r"sequence.txt")
while 1:
    line = fseq.readline()
    if not line:
        break
    sequences.append(line.strip())
fseq.close()
num = 0
for seq in sequences:
    num+=1
    print(seq,num)

# create similarity matrix
score_matrix= [ [ 0 for i in range(len(sequences)) ] for j in range(len(sequences)) ]
dist_matrix= [ [ 0 for i in range(len(sequences)) ] for j in range(len(sequences)) ]
dist_map = {}
for d1 in range(len(sequences)):
    for d2 in range(len(sequences)):
        identity,score,align1,align2 = Needleman.alignment(sequences[d1], sequences[d2], s_matrix)
        Needleman.scorerand(align1, align2, s_matrix)
        dist = Needleman.distance(sequences[d1], sequences[d2], s_matrix)
        score_matrix[d1][d2] = score
        dist_matrix[d1][d2] = dist
        dist_map[(d1,d2)] = dist


# 输出similarity 到文件
file1 = "similarity.txt"
with open(file1,"w") as f:
    for row in score_matrix:
        f.write(','.join(str(e) for e in row))
        f.write('\n')
# 输出距离矩阵到distance 文件
file2 = "distance.txt"
with open(file2,"w") as f:
    for row in dist_matrix:
        f.write(','.join(str(e) for e in row))
        f.write('\n')

# # 打印距离map
# print(dist_map)
#
# 打印upgma树
tree = UPGMA.upgma(dist_map,len(sequences),len(sequences))
MSAseqs = Needleman.MSA(tree,sequences,score_matrix,s_matrix)

# again create msa_dist_map according the msaseqs
msa_dist_map = {}
for d1 in range(len(MSAseqs)):
    for d2 in range(len(MSAseqs)):
        dist = Needleman.distance(MSAseqs[d1], MSAseqs[d2], s_matrix)
        msa_dist_map[(d1,d2)] = dist

finaltree = UPGMA.upgma(msa_dist_map,len(MSAseqs),len(MSAseqs))
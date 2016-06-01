import math
def zeros(shape):
    retval = []
    for x in range(shape[0]):
        retval.append([])
        for y in range(shape[1]):
            retval[-1].append(0)
    return retval


gap_penalty  = - 8 # both for opening and extanding

def match_score(alpha, beta,score_matrix):
    if alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return score_matrix[alpha][beta]

def finalize(align1, align2, score_matrix):
    align1 = align1[::-1]    #reverse sequence 1
    align2 = align2[::-1]    #reverse sequence 2

    i,j = 0,0

    #calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    score = 0
    identity = 0
    for i in range(0,len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            symbol = symbol + align1[i]
            identity = identity + 1
            score += match_score(align1[i], align2[i],score_matrix)

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            score += match_score(align1[i], align2[i],score_matrix)
            symbol += ' '
            found = 0

        #if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '
            score += gap_penalty

    identity = float(identity) / len(align1) * 100

    return identity,score,align1,align2


def alignment(seq1, seq2,score_matrix):
    m, n = len(seq1), len(seq2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = zeros((m+1, n+1))      # the DP table

    # Calculate DP table
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[i-1], seq2[j-1],score_matrix)
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: # end toching the top or the left edge
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + match_score(seq1[i-1], seq2[j-1],score_matrix):
            align1 += seq1[i-1]
            align2 += seq2[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq1[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq2[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq1[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq2[j-1]
        j -= 1

    identity,score,alignment1,alignment2 = finalize(align1, align2,score_matrix)
    return identity,score,alignment1,alignment2


#1987 compute Srand
def scorerand(aftAlignmentSeq1,aftAlignmentSeq2,s_matrix):
    L = len(aftAlignmentSeq1)
    gapnum1 = aftAlignmentSeq1.count('-')
    gapnum2 = aftAlignmentSeq2.count('-')
    gapnum = gapnum1+gapnum2
    set1 = set(aftAlignmentSeq1)-set('-')
    set2 = set(aftAlignmentSeq2)-set('-')
    sum = 0
    for i in set1:
        N_i = aftAlignmentSeq1.count(i)
        for j in set2:
            N_j = aftAlignmentSeq2.count(j)
            sum = sum + s_matrix[i][j]*N_i*N_j
    sum = 1/L*sum+gapnum*gap_penalty
    return sum

# 1987
def distance(seq1,seq2,s_matrix):
    identity1, score1,_,_ = alignment(seq1,seq1,s_matrix)
    identity2, score2,_,_ = alignment(seq2, seq2, s_matrix)
    identity12, score12,aftAlign1,aftAlign2 = alignment(seq1, seq2, s_matrix)
    smax = (score1+score2)/2
    srand = scorerand(aftAlign1,aftAlign2,s_matrix)
    dist = - math.log((score12-srand)/(smax-srand))*100
    return dist

# from guide tree to build multiple alignment
def creatIndexgap(sequence):
    index = {}
    i = 0
    sum = 0
    for ele in sequence:
        if ele == '-':
            sum+=1
        if ele != '-':
            index[i] = sum
            i+=1
            sum = 0
    index[i] = sum
    return index

def inserIndexgap(original,cur):
    tmp = {}
    for k in original:
        tmp[k] = cur[k]-original[k]
    return tmp

def insertGap(sequence,indexgap):
    seq = ""
    i = 0
    for ele in sequence:
         if ele == '-':
            seq += ele
         if ele != '-':
            seq += indexgap[i]*'-'
            seq += ele
            i += 1
    seq += indexgap[i]*'-'
    return seq


# A :[1,2,3]  B: [4,5,6]    对应的seqs的序列会变化
def profileAlign(seqs,A,B,score_matrix,s_matrix):
    maxScore = 0
    for i in A:
        for j in B:
           if  maxScore<score_matrix[i][j]:
               maxScore = score_matrix[i][j]
               x = i
               y = j
    _,_,align1,align2 = alignment(seqs[x],seqs[y],s_matrix)
    indexnow1 = creatIndexgap(align1)
    indexbefore1 = creatIndexgap(seqs[x])
    # insert gap into seqAs
    indexgapA = inserIndexgap(indexbefore1,indexnow1)
    for ka in A:
        seqs[ka] = insertGap(seqs[ka],indexgapA)


    # insert gap into seqBs
    indexnow2 = creatIndexgap(align2)
    indexbefore2 = creatIndexgap(seqs[y])
    indexgapB = inserIndexgap(indexbefore2, indexnow2)
    for kb in B:
        seqs[kb] = insertGap(seqs[kb], indexgapB)
    return seqs


def MSA(guide_tree,seqs,score_matrix,s_matrix):
    ltree = len(guide_tree)
    for k in guide_tree:
        if guide_tree[k][0] == -1:
            continue
        left = guide_tree[k][0]
        right = guide_tree[k][1]
        seqs = profileAlign(seqs,guide_tree[left][3],guide_tree[right][3],score_matrix,s_matrix)

    file = "MSAResult.txt"
    i = 0
    with open(file, "w") as f:
        for ele in seqs:
            f.write("sequence%6d:"%i)
            f.write(ele)
            f.write('\n')
            i+=1
    return seqs

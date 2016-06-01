# 构建 UPGMA guide tree

# test_map = {
#     (0, 0): 0,
#     (0, 1): 6,
#     (0, 2): 4,
#     (0, 3): 2,
#     (1, 0): 6,
#     (1, 1): 0,
#     (1, 2): 1,
#     (1, 3): 0.5,
#     (2, 0): 4,
#     (2, 1): 1,
#     (2, 2): 0,
#     (2, 3): 0.5,
#     (3, 0): 2,
#     (3, 1): 0.5,
#     (3, 2): 0.5,
#     (3, 3): 0
# }

def closest_cluster(dis_map, m, n):
    best_distance = float('inf')
    best_cluster = (-1, -1)
    for key in dis_map:
        (i, j) = key
        if (i == j):
            continue

        dist = dis_map[(i, j)]
        if dist < best_distance:
            best_distance = dist
            best_cluster = (i, j)
    return best_cluster, best_distance


def upgma(dis_map, m, n):
    upgma_tree = {}
    clusters = list(range(n))
    for i in clusters:
        upgma_tree[i] = (-1, -1, 0, [i])

    current_map = dis_map
    new_cluster = n
    while len(clusters) > 2:
        (i, j), distance = closest_cluster(current_map, m, n)
        i_size = len(upgma_tree[i][3])
        j_size = len(upgma_tree[j][3])
        upgma_tree[new_cluster] = (i, j, distance , upgma_tree[i][3] + upgma_tree[j][3])

        clusters.remove(i)
        clusters.remove(j)

        new_dis_map = {}
        for a in clusters:
            for b in clusters:
                new_dis_map[(a, b)] = current_map[(a, b)]

        for a in clusters:
            new_dis_map[(new_cluster, a)] = (current_map[(i, a)] * i_size + current_map[(j, a)] * j_size) / float(
                i_size + j_size)
            new_dis_map[(a, new_cluster)] = new_dis_map[(new_cluster, a)]

        new_dis_map[(new_cluster, new_cluster)] = 0

        clusters.append(new_cluster)
        current_map = new_dis_map
        new_cluster += 1

    (i, j), distance = closest_cluster(current_map, m, n)
    i_size = upgma_tree[i][3]
    j_size = upgma_tree[j][3]
    upgma_tree[new_cluster] = (i, j, distance , i_size + j_size)
    print(upgma_tree)
    return upgma_tree





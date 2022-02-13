f = open('citation-raw.txt')
lines = f.readlines()
id2authors = {}

flag = False

for line in lines:
    line = line.strip('\n')
    if line.startswith("*Vertices"):
        flag = True
        continue
    if line.startswith("*Edges"):
        break
    if flag:
        lst = line.split('\t')
        paperid = lst[0]
        authors = lst[4]
        author_lst = authors.split(',')
        id2authors[paperid] = author_lst

adj_lst_out = {}
edge_flag = False
name_idx = 0
id2name = {}

# 39	98 means paper 39 cites paper 98

for line in lines:
    line = line.strip('\n')
    if line.startswith("*Edges"):
        edge_flag = True
        continue
    if not edge_flag:
        continue
    
    source = line.split('\t')[0]
    target = line.split('\t')[1].split()[0]

    for author in id2authors[source]:
        if author not in id2name:
            # assign a nameID
            id2name[author] = name_idx
            name_idx += 1

    for author in id2authors[target]:
        if author not in id2name:
            # assign a nameID
            id2name[author] = name_idx
            name_idx += 1
    
    for s_author in id2authors[source]:
        for t_author in id2authors[target]:
            s_id = id2name[s_author]
            t_id = id2name[t_author]
            if s_id == t_id:
                continue
            if s_id not in adj_lst_out:
                adj_lst_out[s_id] = []
            if t_id not in adj_lst_out[s_id]:
                adj_lst_out[s_id].append(t_id)

with open('peopleID_d', 'w') as f1:
    for key in id2name:
        f1.write(str(id2name[key]))
        f1.write('\t')
        f1.write(key)
        f1.write('\n')


total_edges = 0

with open('cite_d', 'w') as f1:
    for key in range(0, name_idx):
        if key in adj_lst_out:
            adj_lst_out[key].sort()
            for t in range(len(adj_lst_out[key])):
                total_edges += 1
                f1.write(str(adj_lst_out[key][t]))
                if t != len(adj_lst_out[key])-1:
                    f1.write(' ')
        f1.write('\n')

print("# of authors: ", name_idx)
print("# of edges: ", total_edges)

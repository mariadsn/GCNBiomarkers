# crear ficheros para clusters 

lClusterG = {}
clusterD = open('/home/maria/Documentos/PhD_bired/analysis_Gene_data/tumorProstate_DEGs_table_010_symbolthreshold_0.8_glayCluster_hub.csv','r')
for line in clusterD:
    line = line.strip()
    if not line.startswith('"__glayCluster"'):
        spl = line.split(',')
        clstr = spl[0].replace('"','')
        name = spl[1].replace('"','')
        if not clstr in lClusterG:
            lClusterG[clstr] = []
        if not name in lClusterG[clstr]:
            lClusterG[clstr].append(name)
clusterD.close()
Good_clust = []
for n in lClusterG:
    if len(lClusterG[n]) >10:
        Good_clust.append(n)
print(len(Good_clust))

for elem in lClusterG:
    if len(lClusterG[elem])>10:
        f = open('/home/maria/Documentos/PhD_bired/analysis_Gene_data/tumor/Cluster_' + str(elem) + '.txt','w')
        for l in lClusterG[elem]:
            f.write(l + '\n')
        f.close()

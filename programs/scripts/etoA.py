import os

vertexNum = 410237
edgeNum = 3000000
edgeList = [[0,[-3]]]

source = "amazon0505_adj.tsv"

  f = open(source, "r")
l = f.readlines()
  l2 = [line.rstrip('\n') for line in l]

  for i in range(1,vertexNum+1):
    edgeList.append([i,[]])

    for line in l2:
  graph_Local = [line.split(" ")[0],line.split(" ")[1]]
  edgeList[int(graph_Local[0])][1].append(int(graph_Local[1]))
edgeList[int(graph_Local[1])][1].append(int(graph_Local[0]))


  with open('amazon0505_adj.tsv_Adj','w') as eFile:
  for item in edgeList:
  eFile.write("%s\n" % item[1])

eFile.close()

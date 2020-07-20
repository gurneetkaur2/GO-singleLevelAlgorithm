import sys
"""
usage:
python adjacency-list-to-edges.py [adjacency-list-file] [edge-list-file]
example:
this program will help convert this:
[adjacency-list-file]
1 2 3 4 5
2 1 3 4
3 1 2 4 6 7
into this:
[edge-list-file]
1 2 1
1 3 1
1 4 1
1 5 1
2 1 1
2 3 1
2 4 1
3 1 1
3 2 1
3 4 1
3 6 1
3 7 1
"""

adjacency_list_filename = sys.argv[1]
edge_list_filename = sys.argv[2]

edge_list = []
with open(adjacency_list_filename, 'r') as f:
  for line in f:
    line = line.rstrip('\n').split(' ')
    source = line[0]
  for target in line[1:]:
    edge_list.append("%s %s 1" % (source, target))

with open(edge_list_filename, 'w') as f:
  f.write('%s\n' % ('\n'.join(edge_list)))

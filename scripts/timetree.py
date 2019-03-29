#!/usr/bin/env python3

# The branch length is retrieved from http://www.timetree.org/

from ete3 import Tree

t = Tree()

n1 = t.add_child(dist = 12.9, name = "root")
t.add_child(name = "rheMac", dist = 28.1)

n2 = n1.add_child(dist = 6.59, name = "hg.panTro.ponAbe")
n1.add_child(name = "ponAbe", dist = 15.2)

n3 = n2.add_child(dist = 2.21, name = "hg.panTro")
n2.add_child(name = "gorGor", dist = 8.61)

n3.add_child(name = "hg", dist = 6.4)
n3.add_child(name = "panTro", dist = 6.4)

#print(t)
print(t.write(features=[]))

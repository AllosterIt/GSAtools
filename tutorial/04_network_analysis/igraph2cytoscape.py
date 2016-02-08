#!/usr/bin/env python
'''
$Id: igraph2cytoscape.py 1391 2013-04-16 00:48:16Z apandini $

Copyright (C) 2010-2012  Alessandro Pandini

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import sys

xScale = 1000
yScale = -1000
wScale = 20

def main(args):
        fin = open(args[1])
        lines = fin.readlines()
        fin.close()
        nodeData = dict()
        for line in lines:
                nodeName, x, y = line.split()
                x = float(x)
                y = float(y)
                nodeData[nodeName] = {'x': x, 'y': y}
        fin = open(args[2])
        lines = fin.readlines()
        fin.close()
        fin = open(args[3])
        colorLines = fin.readlines()
        fin.close()
        colorDict = dict()
        for line in colorLines:
                sNodeLabel, eNodeLabel, color = line.split()
                sNode = int(sNodeLabel.split('_')[1]) - 1
                eNode = int(eNodeLabel.split('_')[1]) - 1
                colorDict[(sNode,eNode)] = color
        fout = open(args[2][:-4] + '.cys.gml', 'w')
        for line in lines:
                fout.write(line)
                if line[4:8] == 'name':
                        thisNode = (line.split()[1]).replace('\"', '')
                        fout.write("    graphics\n    [\n        x %8.3f\n        y %8.3f\n        type \"oval\"\n    ]\n" % (nodeData[thisNode]['x'] * xScale, nodeData[thisNode]['y'] * yScale))
                        fout.write("    label \"%s\"\n" % (thisNode))
                        fout.write("")
                if line[4:10] == 'source':
                        sourceId = int(line.split()[1])
                if line[4:10] == 'target':
                        targetId = int(line.split()[1])
                if line[4:10] == 'weight':
                        eWeight = float(line.split()[1])
                        if colorDict.has_key((sourceId, targetId)):
                                eColor = colorDict[(sourceId, targetId)].lower()
                        else:
                                eColor = colorDict[(targetId, sourceId)].lower()
                        fout.write("    graphics\n    [\n        width %8.3f\n        fill \"%s\"\n    ]\n" % ((1 - eWeight) * wScale, eColor))
                        
        fout.close()

if __name__ == '__main__':
        if len(sys.argv) != 4:
                print "Usage: %s <pcalayoutTable> <gmlfile> <colorData>" % sys.argv[0]
                exit()
        main(sys.argv)





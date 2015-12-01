from pprint import PrettyPrinter
from math import sqrt
import networkx as nx
import pygraphviz as pgv
import matplotlib.pyplot as plt
from subprocess import check_output
from sys import argv

# This should be the final version of the concept drawing function. I'm going to leverage
# several different libraries in the name of not messing anything up myself, and to make
# this whole thing slightly more readable. The general idea is as follows:
# 1. Look ahead to the whole lession. Build a graph out of it. Lay it out using a modified sugiyama
#   layout using graphviz. Record the positions of the nodes (maybe edges too?) in the associated
#   pos dict.
#
# 2. Enter listening mode. The program will wait for input, whether it be from the agent or
#     from the user. If agent input, retrieve the location of the soon-to-be displayed node
#     and draw it using matplotlib. If user input, do nothing unless it overlaps, in which case
#     use the force-scan algorithm
# 
# 3. Either way, add a history item after number 2 is completed. Try and keep it uniform, i
#     suppose


GVPIX = 1.388

class ForceScan:
    def __init__(self, graph, pos):
        self._defaultWidth = 0
        self._defaultHeight= 0
        self._nodesep = 20
        self._pos = {}
        self._graph = nx.DiGraph()
        self._initialLayout(graph, pos)
        self._forceScan()

    def _returnNewPos(self):
        pass

    def _initialLayout(self, graph, posd):
        #given a posd dict, create a forcescan class with nodes in their specified positions
        for node in posd:
            width = posd[node]['width']
            height = posd[node]['height']
            center = [ posd[node]['pos'][0] + int(width/2), posd[node]['pos'][1] + int(height/2) ]
            self._pos[node] = {'width': width, 'height': height, 'center':center}

        self._graph = graph

    def _addNode(self, node):
        #add a node to node dict. Here, node is a dict with initial position information
        #ForceScan._addNode(self, {'id'=name, 'pos'=(x,y)})
        #self.graph.add_node(
        pass

    def _removeNode(self, node):
        pass

    def _addEdge(self, node):
        #adds an edge to the edge dict
        pass
    
    def _unitVector(self, node, target):
        #reutrns a 2-element list representing the unit vector in the direnction u -> v
        np = self._pos[node]['center']
        tp = self._pos[target]['center']
        d = self._euclideanDistance(node, target)
        v = [tp[0] - np[0], tp[1] - np[1]]
        uv = [x / d for x in v]
        return uv

    def _euclideanDistance(self, node, target):
        #returns a scalar that is the distance in pixels between node and target
        np = self._pos[node]['center']
        tp = self._pos[target]['center']
        d = sqrt( (np[0] - tp[0])**2 + (np[1] - tp[1])**2)
        return d

    def _calculateForce(self, node, target):
        #returns a 2-element list representing the projections of the forces in the x and y directions
        uv = self._unitVector(node, target)
        duv = self._euclideanDistance(node, target)
        
        u = self._pos[node]
        v = self._pos[target]

        w_u = u['width']
        w_v = v['width']

        h_u = u['height']
        h_v = v['height']

        k_uv_x = (w_u + w_v) / 2 
        k_uv_y = (h_u + h_v) / 2

        f_uv = [((k_uv_x - duv) * uv[0]), ((k_uv_y - duv) * uv[1])]

        return f_uv

    def _forceScan(self):
        #updates all position dict entries
        g = self._graph.succ
        p = self._pos

        V = len(g)

        xs = sorted(p, key=lambda x: p[x]['center'][0])
        ys = sorted(p, key=lambda y: p[y]['center'][1])

        # FORCE-SCAN PSUEDOCODE
        # i <- 1
        # while i < |V| do
        #   suppose x_i = x_i+1 = ... = x_k;
        #   delta <- max ([ f^x(v_m, v_j) for i <= m  <= k <=j <=|V|]
        #   for (j=k+1, |V|) do: x_vj <- x_vj + delta
        #   i <- k+1

        #Horizontal Scan
        i = 0
        while i < V:
            #nodes after i whose centers have the same x value as i's
            xi = [x for x in xs[i:] if self._pos[x]['center'][0] == self._pos[xs[i]]['center'][0]]
            #print xi
            k = len(xi)
            try:
                delta = max([self._calculateForce(x, y)[0] for x in xi for y in xs[i + k:]])
            except ValueError:
                delta = 0

            for x in xs[i+k:]:
                p[x]['center'][0] = int(p[x]['center'][0] + delta)

            i += k

        #Vertical Scan
        i = 0
        while i < V:
            #nodes after i whose centers have the same x value as i's
            yi = [y for y in ys[i:] if self._pos[y]['center'][1] == self._pos[ys[i]]['center'][1]]
            k = len(yi)
            try:
                delta = max([self._calculateForce(x, y)[1] for x in yi for y in ys[i + k:]])
            except ValueError:
                delta = 0

            for y in ys[i+k:]:
                p[y]['center'][1] = int(p[y]['center'][1]+ delta)

            i += k

class ConceptDrawer:
    def __init__(self, graph):
        self._dataGraph = nx.DiGraph()
        self._displayGraph = nx.DiGraph()
        self._history = []
        self._compound = {'pred': {},'succ': {}}
        self._dotfile = ''
        self._dotnode = {}
        self._minNodeDist = 20
        self._pos = {}
        self._processBaseGraph(graph)

    def _buildInitialGraph(self):
        pass

    def _nodeDist(node1, node2):
        pass

    def _findInitialPos(self):
        pass

    def _draw(self):
        pass

    def _updateDrawGraph(self):
        pass

    def _addToHistory(self, histObject):
        pass

    def _forceScan(self):
        pass

    def _toPyplot(self, pos):
        #takes a pos dict in pixels and returns it in an array, just like matplotlib likes it.
        posnx = {}
        for node in pos:
            posnx[node] = [pos[node]['pos'][0],pos[node]['pos'][1]]

        return posnx

    def drawplt(self, pos):
        #G is actual semantic graph
        #pos is the position dict (posnx, use the _toPyplot)
        nx.draw_networkx(self._dataGraph, pos)
        plt.show()
        pass

    def findCompound(self, graph):
        #takes care of all the pretty stuffs for the base graph - namely, compound edges
        #and horizontal relationships
        grouped = []
        groups = []
        groupsh = []            #horizontal group
        labels= []
        for node in graph.succ:
            #new strategy: iterate through all edge labels. find all nodes with edges with that label.
            #check to see if others are compound
            mlabel = [graph.succ[node][value]['label'] for value in graph.succ[node]]
            labels.extend(mlabel)

        for node in graph.pred:
            mlabel = [graph.pred[node][value]['label'] for value in graph.pred[node]]
            labels.extend(mlabel)

        labels = list(set(labels))

        #now that we have the labels in question:
        for label in labels:
            for node in graph.succ:
                #for each label, make key, value pairs for nodes
                sharedlabels = [x for x, y in graph.succ[node].iteritems() if y['label'] == label]
                if len(sharedlabels) > 1:
                    #check if all horiz attributes are True
                    if False not in [graph.succ[node][x]['horiz'] == True for x in sharedlabels]:
                        groupsh.append(sharedlabels)
                    else:
                        groups.append(sharedlabels)
                    for x in sharedlabels:
                        grouped.append(x)

            for node in graph.pred:
                #for each label, make key, value pairs for nodes
                sharedlabels = [x for x, y in graph.pred[node].iteritems() if y['label'] == label]
                if len(sharedlabels) > 1:
                    if False not in [graph.pred[node][x]['horiz'] == True for x in sharedlabels]:
                        groupsh.append(sharedlabels)
                    else: 
                        groups.append(sharedlabels)
                    for x in sharedlabels:
                        grouped.append(x)


        for node in graph.succ:
            if node not in grouped:
                grouped.append(node)
                groups.append([node])

        return groups, groupsh
    
    def addNodeToDotfile(self, node, *args):
        #node is actually a list - we will do different things based on whether this list is of
        #length 1 or not.

        #only called during creation of the dotfile.

        #create property dict for all nodes: keys are node names, values what graphviz will be using
        #to id nodes in record nodes when making edges
        #while iterating through the groups, add all nodes to the dotfile

        #this is where the ids should be stripped of spaces
        if len(node) == 1:
            nodeit = node[0]
            nodestrip = nodeit.replace(' ','').replace('-','')
            self._dotfile += 'node [shape=box];\n'
            self._dotnode[nodeit] = {'id' : nodestrip, 'label':nodeit}  
            self._dotfile += nodestrip + ' [label="' + nodeit + '"];\n'

        elif 'horizontal' in args:
            self._dotfile += 'node [shape=record];\n'
            recordlabel = ''.join(node).replace(' ','').replace('-','')
            self._dotfile += recordlabel + ' [label="{'
            for num, nodeit in enumerate(node):
                nodestrip = nodeit.replace(' ','').replace('-','')
                self._dotnode[nodeit] = {'id' : recordlabel + ':' + nodestrip, 'label': nodeit}
                self._dotfile += '<' + nodestrip + '>' + nodeit
                if len(node) == num + 1:
                    self._dotfile += '}"];\n'
                else:
                    self._dotfile += '|'
        else: 
            self._dotfile += 'node [shape=record];\n'
            recordlabel = ''.join(node).replace(' ','').replace('-','')
            self._dotfile += recordlabel + ' [label="'
            for num, nodeit in enumerate(node):
                nodestrip = nodeit.replace(' ','').replace('-','')
                self._dotnode[nodeit] = {'id' : recordlabel + ':' + nodestrip, 'label': nodeit}
                self._dotfile += '<' + nodestrip + '>' + nodeit
                if len(node) == num + 1:
                    self._dotfile += '"];\n'
                else:
                    self._dotfile += '|'

    def formatBaseGraph(self, graph):
        #this method that takes groups of compound predicates and the initial data graph and
        #makes a dotfile for layout with graphviz. Not quite worried about edge nodes yet - 
        #will work on that for the rest of the month. 

        #iterate through groups list
        groups, groupsh = self.findCompound(graph)

        for group in groups:
            self.addNodeToDotfile(group)

        for group in groupsh:
            self.addNodeToDotfile(group, 'horizontal')
            
        #next, iterate through the successor dict and add all edges to the dotfile, using the
        #property dict to ensure proper edge labels
        
        for source in graph.succ:
            for target in graph.succ[source]:
                self._dotfile += self._dotnode[source]['id'] + '->'
                self._dotfile += self._dotnode[target]['id'] +';\n'
                    
        self._dotfile = 'digraph G {\n'+ self._dotfile + '}'


    def _posFromParsedDot(self):
        #pos is the dict of positions we're going to return. It'll look like this:
        #pos = {node: ("width": width, "height": height, "pos": pos} for however many
        #nodes there may be
        posd = {}

        dotfile = self._dotout
        dotfile = dotfile.replace('\n','').replace('\t','')
        entries = dotfile.split(';')
        entries2 = [entry for entry in entries if '->' not in entry and 'pos=' in entry]

        #record nodes' positions are defined in "GV pixels", from which I can tell are
        #approximately 1/1.388 the size of the pixels on my screen. Height and width are
        #given in inches, and can be converted to GV pixels by multiplying by 96.
        #we want to define a bottom-left hand corner (assuming we're in Quadrant I) 
        for entry in entries2:
            if 'rects=' in entry: 
                #record nodes' rects are so far unknwon. We use their relative widths/heights to
                #figure out how many GV pixels to assign to each record

                #split on the '"'
                names = []
                horizNames = []
                pos = (0,0)
                horiz = False
                ids = []
                width = 0
                height = 0
                current = entry.split('"')
                for i, split in enumerate(current):
                    if 'pos=' in split: pos = tuple(current[i+1].split(","))
                    if 'label=' in split: 
                        names = current[i+1].split("|")
                        if True in ['{' in x or '}' in x for x in names]:
                            horizNames.extend(names)
                    if 'rects=' in split: rects = current[i+1].split(" ")
                    splits = split.split("=")
                    for j, subsplit in enumerate(splits):
                        if 'height' in subsplit: 
                            height = splits[j+1].split(",")[0]
                        if 'width' in subsplit: 
                            width = splits[j+1].strip("]")

                if names:
                    horiz = False
                    for name in names:
                        idtoappend = name.split('>')[1]
                        if name in horizNames:
                            idtoappend = idtoappend.replace('{','').replace('}','')
                            horiz = True
                        ids.append(idtoappend)

                    #begin arduous hacking logic. Going to keep everything in "real pixels"
                    center = tuple([int(float(coord) * GVPIX) for coord in pos])
                    width = int(float(width)*96)
                    height = int(float(height)*96)

                    #get realtive widths and heights
                    ys = [float(x) for x in ','.join(rects).split(',')[1::2]]
                    xs = [float(x) for x in ','.join(rects).split(',')[0::2]]
                    maxw = abs(max(xs) - min(xs))
                    maxh = abs(max(ys) - min(ys))

                    relws = [(xs[i+1] - xs[i]) / maxw for i in range(0, len(xs), 2)]
                    relhs = [(ys[i+1] - ys[i]) / maxh for i in range(0, len(ys), 2)]
                    widths = [int(width * x) for x in relws]
                    heights = [int(height * x) for x in relhs]

                    #get bottom-left corner
                    corner = (center[0] - width/2, center[1] - height/2)

                    #this whole block, beginning with "if names:", does positioning in chunks
                    #Therefore we can assume that all nodes in "ids" are all laid out the
                    #same way (horizontally, vertically) and change the positioning chunk as
                    #we need it

                    if not horiz:
                        for i, myid in enumerate(ids):
                            if i > 0:
                                sumw = sum([widths[j] for j in range(i)])
                                posd[myid] = {'width': widths[i], 'height':height, 'pos': (
                                        corner[0] + sumw, corner[1])}
                            else:
                                posd[myid] = {'width': widths[i], 'height':height, 'pos': corner}

                    #horizontal arrangement
                    else:
                        for i, myid in enumerate(ids):
                            if i > 0:
                                sumh = sum([heights[j] for j in range(i)])
                                posd[myid] = {'width': width, 'height':heights[i], 'pos': (
                                        corner[0], corner[1] + sumh)}
                            else:
                                posd[myid] = {'width': width, 'height':heights[i], 'pos': corner}
                
            #i.e. shape is not record
            else: 
                current = entry.split("=")
                ids = entry.split('[')[0].strip()
                ids = self._nameFromDotLabel(ids)
                #fix space in normal node issue here
                for i, split in enumerate(current):
                    if "height" in split: height = int(float(current[i+1].split(',')[0] ) * 96)
                    if "width" in split: width = int(float(current[i+1].strip(']')) * 96)
                    if "pos" in split: pos = tuple(current[i+1].split('"')[1].split(','))
                pos = tuple([int(float(x)) for x in pos])
                posd[ids] = {'width':width, 'height':height, 'pos':pos}
                
        #pp = PrettyPrinter()
        #pp.pprint(posd)
        return posd

    def _nameFromDotLabel(self, nodeid):
        for node, item in self._dotnode.items():
            if nodeid == item['id']: return item['label']
        return nodeid

    def addNode(self, node, pos):
        self._dataGraph.add_node(node)
        self._pos[node] = {'pos':tuple(pos), 'height':50,'width':90}
        self.plotGraph('forcescan')
        pass
        
    def plotGraph(self, *args):
        if 'forcescan' in args:
            d = ForceScan(self._dataGraph, self._pos)
            self._updatePosDict(d)

        posd = {}
        for node, pos in self._pos.items():
            posd[node] = pos['pos']
        nx.draw(self._dataGraph, pos=posd, node_shape='s', node_size=1200) 
        plt.show()

    def listen(self):
        while True:
            node = input('add a node: name?\n')
            posx = input('add a node: x position?\n')
            posy = input('add a node, y position?\n')
            print 'added node', node, 'at', posx, posy
            #add the node to the conceptDrawer class with left-corner notation and corresponding w/h 
            #create apply force-scan to the new pos dict in conceptdrawer
            self.addNode(node, (posx,posy))

    def _updatePosDict(self, forcescan):
        #take updated dict from force scan and update the drawgraph (self.graph) with new values
        #to get left bottom corner, subtract half the width from x and half the height from y
        fd = forcescan._pos
        sd = self._pos
        for node in sd:
            width = sd[node]['width']
            height = sd[node]['height']
            sd[node]['pos'] = (fd[node]['center'][0] - int(width/2), fd[node]['center'][1] - int(height/2))

    def _processBaseGraph(self, graph):
        self._dataGraph = graph.copy()
        self.formatBaseGraph(graph)

        #write the dotfile and run it with graphviz
        with open('./dotout.txt', 'w+') as f:
            f.write(self._dotfile)

        #actually run dot
        self._dotout = check_output(['dot', './dotout.txt'])


        #write the output of dot for checkin'
        #with open('./conceptDraw.dot','w') as f:
        #    f.write(self._dotout)


        #plot the graph pre-forcescan
        self._pos = self._posFromParsedDot()
        nxpos = self._toPyplot(self._pos)

        self.drawplt(nxpos)

        #create a forescan class and feed it the position dict we harvested from DOT
        d = ForceScan(graph, self._pos)
        self._updatePosDict(d)

        #plot the graph again
        nxpos = self._toPyplot(self._pos)
        self.drawplt(nxpos)

if __name__== '__main__':
    inputgraph = './lessons/14-6-7.py'
    graphInput = nx.DiGraph()
    grGlobals = {'graph': graphInput}
    execfile(inputgraph, grGlobals)
    d = ConceptDrawer(graphInput)
    #d.listen()

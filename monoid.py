import numpy as np

ex = np.array([1,0])
ey = np.array([0,1])

class word:
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix

    def __str__(self):
        return self.prefix + "W" + str(self.suffix)

    def __eq__(self, other):
        return self.prefix == other.prefix and self.suffix == other.suffix
    
    def latex(self, inline = False):
        if self.suffix > 1:
            ret = self.prefix + "W_{" + str(self.suffix) + "}"
        elif self.suffix == 1:
            ret = self.prefix + "V"
        else:
            ret = self.prefix + "U"
        if inline:
            return ret
        else:
            return "$" + ret + "$"

class signed_letter:
    def __init__(self, letter, inv):
        self.letter = letter
        self.inv = inv


class signed_word:
    def __init__(self, letters, signed = False):
        if not signed:
            assert all(not isinstance(l, signed_letter) for l in letters) 
            self.letters = [signed_letter(l, False) for l in letters]
        else:
            self.letters = letters
    #make the signed word subscriptable
    def __getitem__(self, val):
        if isinstance(val, int):
            return self.letters[val]
        else:
            return signed_word(self.letters[val], signed = True)
        
    def __str__(self):
        ret = ""
        for l in self.letters:
            if l.inv:
                ret += "!" + str(l.letter)
            else:
                ret += str(l.letter)
        return ret
    #truncate
    
    def inverted(self):
        return signed_word([signed_letter(l.letter, not l.inv) for l in self.letters[::-1]], signed = True)

    def __add__(self, other):
        return signed_word(self.letters + other.letters, signed = True)

    def latex(self, inline = False):
        ret = ""
        for l in self.letters:
            if l.inv:
                ret += "\\overline{" + str(l.letter) + "}"
            else:
                ret += str(l.letter)
        if inline:
            return ret
        else:
            return "$" + ret + "$"

def rr_start(u,v):
    return signed_word(u).inverted()+ signed_word(v)

class equation:
    def __init__(self, w1, w2):
        self.w1 = w1
        self.w2 = w2

    def __str__(self):
        return str(self.w1) + " = " + str(self.w2)
    
    def latex(self):
        return "$" + self.w1.latex(inline = True) + " \\meq " + self.w2.latex(inline = True) + "$"
    

class trequation:
    def __init__(self, w1, wmid, w2):
        self.w1 = w1
        self.wmid = wmid
        self.w2 = w2

    def __str__(self):
        return str(self.w1) + " = " + str(self.wmid) + " = " + str(self.w2)
    
    def latex(self):
        return "$"+self.w1.latex(inline = True) + " \\meq " + self.wmid.latex(inline = True) + " \\meq " + self.w2.latex(inline = True) + "$"
         
class vertex:
    def __init__(self, equation, color = "black", depth = 0):
        self.equation = equation
        self.up_edge = []
        self.bot_edge = []
        self.depth = depth
        self.color = color
        self.position = (0,0) #position in the graph

class edge:
    def __init__(self, o, t, label, color = "black"):
        self.o = o
        self.t = t
        self.label = label
        o.bot_edge.append(self)
        t.up_edge.append(self)
        self.color = color



class lcm:
    def __init__(self, rco, u, v, show = False):
        self.rco = rco
        self.show = show
        Z1 = word(u, 0)
        Z2 = word(v, 1)
        self.count = 1
        self.free = {} #dictionary of free variables as key and the corresponding vertex as value
        self.graph = [[vertex(equation(Z1,Z2))]]
        self.U = word("", 0)
        self.V = word("", 1)
        self.limit_depth = None
        self.edges = [] 
        self.converged = None

    def add_vertex(self, v, left = False, test = True):
        #call after adding the edge
        if v.depth >= len(self.graph):
            self.graph.append([])
        if left:
            self.graph[v.depth].insert(0, v)
        else:
            self.graph[v.depth].append(v)
        if test:
            self.test_stopping(v)


    def counter(self):
        self.count += 1
        return self.count

    def test_stopping(self, v):
        #to stop the algorithm in the non spherical case
        #test if the pair pref1,pref2 already appeard in a parent vertex of v:
        if self.limit_depth: #if limit_depth is already set, the stopping is already programmed
            return
        #test if the pair pref1,pref2 already appeard in a parent vertex of v:
        pref1 = v.equation.w1.prefix
        pref2 = v.equation.w2.prefix
        ancestors = [[v]]
        edge_list = []
        while ancestors [-1]:
            #if the last layer is empty, stop
            ancestors.append([])
            for a in ancestors[-2]:
                for e in a.up_edge:
                    edge_list.append(e)
                    if e.o not in ancestors[-1]:
                        ancestors[-1].append(e.o)
            if len(ancestors[-1]) == 1:
                apref1 = ancestors[-1][0].equation.w1.prefix
                apref2 = ancestors[-1][0].equation.w2.prefix
            if (apref1 == pref1 and apref2 == pref2) or (apref1 == pref2 and apref2 == pref1):
                #if it is the case, stop developing after the current depth and remember vertex v and its ancestor
                print("stop developing at depth", v.depth, "because of the pair", pref1, pref2)
                self.limit_depth = v.depth
                for e in edge_list:
                    e.color = "red"
                for layer in ancestors:
                    for a in layer:
                        a.color = "red"

    def develop(self, v, depth=0):
        if depth > 20:
            raise ValueError("depth too high")
        if self.limit_depth and v.depth > self.limit_depth:
            return
        if self.show:
            print("developing", v.equation)
        w2 = v.equation.w2
        w1 = v.equation.w1
        if w2.prefix == "":
            n = w2.suffix
            self.new_free(n, v)
                
        elif w1.prefix == "":
            n = w1.suffix
            self.new_free(n, v)
        else:
            #both words have non empty prefix
            x = w1.prefix[0]
            y = w2.prefix[0]
            if x == y:
                #if the first letter are the same, remove it and recurse
                v2 = vertex(equation(word(w1.prefix[1:], w1.suffix), word(w2.prefix[1:], w2.suffix)), depth=v.depth+1)
                #add edge and vertex to the graph
                self.edges.append(edge(v, v2, x))
                self.add_vertex(v2)
                if self.show:
                    print(v2.equation)

                
            else:
                #if the first letter are different, create two new vertices by using the relation 
                n = self.counter()
                pref1 = w1.prefix
                pref2 = w2.prefix
                v2 = vertex(equation(word(pref1[1:], w1.suffix), word(self.rco(x,y), n)), depth=v.depth+1)
                v3 = vertex(equation(word(self.rco(y,x), n), word(pref2[1:], w2.suffix)), depth=v.depth+1)
                if self.show:
                    print(v2.equation, " and ", v3.equation)

                #add edge and vertex to the graph        
                self.edges.append(edge(v, v2, x))
                self.edges.append(edge(v, v3, y))
                self.add_vertex(v2)
                self.add_vertex(v3)
                
    def new_free(self, n, v):
        assert n not in self.free
        self.free[n] = v

    def build_path(self, ver, wor, left = False):
        if self.show: 
            print("growing path from", ver.equation)
        if len(wor.prefix) == 0:
            return
        if len(wor.prefix) == 1:
            if wor.suffix in self.free:
                self.edges.append(edge(ver, self.free[wor.suffix], wor.prefix[0], color = "blue"))
                if self.show:
                    print("free variable", wor.suffix, "found")
            else:
                if not self.converged:
                    new_wor = word("", wor.suffix)
                    new_ver = vertex(new_wor, depth=ver.depth+1, color = "blue")
                    self.edges.append(edge(ver, new_ver, wor.prefix[0], color = "blue"))
                    self.add_vertex(new_ver, left = left, test= False)
                    self.converged = new_ver
                else:
                    self.edges.append(edge(ver, self.converged, wor.prefix[0], color = "blue"))
                    assert self.converged.equation.suffix == wor.suffix 
        else:
            new_wor = word(wor.prefix[1:], wor.suffix)
            new_ver = vertex(new_wor, depth=ver.depth+1, color = "blue")
            self.edges.append(edge(ver, new_ver, wor.prefix[0], color = "blue"))
            self.add_vertex(new_ver, left = left, test= False)
            self.build_path(new_ver, word(wor.prefix[1:], wor.suffix), left=left)

    def fill_U_V(self):
        while self.U.suffix in self.free:
            v = self.free[self.U.suffix]
            wor = v.equation.w2
            self.build_path(v, wor, left = True)
            del self.free[self.U.suffix]
            self.U = word(self.U.prefix + wor.prefix, wor.suffix)
        while self.V.suffix in self.free:
            v = self.free[self.V.suffix]
            wor = v.equation.w1
            self.build_path(v, wor, left = False)
            del self.free[self.V.suffix]
            self.V = word(self.V.prefix + wor.prefix, wor.suffix)
        print("found lcm", self.graph[0][0].equation.w1.prefix + self.U.prefix, "=", self.graph[0][0].equation.w2.prefix + self.V.prefix)
        
    def solve(self, max_depth = 50):
        #print("looking for the lcm of", self.graph[0][0].equation.w1.prefix, "and", self.graph[0][0].equation.w2.prefix, "with the rcomplement matrix", self.rco)
        depth = 0
        while self.graph[depth]:
            if self.limit_depth and depth > self.limit_depth-1:
                if self.show:
                    print("limit depth reached", self.limit_depth)
                break
            if depth == max_depth:
                print("depth too high")
                break
            self.graph.append([])
            for v in self.graph[depth]:
                self.develop(v)
            if self.graph[depth+1]:
                self.graph[depth+1] = self.merging_at_depth(depth+1)
            if self.show:
                print("go to depth", depth+1)
            depth += 1
        
        if not self.limit_depth:
            if depth == max_depth:
                return
            else:
                self.fill_U_V()
                assert self.U.suffix == self.V.suffix
        return self.U.prefix, self.V.prefix



    def merging_at_depth(self, i):
        merged = [self.graph[i][0]]
        if self.show:
            print("merging at depth", i)
        j = 1
        while j < len(self.graph[i]):
            if  (merged[-1].equation.w2 == self.graph[i][j].equation.w1):
                eq = equation(merged[-1].equation.w1, self.graph[i][j].equation.w2)
                new_v = vertex(eq, depth = i)
                while merged[-1].up_edge:
                    e = merged[-1].up_edge.pop()
                    e.o.bot_edge.remove(e)
                    self.edges.append(edge(e.o, new_v, e.label))
                    self.edges.remove(e)
                while self.graph[i][j].up_edge:
                    e = self.graph[i][j].up_edge.pop()
                    e.o.bot_edge.remove(e)
                    self.edges.append(edge(e.o, new_v, e.label))
                    self.edges.remove(e)
                merged[-1] = new_v
            else:
                merged.append(self.graph[i][j])
            j += 1
        return merged

def random_word(n):
    u = ""
    l = np.random.randint(1, n)
    for i in range(l):
        u += chr(np.random.randint(97, 97+n))
    return u

def random_mat(n):
    mat = np.ones((n,n), dtype=int)
    for j in range(1,n):
        for i in range(j):
            mat[i,j] = np.random.randint(2, 10)
            mat[j,i] = mat[i,j]
    return mat 
                
class graph_arrow:
    def __init__(self, letter, start, end):
        self.letter = letter
        self.start = start
        self.end = end

class right_reverse:
    def __init__(self, rco, w):
        self.rco = rco
        self.steplist = [w]
        self.w = w
        self.table = [np.array([0,0])]
        self.arrows = []

    def latex(self):
        content = ""
        for w in self.steplist:
            content += w.latex() + "\\newline\n"
        return content

    def build_graph(self): 
        for i in range(len(self.w.letters)):
            if self.w.letters[i].inv:
                self.table.append(self.table[-1] - ey)
                self.add_arrow(self.w.letters[i].letter, self.table[-1], self.table[-2])
            else:
                self.table.append(self.table[-1] + ex)
                self.add_arrow(self.w.letters[i].letter, self.table[-2], self.table[-1])
        self.w_arrows = self.arrows.copy()
        
    def add_point(self, o):
        for i in range(len(self.table)):
            if o[0] == self.table[i][0] and o[1] == self.table[i][1]:
                return i
        self.table.append(o)
        return len(self.table) - 1
    
    def add_arrow(self, letter, o, t):
        start = self.add_point(o)
        end = self.add_point(t)
        a = graph_arrow(letter, start, end)
        self.arrows.append(a)
        return a

    def step(self, build = False):
        w = self.steplist[-1]
        if len(w.letters) > 100:
            raise ValueError("word too long, please reduce the size of the input words")
        #find a succession of a negative letter followed by a positive letter
        i = 0
        while i+1 < len(w.letters):
            #find a succession of a negative letter followed by a positive letter
            s = w[i]
            t = w[i+1]
            if s.inv and not t.inv:
                #if s is negative and t is positive, we can apply the relation
                s1 = signed_word(self.rco(t.letter, s.letter))
                t1 = signed_word(self.rco(s.letter, t.letter))
                w = w[:i] + t1 + s1.inverted() + w[i+2:]
                if build:
                    n_horiz = len(t1.letters)
                    n_vert = len(s1.letters)
                    s_a = self.w_arrows[i]
                    sao = self.table[s_a.end]
                    sat = self.table[s_a.start]
                    t_a = self.w_arrows[i+1]
                    tao = self.table[t_a.end]
                    tat = self.table[t_a.start]
                    assert (tat == sat).all()
                    sqe = (int(tao[0]), int(sao[1]))
                    for j in range(len(self.table)):
                        pos = self.table[j]
                        if pos[0] >= sqe[0]:
                            pos[0] += n_horiz-1
                        if pos[1] >= sqe[1]:
                            pos[1] += n_vert-1
                    harrows = []
                    for j in range(n_horiz):
                        a = self.add_arrow(t1[j].letter, sao + j*ex, sao + (j+1)*ex)
                        harrows.append(a)

                    varrows = []
                    for j in range(n_vert):
                        a = self.add_arrow(s1[j].letter, tao + j*ey, tao + (j+1)*ey)
                        varrows.append(a)                    
                    self.w_arrows = self.w_arrows[:i] + harrows + varrows[::-1] + self.w_arrows[i+2:]

                self.w = w
                self.steplist.append(self.w)
                return True
            else:
                i += 1 
        return None
    
    def solve(self):
        while self.step()!= None:
            pass
        return self.w

    def solve_and_arrows(self):
        self.build_graph()
        while self.step(build = True) != None:
            pass
        return self.draw_tikz()

    def draw_tikz(self):
        content = "\\begin{tikzpicture}\n"
        for i in range(len(self.table)):
            content += "\\node[scale = 0.5] (P{}) at ({}, {}) {{}};\n".format(i, self.table[i][0], -self.table[i][1])
        for arrow in self.arrows:
            content += self.draw_arrow_tikz(arrow)
        content += "\\end{tikzpicture}\n"
        return content

    def draw_arrow_tikz(self, arrow):
        #draw the arrow
        if arrow.letter == "1":
            arrow.label = ""
            arrow.type = ""
            arrow_param = "[double]"
        else:
            arrow.label = "$"+arrow.letter+"$"
            arrow.type = "-Stealth"
            arrow_param = ""
        return "\\draw[{}] (P{}) edge{} node[auto]{{{}}} (P{});\n".format(arrow.type, arrow.start, arrow_param, arrow.label, arrow.end)

def erco(rco, x, y):
    rr = right_reverse(rco, rr_start(x,y))
    w = rr.solve()
    #return the first letter until the first non poisitive letter
    prefix = ""
    for l in w.letters:
        if l.inv:
            break
        if l.letter == "1":
            continue
        prefix += l.letter
    return signed_word(prefix)
    

    
import numpy as np

def r(a,b,n):
    if n == 0:
        return ""
    else:
        return b + r(b,a,n-1)
    
def cox_abc(ab, bc, ca):
    return np.array([[1,ab,ca],[ab,1,bc],[ca,bc,1]])

def std_dict(n):
    #create a dictionary with n letters
    dict = {}
    for i in range(n):
        dict[chr(i+97)] = i
    return dict

dict_abc = std_dict(3) 


class right_compl:
    def __init__(self, array, dictionary):
        self.array = array
        self.dictionary = dictionary
        
    def __call__(self,x,y):
        #x,y are two letters in the dictionary
        return self.array[(x,y)]
    
    def add_eps(self):
        #copy the dictionary
        dictionary = self.dictionary.copy()
        n = len(dictionary)
        dictionary["1"] = n
        for x in dictionary:
            self.array[(x,x)]="1"
            self.array[(x,"1")]= "1"
            self.array[("1",x)]= x
        self.dictionary = dictionary
    

#a class cox of coxeter matrix inherited from numpy array, assuming the matrix is symmetric
#with integer entries, with diagonal elements 1 and other greater or equal than 2
class coxeter_matrix(right_compl):
    def __init__(self, matrix, dictionary):
        array = {}
        for x in dictionary:
            for y in dictionary:
                array[(x,y)] = r(x,y, matrix[dictionary[x], dictionary[y]] -1)
        super().__init__(array, dictionary)
        self.matrix = matrix
        self.n = self.matrix.shape[0]

    def entry(self, x, y):
        return self.array[self.dictionary[x],self.dictionary[y]]
        

    def __str__(self):
        s = "\n"
        for i in range(self.n):
            s += "["
            for j in range(self.n):
                s += str(self.matrix[i,j]) + ", "
            s+= "],\n"

        return s

class element:
    def __init__(self, cox, s):
        self.cox = cox
        self.lead = s
        self.dev = False
        self.explored = []
        self.to_explore = [s]
        self.cache = s

    def explore(self):
        #find all possible transformations
        while self.to_explore:
            s = self.to_explore.pop()
            #for all pair of distinct letters
            l1 = s[0]
            l2 = None
            k = 1
            i = 1
            while i < len(s):
                if k == 1:
                    if s[i] != l1:
                        l2 = s[i]
                        k = 2
                else:
                    if (k%2 == 0 and s[i] == l1) or (k%2 == 1 and s[i] == l2):
                        k += 1
                    else:
                        k = 2
                        l1 = s[i-1]
                        l2 = s[i]
                if k == self.cox.entry(l1,l2):
                    new_s = s[:i-k+1] + r(l2,l1, k) + s[i+1:]
                    if new_s not in self.explored and new_s not in self.to_explore:
                        self.to_explore.append(new_s)
                    i = i-k+2
                    k = 1
                    l1 = s[i]
                    l2 = None
                i+=1
            self.explored.append(s)    
        self.dev = True

    #equality test between elements
    def __eq__(self, other):
        if self.dev == True:
            if other.lead in self.explored:
                other.explored = self.explored
                other.dev = True
                return True
            else:
                return False
        if other.dev == True:
            return other.__eq__(self)
        self.explore()
        return self.__eq__(other)

    #order test (prefix order)
    def is_prefix(self, other):
        if other.dev == False:
            other.explore()
        for rep in other.explored:
            if rep.startswith(self.lead):
                other.cache = rep
                return True
        return False
    
    #magic method for <
    def __lt__(self, other):
        return self.is_prefix(other)
    
    #subtraction of elements
    def __sub__(self, other):
        if other < self:
            return element(self.cox, self.cache[len(other.lead):])
        else:
            raise ValueError(f"The second element '{other.lead}' is not a prefix of the first one '{self.lead}'")
    
    def __str__(self):
        return self.cache


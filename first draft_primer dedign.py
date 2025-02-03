from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as TM
from Bio.SeqUtils import gc_fraction as gc
from Bio.Restriction import Restriction, RestrictionBatch

# add primer binding in primer design
# or a neat way is to create a general sequence class and inherit from it. 
#hair pin formation visuilazation
# Ta for primer design
class primer:
    
    ''' attributes: TM, GC_content, length, hairpin, binding:pattern matching '''

    def __init__(self,sequence):
        self.sequence=Seq(sequence)
        self.complement_sequence=self.sequence.complement()
        self.length=len(self.sequence)
        self.gc_content=self.calculate_gc()
        self.tm=self.calculate_tm()
        self.patterns_list=self.pattern_complement()

        self.hair_pin, self.double_binding = self.complement_finder()
        return
        


    def calculate_gc(self):
        return round(gc(self.sequence) * 100, 2)
    

    def calculate_tm(self):
        #tm.Tm_NN(my_dna)
        return round(TM.Tm_NN(self.sequence) , 2)
    
    def pattern_complement(self):
        p = []
        complement_seq = self.sequence.complement()
        for i in range(len(complement_seq) - 5 + 1):
            p.append(complement_seq[i:i + 5])
        return p
    
    # RK code is from geeks for geeks
    #@classmethod
    def complement_finder(self):
        d = 256
        q = 101
        M = 5
        N = len(self.sequence)
        i = 0
        j = 0
        p = 0 # hash value for pattern
        t = 0 # hash value for txt
        h = 1
        occurence=[]
        occurence_seq=[]

        # The value of h would be "pow(d, M-1)%q"
        for i in range(M-1):
            h = (h*d) % q

        # Calculate the hash value of pattern and first window
        # of text
        for pa in self.patterns_list:
            for a in range(M):
                p = (d*p + ord(pa[a])) % q
                t = (d*t + ord(self.sequence[a])) % q

            # Slide the pattern over text one by one
            for i in range(N-M+1):
                # Check the hash values of current window of text and
                # pattern if the hash values match then only check
                # for characters one by one
                if p == t:
                    mismatch_count=0
                    # Check for characters one by one
                    for j in range(M):
                        if self.sequence[i+j] != pa[j]:
                            mismatch_count += 1
                        if mismatch_count > 1:
                            break
                        else:
                            j += 1

                    # if p == t and pat[0...M-1] = txt[i, i+1, ...i+M-1]
                    if j == M:
                        occurence_seq.append((pa, self.sequence[i:i+5]))
                        occurence.append((self.patterns_list.index(pa), i))
                #print(occurence,occurence_seq)

                # Calculate hash value for next window of text: Remove
                # leading digit, add trailing digit
                if i < N-M:
                    t = (d*(t-ord(self.sequence[i])*h) + ord(self.sequence[i+M])) % q

                    # We might get negative values of t, converting it to
                    # positive   occurence.append((pattern.index(pa), i))

                    if t < 0:
                        t = t+q
        #processing the matches
        hairpin = []
        double_bind = []
        for matches in occurence:
            dist=max(matches)-(min(matches)+5)
            if dist>=3 and dist<=10:
                hairpin.append(matches)
            else:
                double_bind.append(matches)

        return ("hairpin at positon", hairpin), ("double bind at position", double_bind)


my_dna = Seq("ATTCGGGGAAAAAAATCGAAATGAATAAGCTCCCCCCCGA")
my_dna=primer(my_dna)
print(my_dna.gc_content)
print(my_dna.tm)
print(my_dna.hair_pin)
print(my_dna.double_binding)


#class primer_design:

# primer_1=primer(my_dna)
# print(primer_1.gc_content, primer_1.tm)
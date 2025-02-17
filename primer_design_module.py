from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as TM
from Bio.SeqUtils import gc_fraction as gc
from prettytable import PrettyTable 
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
        self.GC_clamp=self.GC_clamp_finder()
        self.hair_pin, self.double_binding = self.complement_finder()
        return
        
    def GC_clamp_finder(self):
        n= self.length
        A=Seq("A")
        T=Seq("T")
        b=[]
        seqs=(self.sequence[0:2],self.sequence[n-2:n])
        for p in seqs :
            if (A not in p) and (T not in p) :
                b.append(1)
            else:
                b.append(0)

        if b==[1,1]:
            result="GC clamp at both positions"
        elif b==[0,0]:
            result="GC clamp not seen at both positions"
        else:
            result="GC clamp at one position"
        return result,seqs



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

class primer_design(primer):
    rev_p=False
    for_p=False
    def __init__(self, sequence, forward_primer, reverse_primer):
        #super().__init__(sequence)
        self.sequence=Seq(sequence)
        self.length=len(self.sequence)
        self.forward_primer = primer(forward_primer)
        self.reverse_primer = primer(reverse_primer)
        _,self.forward_primer_position =self.ــforward_primer_finder()
        _,self.reverse_primer_position =self.ــreverse_primer_finder()
        self.forward_primer_bool,_ =self.ــforward_primer_finder()
        self.reverse_primer_bool,_ =self.ــreverse_primer_finder()
        self.amplicon= self.sequence[self.forward_primer_position:self.reverse_primer_position]
        self.amplicon_tm=self.calculate_tm()

        
    def ــforward_primer_finder(self):

        d = 256
        q = 101
        M = self.forward_primer.length
        N = self.length
        i = 0
        j = 0
        p = 0 # hash value for pattern
        t = 0 # hash value for txt
        h = 1
        occurence=[]

        # The value of h would be "pow(d, M-1)%q"
        for i in range(M-1):
            h = (h*d) % q

        # Calculate the hash value of pattern and first window
        # of text
        for i in range(M):
            p = (d*p + ord(self.forward_primer.sequence[i])) % q
            t = (d*t + ord(self.sequence[i])) % q

        # Slide the pattern over text one by one
        for i in range(N-M+1):
            # Check the hash values of current window of text and
            # pattern if the hash values match then only check
            # for characters one by one
            if p == t:
                # Check for characters one by one
                for j in range(M):
                    if self.sequence[i+j] != self.forward_primer.sequence[j]:
                        break
                    else:
                        j += 1

                # if p == t and pat[0...M-1] = txt[i, i+1, ...i+M-1]
                if j == M:
                    occurence.append(i)
                    

            # Calculate hash value for next window of text: Remove
            # leading digit, add trailing digit
            if i < N-M:
                t = (d*(t-ord(self.sequence[i])*h) + ord(self.sequence[i+M])) % q

                # We might get negative values of t, converting it to
                # positive
                if t < 0:
                    t = t+q
        if not occurence:
            print("Forwrad primer doesn't align with DNA")
            return self.for_p, -1

        elif len(occurence)>1:
            print("forward primer aligns with more than one position, PCR not possible")
            return self.for_p , -1
        else:
            for_p_start_position=occurence[0]
            self.for_p=True
            return self.for_p , for_p_start_position
    
    def ــreverse_primer_finder(self):


        d = 256
        q = 101
        M = self.reverse_primer.length
        N = self.length
        i = 0
        j = 0
        p = 0 # hash value for pattern
        t = 0 # hash value for txt
        h = 1
        occurence=[]
        reverse_primer_comp=self.reverse_primer.complement_sequence
        # The value of h would be "pow(d, M-1)%q"
        for i in range(M-1):
            h = (h*d) % q

        # Calculate the hash value of pattern and first window
        # of text
        for i in range(M):
            p = (d*p + ord(reverse_primer_comp[i])) % q
            t = (d*t + ord(self.sequence[i])) % q

        # Slide the pattern over text one by one
        for i in range(N-M+1):
            # Check the hash values of current window of text and
            # pattern if the hash values match then only check
            # for characters one by one
            if p == t:
                # Check for characters one by one
                for j in range(M):
                    if self.sequence[i+j] != reverse_primer_comp[j]:
                        break
                    else:
                        j += 1

                # if p == t and pat[0...M-1] = txt[i, i+1, ...i+M-1]
                if j == M:
                    occurence.append(i)


            # Calculate hash value for next window of text: Remove
            # leading digit, add trailing digit
            if i < N-M:
                t = (d*(t-ord(self.sequence[i])*h) + ord(self.sequence[i+M])) % q

                # We might get negative values of t, converting it to
                # positive
                if t < 0:
                    t = t+q
        if not occurence:
            print("reverse primer doesn't align with DNA")
            return self.rev_p , -1
        elif len(occurence)>1:
            print("reverse primer aligns with more than one position, PCR not possible")
            return self.rev_p , -1
        else:
            rev_p_start_positin=occurence[0]
            self.rev_p=True

            return self.rev_p , rev_p_start_positin


    def calculate_tm(self):
    #tm.Tm_NN(my_dna)
        return round(TM.Tm_NN(self.sequence) , 2)
    
    def primer_analysis(self):
        analysis_table = PrettyTable(["Primer", "Tm","GC content", "Hairpin","Double Bind", "GC clamp"]) 

        if self.reverse_primer_bool and self.forward_primer_bool:
            pass
        else:
            print("primers not valids")
            return

        Ta = (0.3*(min(self.forward_primer.tm,self.reverse_primer.tm))) + 0.7*(self.amplicon_tm) - 14.9
        print(Ta)
        analysis_table.add_row(["Forward Primer", self.forward_primer.calculate_tm(), self.forward_primer.calculate_gc() , self.forward_primer.hair_pin, self.forward_primer.double_binding, self.forward_primer.GC_clamp ])
        analysis_table.add_row(["Forward Primer", self.reverse_primer.calculate_tm(), self.reverse_primer.calculate_gc() , self.reverse_primer.hair_pin, self.reverse_primer.double_binding, self.reverse_primer.GC_clamp ])
        analysis_table.get_html_string()
        print(analysis_table)



class lysis_analysis( Seq):
     
    def __init__(self,sequence):
        self.sequence=Seq(sequence)
        self.lysis_map=self.lysis()
        self.single_cut_enzymes,_=self.enzymes()
        _,self.multi_cut_enzymes=self.enzymes()


    def lysis(self):
        rb = RestrictionBatch(Restriction.AllEnzymes)
        ana=Restriction.Analysis(rb , self.sequence)
        ana.print_as('map')
        ana_pt=ana.print_that
        return ana_pt(ana.with_sites())
    

    def enzymes(self):
        rb = RestrictionBatch(Restriction.AllEnzymes)
        ana=Restriction.Analysis(rb , self.sequence)
        enzymes_cut=ana.with_sites()# gives a dict
        enzyme_single_cut={}
        enzyme_multi_cut={}

        for enzyme , positions in enzymes_cut.items():
            num_cuts=len(positions)
            if num_cuts==1:
                enzyme_single_cut[enzyme]=(num_cuts,positions)
            else:
                enzyme_multi_cut[enzyme]=(num_cuts,positions)

        return  enzyme_single_cut , enzyme_multi_cut



    

'''
The optimal annealing temperature (Ta Opt) for a given primer pair on a particular target
 can be calculated as follows: Ta Opt = 0.3 x (Tm of primer) + 0.7 x (Tm of product) – 14.9
where Tm of primer is the melting temperature of the less stable primer-template pair

Length of 18-24 bases
40-60% G/C content
Start and end with 1-2 G/C pairs
Melting temperature (Tm) of 50-60°C
Primers with melting temperatures in the range of 52-58 oC generally produce the best 
results. Primers with melting temperatures 
above 65oC have a tendency for secondary annealing. 
Primer pairs should have a Tm within 5°C of each other
Primer pairs should not have complementary regions
'''
        









# # primer_1=primer(my_dna)
# # print(primer_1.gc_content, primer_1.tm)




from Bio.Seq import Seq
import primer_design_module as pr_d






my_dna = Seq("ATTCGGGGAAAAAAATCGAAATGAATAAGCTCCCCCCCGA")


#primer analysis
pr_design=pr_d.primer_design(my_dna, "ATTCGG", "CTTA")
print(pr_design.primer_analysis())
print(pr_design.amplicon)
amplicon=pr_design.amplicon


#lysis analysis

lysis=pr_d.lysis_analysis(amplicon)
lysis.lysis_map
single_cut=lysis.single_cut_enzymes
multi_cut=lysis.multi_cut_enzymes






from Bio import SeqIO
from pyteomics import parser

def digest_peptides(filename, digest):
    list_peptides = []
    set_peptides = set()
    dict_peptides = dict()
    for record in SeqIO.parse(filename, "fasta"):
        protein_id = str(record.id)
        sequence = str(record.seq)
    # use pyteomics to generate digested peptides
        digested_peptides = parser.cleave(filename, digest)
        for pep in digested_peptides:
            list_peptides.append(pep)
            set_peptides.add(pep)
            if pep in dict_peptides.keys():
                dict_peptides[pep] += 1
            else:
                dict_peptides[pep] = 1

    print("number of peptides overall:", len(list_peptides))
    print("number of peptide sequences", len(set_peptides))
    print("number of peptide sequences in dict", len(dict_peptides.keys()))

    counter = 0
    for pep in dict_peptides.keys():
        if dict_peptides[pep] > 1:
            counter += 1

    print(counter)

    #defines value for new key(retention time) as a set
times_dict = ddict(set)
for n, pep in enumerate(peptides):
    norm_rt = normalized_times[n]
    int_rt = int(round(norm_rt * 60, 0))
    for t in range(int_rt - 15, int_rt + 15): #change to 
        times_dict[t].add(pep)
print("times dict done")
#extracts values from predicted_tr and stores them in list times
#extracts peptide seq from column seq and stores them in list peptides
#ddict store peptiedes corresponding to each retention time

#plotting bar graph
x_list = []
y_list = []
for t, pep_set in times_dict.items():
    x_list.append(t)
    y_list.append(len(pep_set))

fig2 = go.Figure()
fig2.add_trace(go.Bar(x=x_list, y=y_list, marker=dict(color=[0, 255, 255])))
fig2.update_xaxes(title_text="Retention Time (s)")
fig2.update_yaxes(title_text="Number of co-eluting peptides")
fig2.show()

print("plotting done")

import streamlit as st
import pydna, requests, sys, itertools
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Restriction import *
from pydna.dseqrecord import Dseqrecord
from pydna.seqrecord import SeqRecord
from pydna import tm
from pydna.design import primer_design
from pydna.gel import gel
from pydna.ladders import *
from pydna.ladders import HI_LO_DNA_MARKER #need
from bs4 import BeautifulSoup


# check all seq functions, work on pro_seq (ncbi_pro_seq integration)

# expand or import get_url
# import seq search functions?

# !sequence as input need, primer also has sequence
# primer site (need), primer sequence
# restriction site, restriction sequence

# Common functions 

def get_url(url, **kwargs): # need and integration possibility
  response = requests.get(url, **kwargs);

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response

def try_except(word,data):
  try:
    i = data.index(word)
  except ValueError:
    i = -1
  return i

def partial(n,i):
  d1 = ""
  d2 = ""
  if data[i-n].find("-") != -1:
    ss_text = data[i-n].split('-')
    n += 1
    for ss in ss_text:
      if ss.isdigit() == True:
        if d1 == "":
          d1 = ss
        elif d2 == "":
          d2 = ss
  return n,d1,d2

# Sequence functions 

#(working)
def gen_seq(hum_i,seq_i): 
    def gene_term_quiry(tax,gene): #import to gen_seq
        a_tax = [["Homo sapiens","human"], ["Oryctolagus cuniculus","rabbit"],["Mus musculus", "mouse"]]
        gene_term = ""
        for i in a_tax:
            if i[0] == tax or i[1] == tax:
                tax_term = '(%s [Organism] OR %s [All Fields]) AND ' % (i[0],i[1])
                tax_sci = i[0]
                gene_term = tax_term + gene + "[All Fields] AND gene[All Fields]"# complete[All Fields]" # AND cds[All Fields]"
        return gene_term,tax_sci 
    # work on taxons
    n,d1,d2 = partial(2,seq_i) # right place?
    if hum_i == -1:
      opt = "acc"; gene = ""; gene_term = data[0]; tax_sci = ""; 
    else:
      opt = "seq"; gene = data[seq_i-n]; gene_term, tax_sci = gene_term_quiry(tax,gene)
    ncbi_gen_seq(gene_term,tax_sci,d1,d2,gene,opt,seq_i,n) #change the order, check if all is needed, what can be global?


#(working)
def pla_seq(hum_i,seq_i): # need to split into pla_seq and addgene_pla_seq
    n,d1,d2 = partial(2,seq_i)
    plasmid = data[0]
    a_id = []
    a_pla = []
    if plasmid.isdigit() == False:
        path_template = "http://addgene.org/search/catalog/plasmids/?page_number=1&page_size=10&q={}"
        path = path_template.format(plasmid)
        page = requests.get(path)
        parser = BeautifulSoup(page.text, 'html.parser')
        text = str(parser)
        #st.success(text)
        n = 0
        for line in text.split('\n'):
            #st.success(line)
            if '<div class="col-xs-10">#' in line:
               lin = line.strip()
               ## example: <div class ="col-xs-10" >  # 107251</div>
               id = lin.split('#')[1].split('</div>')[0]
               a_id.append(id)
            if line == '<h3 class="search-result-title">':
                n += 1
            if n > 0:
                n += 1
            if n == 4:
                n = 0
                start = line.find(">")
                finish = line[start:].find("<")
                pla = line[start+1:start+finish]
                a_pla.append(pla)
        path_template = 'http://www.addgene.org/{}/sequences/'
        # find name 
        for p,i in zip(a_pla,a_id):
            path = path_template.format(i)
            page = requests.get(path)
            parser = BeautifulSoup(page.text, 'html.parser')
            text = str(parser)
            list_of_attributes = {"class": "copy-from form-control"}
            tags = parser.findAll('textarea', attrs=list_of_attributes)
            tag = tags[0].text
            lines = tag.split('\n')
            #ref = lines[0].strip()
            lines = lines[1:]
            dna = ('').join(lines).strip()
            with tab1:
                with st.expander("%s (%s)" % (p,i)):
                    html_string = '<a href=http://www.addgene.org/%s/sequences/><img alt="%s" src="https://static.addgene.org/addgene-core/4df3f5733a/images/common/svg/logo-addgene.svg" width="100" ></a>' % (i,p)
                    st.markdown(html_string, unsafe_allow_html=True)
                    if d1 != "" and d2 != "":
                        dna = dna[int(d1):int(d2)+1]
                    st.code(dna.upper())
                    st.code(i)
    elif plasmid.isdigit() == True:
        path_template = 'http://www.addgene.org/{}/sequences/'
        path = path_template.format(plasmid)
        page = requests.get(path)
        parser = BeautifulSoup(page.text, 'html.parser')
        text = str(parser)
        for line in text.split('\n'):
            #st.success(line)
            if "<title>Addgene:" in line:
                lin = line.strip()
                li = lin.split(' ')
                name = li[1]
        list_of_attributes = {"class": "copy-from form-control"}
        tags = parser.findAll('textarea', attrs=list_of_attributes)
        tag = tags[0].text
        lines = tag.split('\n')
        #ref = lines[0].strip()
        lines = lines[1:]
        dna = ('').join(lines).strip()
        with tab1:
            with st.expander("%s (%s)" % (name,plasmid)):
                html_string = '<a href=http://www.addgene.org/%s/sequences/><img alt="%s" src="https://static.addgene.org/addgene-core/4df3f5733a/images/common/svg/logo-addgene.svg" width="100" ></a>' % (plasmid,name)
                st.markdown(html_string, unsafe_allow_html=True)
                if d1 != "" and d2 != "":
                    dna = dna[int(d1):int(d2)+1]
                st.code(dna.upper())   
    return dna

#(working)
def pro_cod_seq(hum_i,seq_i):
    n,d1,d2 = partial(3,seq_i)
    if data[seq_i-2] == "protein": 
        opt = "acc_cds"; prot_term = data[0]; tax_sci = "" # prot = ?
        ncbiprot_seq(prot_term,tax_sci,d1,d2,opt,seq_i,n)
    return  
  
#(not necessary now, not working)
def pro_seq(hum_i,seq_i): #import ncbipro_seq
    def prot_term_query(tax,prot): 
        st.success("prot_term_query")
        a_tax = [["Homo sapiens","human"], ["Oryctolagus cuniculus","rabbit"],["Mus","Mus musculus", "mouse"]]
        prot_term = ""
        for i in a_tax:
            if i[0] == tax or i[1] == tax:
                #tax_term = '(%s [Organism] OR %s [All Fields]) AND ' % (i[0],i[1])
                tax_term = '(%s [Organism] OR %s [All Fields]) AND ' % (i[0],i[1])
                tax_sci = i[0]
                prot_term = tax_term + prot + "[Protein Name]"
        return prot_term,tax_sci  
    def prot_name(i,n): # import to pro_seq
        # prot_name between tax and i_start
        st.success("prot_name")
        i_start = i
        # delete human from name (check human appearance in names, if yes check if prot starts from upper case)
        r_data = data[::-1]
        r_start = r_data.index(i)
        r_current = r_start
        r_finish = None
        #st.success(r_data[r_start+1:])
        for i in r_data[r_start+1:]:
            if i[0].isupper() == True:
              r_finish = r_current
              #st.success(i)
            elif i[0].isupper() == False and r_finish != None:
              #st.error(i)
              break
            r_current += 1
        #st.success("r_i and r_x")
        #st.success(r_i)
        #st.success(r_x)
        if i_start == "sequence" or i_start == "primers":
            r_protein = r_data[r_start+1:r_start+1+r_finish+1]
        elif i_start == "protein":
            r_protein = r_data[r_start+1:r_start+r_finish+1]
        #st.success(r_protein)
        prot = r_protein[::-1]
        return prot
    return

# Sequence search functions (import?)

def unipro_seq(protein,tax_id): # 
  WEBSITE_API = "https://rest.uniprot.org"
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{protein}) AND (taxonomy_id:{tax_id})&fields=protein_name,gene_names,accession&size=1", headers={"Accept": "text/plain; format=fasta"})
  return r.text

def ncbiprot_seq(prot_term,tax_sci,d1,d2,opt,seq_i,n): # need to split to ncbiprot_seq and pri#optimise and apply all taxons
  Entrez.tool = 'Essequery'
  Entrez.email = ''
  h_search = Entrez.esearch(db="protein", term=prot_term, retmax=20,sort = "relevance")
  records = Entrez.read(h_search)
  h_search.close()
  identifiers = list(records['IdList'])
  num = 0
  if len(identifiers) == 0:
      st.warning("NCBI first 20 results do not match")
      prot_plus = "+".join(prot)
      link = "https://www.ncbi.nlm.nih.gov/protein/?term=%s+%s" % (tax,prot_plus)
      html_string = '<a href="%s"><img alt="" src="https://serratus.io/ncbi.png" width="100" ></a>' % link
  for i in identifiers:
        percent = round(num/len(identifiers)*100)
        with st.spinner('Filtering NCBI protein results (%s %%)' % percent):
            num += 1
            h_fetch = Entrez.efetch(db = 'protein', id =i, rettype = 'gb')
            recs = list(SeqIO.parse(h_fetch,'gb'))
            v = recs[0].description.lower()
            s = recs[0].seq
            source = "Genbank Accession no %s" % recs[0].id
            if opt == "seq":
                if (v.find(tax_sci.lower()) != -1 or v.find(tax.lower()) != -1) and v.find(protein.lower()) != -1:
                    if d1 != "" and d2 != "":
                         s = s[int(d1)-1:int(d2)]
                    with tab1:
                          with st.expander("%s (%s)" % (v,recs[0].id)):
                              html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                              st.markdown(html_string, unsafe_allow_html=True)
                              st.code(s)
                              st.code(recs[0].id)
                    with tab2:
                          with st.expander("%s (%s)" % (v,recs[0].id)):
                              html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                              st.markdown(html_string, unsafe_allow_html=True)
                              st.code(s)
                              st.code(recs[0].id)
                    with tab3:
                          exp_text = "%s %s sequence (%s)" % (tax_sci,protein,source)
                          with st.expander(exp_text):
                              st.code(s)
                              st.code(recs[0].id)
            elif opt == "cds":
                st.success("cds")
            elif opt == "acc":
                if d1 != "" and d2 != "":
                         s = s[int(d1)-1:int(d2)]
                with tab1:
                      with st.expander("%s (%s)" % (v,recs[0].id)):
                          html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                          st.markdown(html_string, unsafe_allow_html=True)
                          st.code(s)
                with tab2:
                      with st.expander("%s (%s)" % (v,recs[0].id)):
                          html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                          st.markdown(html_string, unsafe_allow_html=True)
                          st.code(s)
                with tab3:
                      exp_text = "%s %s sequence (%s)" % (tax_sci,protein,source)
                      with st.expander(exp_text):
                          st.code(s)
            elif opt == "acc_cds":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #rename to gene_term (combine all below?)
                gene_acc_split = gene_acc.split() #rename to gene_term
                gene_term = gene_acc_split[-1]
                #n,d1,d2 = partial(2)
                #if d1 != -1:
                    #acc = data[seq_i-n] 
                #else:
                    #acc = data[seq_i-n-3] 
                     #s = s[((int(d1)-1)*3):(int(d2)*3)]
                tax_sci = ""; gene = ""; opt = "prot_cds";
                ncbi_gen_seq(gene_term,tax_sci,d1,d2,gene,"prot_cds",seq_i,n)
                #ncbi_gen_seq(gene_term,tax_sci,d1,d2,gene,opt,seq_i,n):
            elif opt == "primer":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = gene_seq(gene_acc,"",d1,d2,"primer")
                t=Dseqrecord(dna)
                ampl = primer_design(t, target_tm=55.0)
                if d1 != "" and d2 != "":
                    d1_d2 = d1 + "-" + d2
                else:
                    d1_d2 = ""
                with tab1:
                    st.code(ampl.forward_primer.seq)
                    st.code(ampl.reverse_primer.seq)
                with tab2:
                    st.success("%s %s %s sequence " % (tax,protein,d1_d2))
                    st.code(s)
                    #coding sequence
                    st.success("forward primer of %s %s %s" % (tax,protein,d1_d2)) # add temperature
                    st.code(ampl.forward_primer.seq)
                    st.success("reverse primer of %s %s %s" % (tax,protein,d1_d2)) # add temperature
                    st.code(ampl.reverse_primer.seq)
                with tab3:
                    st.code(ampl.forward_primer.seq)
                    st.code(ampl.reverse_primer.seq)
                    #st.write("%s %s coding gene is %s with coding sequence %s (UniProt Consortium, 2021; NCBI Resource Coordinators ,2016)." % (tax.capitalize(),recs[0].description,gene,s))
                    #st.write("%s %s (NCBI Resource Coordinators, 2016) primers are %s and %s." 
                             #% (tax.capitalize(),gene, ampl.forward_primer.seq, ampl.reverse_primer.seq))
            elif opt == "res_primer":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = gene_seq(gene_acc,"",d1,d2,"primer")
                t=Dseqrecord(dna)
                ampl = primer_design(t, target_tm=55.0)
                and_i = try_except("and",data)
                if and_i != -1:
                    e1 = data[and_i+1]
                    e2 = data[and_i-1]
                else:
                    e1 = data[pri_i+2]
                    e2 = e1
                re_s = eval(e1+ ".site")
                fe_s = eval(e2+".site")
                rp_e = re_s + ampl.reverse_primer.seq
                fp_e = fe_s + ampl.forward_primer.seq
                st.success("forward primer with %s site" % e2)
                st.code(fp_e)
                st.success("reverse primer with %s site" % e1)
                st.code(rp_e)
            elif opt == "over_primer":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = gene_seq(gene_acc,"",d1,d2,"primer")
                t=Dseqrecord(dna)
                ampl = primer_design(t, target_tm=55.0)
                and_i = try_except("and",data)
                if and_i != -1:
                    e1 = data[and_i+1]
                    e2 = data[and_i-1]
                else:
                    e1 = data[pri_i+2]
                    e2 = e1
                re_s = eval(e1+ ".site")
                fe_s = eval(e2+".site")
                DNA = dna.upper()
                st.success(DNA)
                Reverse_DNA = DNA[::-1]
                overhang_size = 4 # take from text and use
                overhang_parts = ["A","T", "G", "C"]
                perm = itertools.product(overhang_parts, repeat=4)
                num_seq = []
                
                for i in list(perm):
                    j = "".join(i)
                    n = DNA.count(j)
                    n_r = Reverse_DNA.count(j)
                    a = [n+n_r,j]
                    num_seq.append(a)
                    #print(j,":",n+n_r,"(",n,"+",n_r,")")
                
                num_seq.sort(key=lambda x: int(x[0]))
                st.success(num_seq)
                nucl_n = 0
                st.success(num_seq[0])
                rp_e = re_s + ampl.reverse_primer.seq
                fp_e = fe_s + ampl.forward_primer.seq
                # 4 different and 4th different from first of rp_e and fp_e
                for i in num_seq: 
                    st.success(i[1])
                    for nucl in overhang_parts:
                        nucl_c = i[1].count(nucl)
                        if nucl_c > nucl_n:
                            nucl_n = nucl_c
                            st.success(nucl_n)
                    if nucl_n == 1:
                        overhang = i[1]
                        orp_e = overhang + re_s + ampl.reverse_primer.seq
                        ofp_e = overhang + fe_s + ampl.forward_primer.seq
                        st.success("forward primer with %s site" % e2)
                        st.code(fp_e)
                        st.success("reverse primer with %s site" % e1)
                        st.code(rp_e)
                        break
                    nucl_n = 0
            elif opt == "res_all" or opt =="res" or opt == "res_seq" or opt =="res_gel" or opt == "res_clo":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = gene_seq(gene_acc,"",d1,d2,"primer")
                t= Seq(dna)
                #check len(t)
                a_site = []
                a_c = []
                #st.success(opt)
                if opt == "res" or opt == "res_seq" or opt == "res_gel":
                    wit_i = try_except("with",data)
                    com_i = try_except(",",data)
                    if com_i == -1:
                        rb = RestrictionBatch([data[wit_i+1]])
                    else:
                        rb = RestrictionBatch(' '.join(data[wit_i+1:]).split(","))
                    test = rb.search(t)
                    for a,b in test.items():
                        if opt == "res":
                            st.success("%s %s" %(a,b))
                        if opt == "res_seq" or "res_gel":
                            a_site.append(b)
                if opt == "res_clo":
                    a_clo = []
                    rb = AllEnzymes.search(t)
                    for a,b in rb.items():
                        if len(b) != 0:
                            a_clo.append(str(a))
                    return a_clo
                if opt == "res_all":
                    rb = AllEnzymes.search(t)
                    for a,b in rb.items():
                        if len(b) != 0:
                            #st.success("%s %s" %(a,b))
                            c = str(a)+str(b)
                            a_clo.append(a)
                            a_c.append(c)
                    for s in sorted(a_c):
                        st.success(s)
                if len(a_site) != 0 and opt != "res_all":
                    # if circular (add)
                    # ...
                    # if linear 
                    a_site = list(itertools.chain.from_iterable(a_site))
                    a_site.sort()
                    site_0 = 1
                    a_at = []
                    a_att = []
                    for site_1 in a_site:
                        at = t[site_0-1:site_1-1]
                        a_at.append(at)
                        site_0 = site_1
                        att = Dseqrecord(at)
                        a_att.append(att)
                        if opt != "res_gel":
                            st.success(at)
                    at = t[site_1-1:]
                    a_at.append(at)
                    att = Dseqrecord(at)
                    a_att.append(att)
                    if opt != "res_gel":
                        st.success(at)
                        st.success(len(att))
                    else:
                        a_max = len(max(a_att, key=len))
                        a_min = len(min(a_att, key=len))
                        and_i = try_except('and',data)
                        ruler = ""
                        if and_i == -1:
                            if a_max <= 10000:
                                if a_min >= 250:
                                    ruler = "GeneRuler_1kb"
                                    st.success(ruler)
                                elif a_min < 250:
                                    ruler = "GeneRuler_1kb_plus"
                            #elif a_max > 10000:
                                #ruler = "HI_LO_DNA_MARKER" doesn't work
                        else:
                            ruler_input = data[and_i+1:]
                            ruler = "_".join(ruler_input)
                        try:
                            exec("st.image(gel([%s,a_att]))" % ruler)
                        except NameError:
                            st.info("Please use one of these ladders: GeneRuler 1kb, GeneRuler 1kb plus or PennStateLadder")
                     
def prot_gene_seq(protein,tax_id,tax):
  WEBSITE_API = "https://rest.uniprot.org"
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{protein}) AND (taxonomy_id:{tax_id})&fields=gene_names", headers={"Accept": "text/plain; format=tsv"})
  entities = r.text
  a_entity = entities.split('\n')
  a_gene = a_entity[1].split(' ')
  gene = a_gene[0] #optimize
  return gene

def ncbi_gen_seq(gene_term,tax_sci,d1,d2,gene,opt,seq_i,n): #check if all is needed, what can be global?
    Entrez.tool = 'Essequery'
    Entrez.email = ''
    h_search = Entrez.esearch(db="nucleotide", term=gene_term, retmax=20, sort = "relevance", api_key = "7c824eac4588c5996739d4a6c136c3f5d808")
    records = Entrez.read(h_search)
    h_search.close()
    identifiers = list(records['IdList'])
    num = 0
    if len(identifiers) == 0:
        st.warning("NCBI first 20 results do not match")
    for i in identifiers:
          percent = round(num/len(identifiers)*100)
          with st.spinner('Filtering NCBI nucleotide results (%s %%)' % percent):
              num += 1
              h_fetch = Entrez.efetch(db = 'nucleotide', id =i, rettype = 'gb', api_key = "7c824eac4588c5996739d4a6c136c3f5d808")
              recs = list(SeqIO.parse(h_fetch,'gb'))
              v = recs[0].description.lower()
              s = recs[0].seq
              source = "Genbank Accession no %s" % recs[0].id
              if opt == "seq":
                if v.find(tax_sci.lower()) != -1 and v.find(gene.lower()) != -1:
                      if d1 != "" and d2 != "":
                              s = s[int(d1):int(d2)+1]
                      with tab1:
                            with st.expander("%s (%s)" % (v,recs[0].id)):
                                html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                                st.markdown(html_string, unsafe_allow_html=True)
                                st.code(s)
                      with tab2:
                            with st.expander("%s (%s)" % (v,recs[0].id)):
                                html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                                st.markdown(html_string, unsafe_allow_html=True)
                                st.code(s)
                      with tab3:
                            exp_text = "%s %s sequence (%s)" % (tax_sci,data[seq_i-n],source)
                            with st.expander(exp_text):
                                st.code(s)     
                              # speed up by query_list[0] try to show up, more button
                              #query_list.append(recs)
              if opt == "acc":
                    if d1 != "" and d2 != "":
                              s = s[int(d1):int(d2)+1]
                    with tab1:
                          with st.expander("%s (%s)" % (v,recs[0].id)):
                              html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                              st.markdown(html_string, unsafe_allow_html=True)
                              st.code(s)
                    with tab2:
                          with st.expander("%s (%s)" % (v,recs[0].id)):
                              html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                              st.markdown(html_string, unsafe_allow_html=True)
                              st.code(s)
                          with tab3:
                              exp_text = "%s %s sequence (%s)" % (tax_sci,data[seq_i-n],source)
                              with st.expander(exp_text):
                                 st.code(s)
              if opt == "prot_cds":
                a_f = recs[0].features
                for f in a_f:
                    if f.type == "CDS":
                        s = f.location.extract(recs[0]).seq
                        if d1 != "" and d2 != "":
                          if d1 == "1":
                             s = s[:(int(d2)*3)]
                          else:
                             s = s[((int(d1)-1)*3):(int(d2)*3)]
                        with tab1:
                              with st.expander("%s (%s)" % (v,recs[0].id)):
                                  html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                                  st.markdown(html_string, unsafe_allow_html=True)
                                  st.code(s)
                        with tab2:
                              with st.expander("%s (%s)" % (v,recs[0].id)):
                                  html_string = '<a href="https://www.ncbi.nlm.nih.gov/nuccore/%s"><img alt="%s" src="https://serratus.io/ncbi.png" width="100" ></a>' % (i,recs[0].id)
                                  st.markdown(html_string, unsafe_allow_html=True)
                                  st.code(s)
                        with tab3:
                              exp_text = "%s %s coding sequence (%s)" % (tax_sci,data[seq_i-n],source)
                              with st.expander(exp_text):
                                  st.code(s)

# Primer functions (make class)

def pro_pri():
    return

### main

# title
st.title('ClonePot') 
# subheader
st.subheader('Cloning One-Pot design assistant') 
# input
task = st.text_input("Please enter the task", "") 

# button click
if(st.button('Submit')):
    global tab1,tab2,tab3
    tab1, tab2, tab3 = st.tabs(["Result-only", "Step-by-step", "Article"])
    data = task.split()
    
    # biological host (to enhance)
    hum_i = try_except('human',data)
    if hum_i != -1:
        tax = data[hum_i]
        if tax == 'human':
            tax_id = 9606
    else:
        tax = 'human'
        tax_id = 9606

    #compulsary keywords
    a_try_except = ['gene','plasmid','protein','coding','restriction','sequence','primers','cloning'] #cds instead of coding sequence?
    #optional keywords
    aa_try_except = ['for','with',',','and']
    
    #find compulsary keywords in input and generate function
    func = ""
    for ate in a_try_except:
        abbr = ate[:3]+ "_i"
        exec("%s = try_except('%s',data)" % (abbr,ate))
        exec("if %s != -1: func += ate[:3] + '_'" % abbr)
    func = func[:len(func)-1]
    
    #check function
    st.success(func)
    
    #execute function
    exec('%s(hum_i,seq_i)' % func)
 
### main end 
 
    #work on primer functions
    #look architecture difference between function "primers" and "primers with"
    #elif pri_i != -1: 
        #n,d1,d2 = partial(2,seq_i)
        #if  data[pri_i-1] == "protein":
            #protein = data[0]
            #if data[pri_i] == data[-1]:
                #ncbiprot_seq(protein,"",d1,d2,"primer") 
            #elif data[-1] == "overhang":
                # option of custom size and custom aminoacids
                #ncbiprot_seq(protein,"",d1,d2,"over_primer")
                # work on over_primer
            #elif data[pri_i+1] == "with":
                #ncbiprot_seq(protein,"",d1,d2,"res_primer")
        # add overhang    
        #elif data[pri_i-1] == "gene":
            #st.success("gene primers")
    

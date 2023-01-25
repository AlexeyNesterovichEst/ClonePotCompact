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
from pydna.amplify import pcr
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np


# check all seq functions, work on pro_seq (ncbi_pro_seq integration)

# expand or import get_url
# import seq search functions?

# !sequence as input need, primer also has sequence
# primer site (need), primer sequence
# restriction site, restriction sequence

# Common functions 

global protein

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
    ncbi_gen_seq(n,d1,d2,opt,gene,gene_term,tax_sci)

#(working)
def pla_seq(opt,seq_i):
    n,d1,d2 = partial(2,seq_i)
    if opt == "clo":
        plasmid0 = data[0].split("-")
        plasmid = plasmid0[0]
    else:
        plasmid = data[0]
    a_id = []
    a_pla = []
    #addgene starts here
    if plasmid.isdigit() == False:
        path_template = "http://addgene.org/search/catalog/plasmids/?page_number=1&page_size=10&q={}"
        path = path_template.format(plasmid)
        page = requests.get(path)
        parser = BeautifulSoup(page.text, 'html.parser')
        text = str(parser)
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
                    html_string = '<a href=http://www.addgene.org/%s/sequences/><img alt="%s" src="https://stemadvocacy.org/wp-content/uploads/2019/07/addgene-logo.png" width="100" ></a>' % (i,p)
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
        if opt != "clo":
            with tab1:
                with st.expander("%s (%s)" % (name,plasmid)):
                    html_string = '<a href=http://www.addgene.org/%s/sequences/><img alt="%s" src="https://stemadvocacy.org/wp-content/uploads/2019/07/addgene-logo.png" width="100" ></a>' % (plasmid,name)
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
        ncbiprot_seq(prot_term,tax_sci,d1,d2,opt,seq_i,n) #check seq_i and n need, change order
    return  
  
#(working)
def pro_seq(hum_i,seq_i): #import ncbipro_seq
    def prot_term_query(tax,prot): 
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
    n,d1,d2 = partial(2,seq_i)
    if hum_i == -1: #(working)
        global protein
        protein = data[0]
        ncbiprot_seq(protein,"",d1,d2,"acc")  
    else:
        #check if prot_name and prot_term_query can be combined?
        prot = prot_name("protein",3)
        prot_s = ' '.join(prot) #check if necessary
        if d1 != "" and d2 != "":
            prot = prot[:-1]
        protein = ' '.join(prot)
        prot_term,tax_sci = prot_term_query(tax,protein)
        ncbiprot_seq(prot_term,tax_sci,d1,d2,"seq")
    return

# Sequence search functions (import?)

def unipro_seq(protein,tax_id): # 
  WEBSITE_API = "https://rest.uniprot.org"
  r = get_url(f"{WEBSITE_API}/uniprotkb/search?query=(protein_name:{protein}) AND (taxonomy_id:{tax_id})&fields=protein_name,gene_names,accession&size=1", headers={"Accept": "text/plain; format=fasta"})
  return r.text

def ncbiprot_seq(prot_term,tax_sci,d1,d2,opt): # need to split to ncbiprot_seq and pri#optimise and apply all taxons
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
            elif opt == "primer": # "pcr"
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = ncbi_gen_seq(3,d1,d2,"primer","",gene_acc,"")
                t=Dseqrecord(dna)
                ampl = primer_design(t, target_tm=55.0)
                f_p = ampl.forward_primer.seq
                r_p = ampl.reverse_primer.seq
                if d1 != "" and d2 != "":
                    d1_d2 = d1 + "-" + d2
                else:
                    d1_d2 = ""
                return f_p,r_p,dna
                
            elif opt == "res_primer":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = ncbi_gen_seq(3,d1,d2,"primer","",gene_acc,"")
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
                return fp_e,rp_e,dna,e1,e2
            elif opt == "over_primer":
                handle = Entrez.efetch(db="protein", id=recs[0].id,rettype="gb")
                seq_record = SeqIO.read(handle, "genbank")
                seqAnn = seq_record.annotations
                gene_acc = seqAnn['db_source'] #test different proteins
                gene_acc_split = gene_acc.split()
                gene_acc = gene_acc_split[-1]
                dna = ncbigen_seq(gene_acc,"",d1,d2,"primer")
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
                dna = ncbi_gen_seq(3,d1,d2,"primer","",gene_acc,"")
                t= Seq(dna) # check len(t)
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
                    a_clo = []
                    b_clo = []
                    c_clo = []
                    for a,b in test.items():
                        if opt == "res":
                            a_clo.append(str(a))
                            b_clo.append(b)
                            c_clo.append(len(b))
                        if opt == "res_seq" or "res_gel":
                            a_site.append(b)
                    d = {'name': a_clo, 'site': b_clo, 'n': c_clo}
                    df = pd.DataFrame(data=d)
                    st.dataframe(data=df)
                if opt == "res_clo":
                    a_clo = []
                    rb = CommOnly.search(t)
                    for a,b in rb.items():
                        if len(b) != 0:
                            a_clo.append(str(a))
                    return a_clo
                if opt == "res_all":
                    a_clo = []
                    b_clo = []
                    c_clo = []
                    rb = CommOnly.search(t)
                    for a,b in rb.items():
                        if len(b) != 0:
                            #st.success("%s %s" %(a,b))
                            a_clo.append(str(a))
                            b_clo.append(b)
                            c_clo.append(len(b))
                    return a_clo,b_clo,c_clo
                if len(a_site) != 0 and opt != "res_all":
                    a_site = list(itertools.chain.from_iterable(a_site))
                    a_site.sort()
                    site_0 = 1
                    a_at = [] 
                    for site_1 in a_site:
                        at = t[site_0-1:site_1-1]
                        a_at.append(at)
                        site_0 = site_1
                    at = t[site_1-1:]
                    a_at.append(at)
                    #a_at[0] = at + a_at[0] #for circular
                    if opt == "res_gel":
                        a_att = []
                        for i in a_at:
                            att = Dseqrecord(i)
                            a_att.append(att)
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

def ncbi_gen_seq(n,d1,d2,opt,gene,gene_term,tax_sci): #check if all is needed, what can be global?
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
              if opt == "primer":
                a_f = recs[0].features
                for f in a_f:
                    if f.type == "CDS":
                        s = f.location.extract(recs[0]).seq
                        if d1 != "" and d2 != "":
                          if d1 == "1":
                             s = s[:(int(d2)*3)]
                          else:
                             s = s[((int(d1)-1)*3):(int(d2)*3)]
                        return s

# PCR functions

def pro_pri(hum_i,pri_i): # protein primers
    n,d1,d2 = partial(2,pri_i)
    protein = data[0]
    if data[pri_i] == data[-1]:
        f_p, r_p, dna = ncbiprot_seq(protein,"",d1,d2,"primer") 
        text1 = ""
        text2 = ""
    elif data[pri_i+1] == "with":
        f_p, r_p, dna,e1,e2 = ncbiprot_seq(protein,"",d1,d2,"res_primer") 
        text1 = "with " + e2 + " site"
        text2 = "with " + e1 + " site"
    #elif data[-1] == "overhang":
        # option of custom size and custom aminoacids
        #ncbiprot_seq(protein,"",d1,d2,"over_primer")
        # work on over_primer
    # add overhang    
    #elif data[pri_i-1] == "gene":
        #st.success("gene primers")
    #viz_pro_pri()
    st.success("template")
    st.code(dna)
    st.success("forward primer %s" % text1) 
    st.code(f_p)
    st.success("reverse primer %s" % text2) 
    st.code(r_p)
    return f_p,r_p,dna

def pro_pcr(hum_i,pcr_i): # protein pcr
    f_p,r_p,dna = pro_pri(hum_i,pcr_i)
    pcr_p = pcr(f_p, r_p, dna)
    st.success("PCR product")
    st.code(pcr_p.seq)
    return pcr_p

# Restriction functions

def pro_res(hum_i,res_i): # protein restriction (working)
    n,d1,d2 = partial(2,res_i)
    protein = data[0]
    if data[res_i] == data[-1]: 
        a_clo,b_clo,c_clo = ncbiprot_seq(protein,"",d1,d2,"res_all")
        d = {'name': a_clo, 'site': b_clo, 'n': c_clo}
        df = pd.DataFrame(data=d)
        st.dataframe(data=df)
    elif data[res_i+1] == "with": #(not working yet, check gen(e)_seq)
        ncbiprot_seq(protein,"",d1,d2,"res")
        #change to ->
        #ncbiprot_seq
        #res_pre()
    return a_clo
    
def pla_res(hum_i,res_i): # plasmid restriction (working)
    dna = pla_seq(hum_i,res_i)
    if data[res_i] == data[-1]:
        a_clo,b_clo,c_clo = res_all(dna)
        #pla_res_viz
        d = {'name': a_clo, 'site': b_clo, 'n': c_clo}
        df = pd.DataFrame(data=d)
        st.dataframe(data=df)
        #
    elif data[res_i+1] == "with":
        res_pre(dna)
    #elif single-site restriction
        #one_clo = []
        #if c_clo == 1:
            #one_clo.append(c_clo)
            
def clo_res(hum_i,res_i): #(working)
    parts = data[0].split("-")
    pla = parts[0]
    ins = parts[1:]
    if len(ins) == 1:
        ins = ins[0]
        n,d1,d2 = partial(2,res_i)
        ins_res,b_clo,c_clo = ncbiprot_seq(ins,"",d1,d2,"res_all") # insert restrictases
        pla_dna = pla_seq("clo",res_i)
        pla_res_one, pla_res_one_pos = res_one(pla_dna) # plasmid restrictases
        pla_res_one_pos = [j for i in pla_res_one_pos for j in i]
        clo_res_one = [x for x in pla_res_one if x not in ins_res]
        clo_res_one_pos = []
        for x in clo_res_one:
           pos = pla_res_one.index(x)
           clo_res_one_pos.append(pla_res_one_pos[pos])
        d = {'name': clo_res_one, 'site': clo_res_one_pos}
        df = pd.DataFrame(data=d)
        st.dataframe(data=df)
        
#instead of res_one and res_all -> res_n
def res_one(dna):
    a_clo = []
    b_clo = []
    t= Seq(dna)
    rb = CommOnly.search(t)
    for a,b in rb.items():
        if len(b) == 1:
            a_clo.append(str(a))
            b_clo.append(b)
    return a_clo,b_clo


def res_all(dna): # restrictases all   
    a_clo = []
    b_clo = []
    c_clo = []
    t= Seq(dna)
    rb = CommOnly.search(t)
    for a,b in rb.items():
        if len(b) != 0:
            a_clo.append(str(a))
            b_clo.append(b)
            c_clo.append(len(b))
    return a_clo,b_clo,c_clo
    
def res_pre(dna): #restrictases predefined
    wit_i = try_except("with",data)
    com_i = try_except(",",data)
    if com_i == -1:
        rb = RestrictionBatch([data[wit_i+1]])
    else:
        rb = RestrictionBatch(' '.join(data[wit_i+1:]).split(","))
    t= Seq(dna)
    test = rb.search(t)
    a_clo = []
    b_clo = []
    c_clo = []
    for a,b in test.items():
        a_clo.append(str(a))
        b_clo.append(b)
        c_clo.append(len(b))
    d = {'name': a_clo, 'site': b_clo, 'n': c_clo}
    df = pd.DataFrame(data=d)
    st.dataframe(data=df)

def mcs_res(hum_i,mcs_i): # (not working)
    dna = pla_seq(hum_i,mcs_i)
    st.success(dna)
    t= Seq(dna)
    y_clo = []
    z_clo = []
    rb = CommOnly.search(t)
    for y,z in rb.items():
        if len(z) == 1:
            for i in z:
                z_clo.append(i)
    x_site = sorted(z_clo) # BioEng: is mcs bigger?
    x_site_dif = []

    n = 0
    while True:
        try:
            x_site_dif.append(x_site[n+1]-x_site[n])
        except:
            break
        n += 1
    
    ##first high-distribution group
    mcs_0 = []
    mcs_1 = []
    for x,xd in zip(x_site,x_site_dif): #41,40
        if xd <= 10: # think how to make it changeable
            mcs_0.append(x)
            if len(mcs_0) > len(mcs_1):
                mcs_1.clear()
                mcs_1.extend(mcs_0)
        else:
            mcs_0.clear()
    
    #find first mcs index        
    mcs_start = x_site.index(mcs_1[0])
    
    #find mcs_positions
    mcs_pos = x_site[mcs_start:mcs_start+len(mcs_1)+1]
    st.success(mcs_pos)
    # make a table of restrictases
    
    ##check if there's a second high-distribution group
    
    ##if yes, 
    
def mcs_clo_res():
    #mcs_res
    #clo_res
    return

def clo():
    #mcs_clo_res()
    #pcr()
    #ligation in pcr.py 
    return
    
#option = st.selectbox('Please choose the topic', ['Sequence', 'PCR', 'Restriction', 'Cloning', 'Gel electrophoresis'])

### main

# title
st.title('ClonePot') 
# subheader
st.subheader('Cloning one-pot design') 
# input
task = st.text_input("Please enter the task", "")  # change to textbox?

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
    a_try_except = ['gene','plasmid','protein','coding','cloning','restriction','sequence','primers','pcr','mcs','test'] #cds instead of coding sequence?
    
    #find compulsary keywords in input and generate function
    func = ""
    for ate in a_try_except:
        abbr = ate[:3]+ "_i"
        exec("%s = try_except('%s',data)" % (abbr,ate))
        exec("if %s != -1: func += ate[:3] + '_'" % abbr)
    func = func[:len(func)-1]
    
    #check function
    st.success(func)
    
    #execute function (check x_i usage and appearance in "for" loop of a_try_except)
    if seq_i != -1:
        exec('%s(hum_i,seq_i)' % func)
    elif pri_i != -1:
        exec('%s(hum_i,pri_i)' % func)
    elif res_i != -1:
        exec('%s(hum_i,res_i)' % func)
    elif mcs_i != -1:
        exec('%s(hum_i,mcs_i)' % func)
    elif pcr_i != -1:
        exec('%s(hum_i,pcr_i)' % func)
    elif clo_i != -1:
        exec('%s(hum_i,clo_i)' % func)
    while tes_i != -1:
        res_1 = st.selectbox('Please choose the first restrictase', ['One', 'Two', 'Three', 'Four', 'Five'], key = "selectbo")
        st.success(res_1)
            #function with list of restrictases return
            #if option:
                #st.write("You entered: ", option)
            #if option0 != option1:
                #st.success(option)
            #option2 = st.selectbox('Please choose the second restrictase', ['1', '2', '3', '4', '5'])

### main end 
    
# single-sequence functions (instead of expander, container)

# update for sequence functions:
    # data scratching part - find main parameters and standardize names
    # other parts - find main parameters and standardize names

#update for restriction functions:
     #change to ->
     #ncbiprot_seq
     #res_all()/res_pre()
     
    #option = st.selectbox('Please choose restrictases', mcs_res)
    # multiselect with compatible restrictases update/two select boxes
    
#update primer vizualization:
    # biotite alignment-based primer alignment figure (https://github.com/biotite-dev/biotite/blob/master/src/biotite/sequence/graphics/alignment.py)
    # streamlit pyplot (https://docs.streamlit.io/library/api-reference/charts/st.pyplot)
    
#move all vizualizations to main functions -> define all vizualizations (example def pro_pri)

# update for sequence functions:
    # upload custom sequences

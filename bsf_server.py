#!/anaconda/bin/python

"""
Very simple HTTP server in python.
Usage::
    ./bsf_server.py [<port>]
Send a GET request::
    curl http://localhost
Send a HEAD request::
    curl -I http://localhost
Send a POST request::
    curl -d "up=PIK3C3&dn=DFFB" http://localhost
"""

import bsf
import numpy as np
import pandas as pd
import gzip
from http.server import BaseHTTPRequestHandler, HTTPServer
import time
import urllib.parse
from numba import jit
import os
import json
import scipy.stats as stats

####################################################################################################
# DATA FILE PATHS
####################################################################################################
working_dir = './'
# gene ids (pr_gene_id) of rows
path_rids = working_dir + 'lincs/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.row_meta.txt'
# signature ids (sig_id) of columns
path_cids = working_dir + 'lincs/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.col_meta.txt'
# gene info
path_gene_info = working_dir + 'lincs/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt'
# signature info
path_sig_info = working_dir + 'lincs/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt'
path_sig_metrics = working_dir + 'lincs/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt'
# up/dn
path_up = working_dir + 'lincs/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06_up_2.0.bin.gz'
path_dn = working_dir + 'lincs/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06_dn_2.0.bin.gz'

# msigdb
path_msigdb_bin = working_dir + 'msigdb/msigdb.v6.1.symbols.bin.gz'
path_msigdb_genes = working_dir + 'msigdb/msigdb.v6.1.symbols.row_meta.txt'
path_msigdb_genesets = working_dir + 'msigdb/msigdb.v6.1.symbols.col_meta.txt'

# kegg
path_kegg_bin = working_dir + 'kegg/kegg.pathway.10-20-2017.bin.gz'
path_kegg_genes = working_dir + 'kegg/kegg.pathway.10-20-2017.row_meta.txt'
path_kegg_genesets = working_dir + 'kegg/kegg.pathway.10-20-2017.col_meta.txt'

####################################################################################################
# READ METADATA FILES
####################################################################################################
rids = pd.read_csv(path_rids, sep='\t')
cids = pd.read_csv(path_cids, sep='\t')
gene_info = pd.read_csv(path_gene_info, sep='\t')
sig_info = pd.read_csv(path_sig_info, sep='\t', index_col=0)
sig_metrics = pd.read_csv(path_sig_metrics, sep='\t', index_col=0)

lincs_genes = rids.merge(gene_info, left_on="rid", right_on="pr_gene_id", how='left')
lincs_sigs = cids.merge(sig_info, left_on="cid", right_on="sig_id", how='left')

msigdb_genes = pd.read_csv(path_msigdb_genes, sep='\t', index_col=0)
msigdb_genesets = pd.read_csv(path_msigdb_genesets, sep='\t', index_col=0)

kegg_genes = pd.read_csv(path_kegg_genes, sep='\t', index_col=0)
kegg_genesets = pd.read_csv(path_kegg_genesets, sep='\t', index_col=0)
####################################################################################################

@jit
def bool2int(x):
    y = 0
    for i, j in enumerate(x):
        if j:
            y += 1 << i
    return y


def bool2ints(fn, bits, _size):
    newsize = _size * 64
    r = np.zeros(_size, dtype=np.uint64)
    start = time.time()
    for j in range(_size):
        _end = min((j + 1) * 64, newsize)
        r[j] = fn(bits[j * 64:_end][::-1])  # reverse)

    done = time.time()
    elapsed = done - start
    print(fn, elapsed)
    return r

# read matrix (from gz)
def readLINCSGzip(filename):
    input_file = gzip.open(filename, 'rb')
    try:
        nrows = int.from_bytes(input_file.read(4), byteorder='little')
        ncols = int.from_bytes(input_file.read(4), byteorder='little')
        dt = np.dtype(np.uint64)
        dt = dt.newbyteorder('L')
        print('length of signatures:', ncols, ', number of signatures:', nrows)
        t = np.frombuffer(input_file.read(), dtype=np.uint64).reshape((nrows, ncols))
    finally:
        input_file.close()
    return ncols, nrows, t

def bin2genes(binarr, gene_symbols):
    """convert an array of 64bit binarized vector into bits and extract indexs of genes
        binarr: uint64 1d numpy array
    """
    genes = []
    for i, bit64 in enumerate(binarr):
        bitarr = np.binary_repr(bit64, 64)
        for j, bit in enumerate(bitarr):
            if bit == '1':
                genes.append(gene_symbols[i*64+j])
    return genes

def bin64_2_bin(binarrs, gene_info, score_for_ones):
    """convert an array of 64bit binarized vector into binarized vector
        binarrs: uint64 2d numpy array of (num_sigs, gene_length/64)
    Return
        _binarrs 2d numpy array of (num_genes, num_sigs)
    """
    num_sigs, leng = binarrs.shape
    num_genes, ncols = gene_info.shape
    _binarrs = np.zeros((num_genes, num_sigs), dtype=np.int)
    for k in range(num_sigs):
        binarr = binarrs[k, :]
        for i, bit64 in enumerate(binarr):
            bitarr = np.binary_repr(bit64, 64)
            for j, bit in enumerate(bitarr):
                if bit == '1':
                    _binarrs[i*64+j, k] = score_for_ones
    return _binarrs

def fisher_exact(row):
    a = row['overlap']
    b = row['DE_genes'] - a
    c = row['ngenes'] - a
    d = row['backgrounds'] - row['DE_genes'] - c
    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    return pvalue

# Load all the data files
print('1. Loading all the data files...')
print('   Loading up/down bits...')
length, nsignatures, lincsup = readLINCSGzip(path_up)
length, nsignatures, lincsdn = readLINCSGzip(path_dn)

n_msigdb_genes, n_msigdb_genesets, msigdb_bin = readLINCSGzip(path_msigdb_bin)
n_kegg_genes, n_kegg_genesets, kegg_bin = readLINCSGzip(path_kegg_bin)
print('Finish loading...')



###################################################
# TESTING
###################################################
num_test = 1000

upbits = lincsup[np.arange(num_test),:]

start = time.time()
up_binmat = bin64_2_bin(upbits, lincs_genes, 1)
done = time.time()

temp_sum = np.sum(up_binmat, axis=1)

print(up_binmat.shape)
print(temp_sum.shape)
print("non-zeros:{0}".format(np.sum(temp_sum>0)))
print('elapsed', (done - start))

start = time.time()
for i in range(num_test):
    up_genes = bin2genes(upbits[i,:], lincs_genes.pr_gene_symbol)
    # print(len(up_genes))
    assert len(up_genes) == np.sum(up_binmat[:,i])
done = time.time()
print('elapsed', (done - start))
###################################################

class BSFServer(BaseHTTPRequestHandler):
    # def get_gene_matrix(self, sig_ids, geneset_info, gene_info, bin64matrix):
    #     if len(sig_ids) > 10000: 
    #         print ("[Warning] Too large to send. Please use less than 10k")
    #         return None

    #     sig_info = geneset_info[geneset_info.sig_id.isin(sig_ids)]
    #     sig_index = sig_info.index
        
    #     upbits = bin64matrix[sig_index,:]
    #     dnbits = bin64matrix[sig_index,:]

    #     _binmat = bin64_2_bin(upbits, gene_info, 2)
    #     _binmat += bin64_2_bin(dnbits, gene_info, 1)

    #     return _binmat

    def get_lincs_gene_matrix(self, query):
        if 'sids' in query: q_sids = query['sids']
        elif 'sids[]' in query: q_sids = query['sids[]']
        else: return bytes("Not Found", "utf-8")

        sig_info = lincs_sigs[lincs_sigs.sig_id.isin(q_sids)]
        sig_index = sig_info.index

        upbits = lincsup[sig_index,:]
        dnbits = lincsdn[sig_index,:]
        _binmat = bin64_2_bin(upbits, lincs_genes, 2)
        _binmat += bin64_2_bin(dnbits, lincs_genes, 1)

        df = pd.DataFrame(_binmat, columns=list(sig_info.sig_id), index=list(lincs_genes.pr_gene_symbol)) 

        return bytes("lincs_gene_matrix("+df.to_json(orient='split')+")", "utf-8")


    def find_genesets(self, q_genes, genes_df, colname, bin_genesets, fout="test"):
        """find relevant genesets by filtering via BSF
        """
        n_genesets, n_genes = bin_genesets.shape
        scores = np.zeros((n_genesets, 1), dtype=np.uint32)
        if len(q_genes) > 0:
            gidx = np.unique(genes_df[genes_df[colname].isin(q_genes)].index)
            qbin = np.zeros(n_genes*64, dtype=bool)
            qbin[gidx] = True
            q = bool2ints(bool2int, qbin, n_genes).reshape((-1, n_genes))
            # run bsf
            bsf.analysis_with_query(bin_genesets, q, fout, '')
            # read bsf
            with open("bin_"+fout+".bin", "rb") as f:
                scores = np.frombuffer(f.read(), dtype=np.uint32).reshape((n_genesets, 1))
        return scores

    def find_lincs(self, query, genes_df, colname, sigs_df):
        # set the limit
        if 'limit' in query:
            print(query['limit'])
            limit = int(query['limit'][0])
        else:
            limit = 100

        q_upgenes = q_dngenes = []
        if 'up[]' in query: q_upgenes = query['up[]']
        if 'dn[]' in query: q_dngenes = query['dn[]']
        
        qupup = self.find_genesets(q_upgenes, genes_df, colname, lincsup, fout="qupup")
        qdndn = self.find_genesets(q_dngenes, genes_df, colname, lincsdn, fout="qupup")
        qupdn = self.find_genesets(q_dngenes, genes_df, colname, lincsup, fout="qupup")
        qdnup = self.find_genesets(q_upgenes, genes_df, colname, lincsdn, fout="qupup")

        sigs_df['upup'] = qupup.flatten()
        sigs_df['dndn'] = qdndn.flatten()
        sigs_df['updn'] = qupdn.flatten()
        sigs_df['dnup'] = qdnup.flatten()

        pos_scores = (qupup + qdndn).flatten()
        neg_scores = (qdnup + qupdn).flatten()

        # pos_jar = qupup/np.sum(qbinup) + qdndn/np.sum(qbindn)
        # neg_jar = qdnup/np.sum(qbinup) + qupdn/np.sum(qbindn)

        pos_scores_idx = np.argsort(pos_scores)[-limit:][::-1]
        neg_scores_idx = np.argsort(neg_scores)[-limit:][::-1]

        # pos_jar_idx = np.argsort(pos_jar, axis=0)
        # neg_jar_idx = np.argsort(neg_jar, axis=0)

        print(pos_scores_idx[:5])
        print(pos_scores[pos_scores_idx[:5]])

        print(neg_scores_idx[:5])
        print(neg_scores[neg_scores_idx[:5]])

        pos = sigs_df.iloc[pos_scores_idx]
        neg = sigs_df.iloc[neg_scores_idx]

        # to compute p_values
        # pos['overlap'] = pos['upup'] + pos['dndn']
        # neg['overlap'] = neg['updn'] + neg['dnup']
        # pos['DE_genes'] = len(q_upgenes) + len(q_dngenes)
        # pos['backgrounds'] = genes_df.shape[0]
        # pos['pvalue'] = top.apply(fisher_exact, axis=1)
        # neg['DE_genes'] = len(q_upgenes) + len(q_dngenes)
        # neg['backgrounds'] = genes_df.shape[0]
        # neg['pvalue'] = top.apply(fisher_exact, axis=1)

        upgene_symbols = np.unique(genes_df[genes_df[colname].isin(q_upgenes)][colname])
        dngene_symbols = np.unique(genes_df[genes_df[colname].isin(q_dngenes)][colname])
        # qgenes = dict()
        qgenes = {
            'upgenes': upgene_symbols,
            'dngenes': dngene_symbols
        }
        return bytes("lincs_table("+
            json.dumps({
                'pos': pos.to_json(orient='records'),
                'neg': neg.to_json(orient='records')})+")", "utf-8")
                # 'genes':genes})+")", "utf-8")
                # 'upgenes': upgene_symbols,
                # 'dngenes': dngene_symbols})+")", "utf-8")
    
    def find_msigdb(self, query, genes_df, colname, sigs_df):
        if 'genes' not in query: return bytes("")
        # set the limit
        if 'limit' in query:
            print(query['limit'])
            limit = int(query['limit'][0])
        else:
            limit = 100
        q_genes = query['genes']
        scores = self.find_genesets(q_genes, genes_df, colname, msigdb_bin, fout="msigdb")
        scores = scores.flatten()
        sigs_df['overlap'] = scores
        nonzeros = np.sum(scores>0)
        if  nonzeros < limit:
            limit = nonzeros

        if limit == 0:
            return bytes("logResults("+
                json.dumps({
                    'top': []})+")", "utf-8")
        else:
            scores_idx = np.argsort(scores)[-limit:][::-1]
            top = sigs_df.iloc[scores_idx]

            # to compute p_values
            top['DE_genes'] = len(q_genes)
            top['backgrounds'] = genes_df.shape[0]
            top['pvalue'] = top.apply(fisher_exact, axis=1)
            return bytes("logResults("+
                json.dumps({
                    'top': top.to_json(orient='records')})+")", "utf-8")

    def find_kegg(self, query, genes_df, colname, sigs_df):
        if 'genes' in query: q_genes = query['genes']
        elif 'genes[]' in query: q_genes = query['genes[]']
        else: return bytes("Not Found", "utf-8")

        # set the limit
        if 'limit' in query:
            print(query['limit'])
            limit = int(query['limit'][0])
        else:
            limit = 100
        
        scores = self.find_genesets(q_genes, genes_df, colname, kegg_bin, fout="kegg")
        scores = scores.flatten()
        sigs_df['overlap'] = scores
        nonzeros = np.sum(scores>0)
        if  nonzeros < limit:
            limit = nonzeros
        
        if limit == 0:
            scores_idx = []
            top = sigs_df.iloc[scores_idx]
        else:
            scores_idx = np.argsort(scores)[-limit:][::-1]
            top = sigs_df.iloc[scores_idx]
            # to compute p_values
            top['DE_genes'] = len(q_genes)
            top['backgrounds'] = genes_df.shape[0]
            top['pvalue'] = top.apply(fisher_exact, axis=1)
        return bytes("kegg_results("+
            # json.dumps({'top': top.to_json(orient='records')})+")", "utf-8")
            top.to_json(orient='records')+")", "utf-8")

    def find_kegg_ids(self, query, genes_df, src_colname, target_colname):
        if 'genes' in query: q_genes = query['genes']
        elif 'genes[]' in query: q_genes = query['genes[]']
        else: return bytes("Not Found", "utf-8")

        kegg_ids = []
        if len(q_genes) > 0:
            kegg_ids = genes_df[genes_df[src_colname].isin(q_genes)]
        return bytes("keggid_results("+
            # json.dumps({'kegg_ids': kegg_ids.to_json(orient='records')})+")", "utf-8")
            kegg_ids.to_json(orient='records')+")", "utf-8")

    def do_GET(self):
        if self.path == "/":
            self.path = "/index.html"
        if self.path.startswith("/lincs/signature?"):
            print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()

            query = urllib.parse.parse_qs(post_data)
            print(query)

            if '/lincs/signature?sig_id' in query:
                # cid is a list
                cid = query['/lincs/signature?sig_id']
                sig_info = lincs_sigs[lincs_sigs.sig_id==cid[0]]
                if sig_info.shape[0] == 0:
                    self.wfile.write(bytes("logResults("+
                    json.dumps({
                        'sigs': [],
                        'genes': []})+")", "utf-8"))
                    return
                
                sig_index = sig_info.index[0]
                upbits = lincsup[sig_index,:]
                dnbits = lincsdn[sig_index,:]
                up_genes = bin2genes(upbits, lincs_genes.pr_gene_symbol.tolist())
                dn_genes = bin2genes(dnbits, lincs_genes.pr_gene_symbol.tolist())
                up_ezids = bin2genes(upbits, lincs_genes.pr_gene_id.tolist())
                dn_ezids = bin2genes(dnbits, lincs_genes.pr_gene_id.tolist())

                g_info = dict()
                g_info[cid[0]] = {
                    'up': up_genes,
                    'dn': dn_genes,
                    'up_ezids': up_ezids,
                    'dn_ezids': dn_ezids
                }

                self.wfile.write(bytes("logResults("+
                    json.dumps({
                        'sigs': sig_info.to_json(orient='records'),
                        'genes': g_info})+")", "utf-8"))
            return

        if self.path.startswith("/msigdb?"):
            print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()
            query = urllib.parse.parse_qs(post_data)
            print(query)

            self.wfile.write(self.find_msigdb(query, msigdb_genes, "gene_symbol", msigdb_genesets))
            return

        if self.path.startswith("/kegg?"):
            # print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()
            query = urllib.parse.parse_qs(post_data)
            # print(query)

            self.wfile.write(self.find_kegg(query, kegg_genes, "symbol", kegg_genesets))
            return

        if self.path.startswith("/kegg/ids?"):
            # print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()
            query = urllib.parse.parse_qs(post_data)
            # print(query)

            self.wfile.write(self.find_kegg_ids(query, kegg_genes, "symbol", "gene"))
            return

        if self.path.startswith("/lincs/gene_matix?"):
            print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()
            query = urllib.parse.parse_qs(post_data)
            # print(query)

            self.wfile.write(self.get_lincs_gene_matrix(query))
            return

        if self.path.startswith("/lincs/query?"):
            ############################################
            print("incomming http: ", self.path)
            # Gets the size of data
            print(self.headers['Content-Length'])
            post_data = self.path
            # content_length = int(self.headers['Content-Length'])
            # Gets the data itself
            # post_data = self.rfile.read(content_length).decode("utf-8")
            self.send_response(200)
            self.send_header('Content-type', 'application/x-javascript')
            self.end_headers()

            # print(urlparse.parse_qs(post_data))
            query = urllib.parse.parse_qs(post_data)
            self.wfile.write(self.find_lincs(query, lincs_genes, "pr_gene_symbol", lincs_sigs))
            return
            ############################################

        try:
            # Check the file extension required and
            # set the right mime type
            sendReply = False
            if self.path.endswith(".html"):
                mimetype = 'text/html'
                sendReply = True
            if self.path.endswith(".jpg"):
                mimetype = 'image/jpg'
                sendReply = True
            if self.path.endswith(".gif"):
                mimetype = 'image/gif'
                sendReply = True
            if self.path.endswith(".js"):
                mimetype = 'application/javascript'
                sendReply = True
            if self.path.endswith(".css"):
                mimetype = 'text/css'
                sendReply = True
            if self.path.endswith(".txt"):
                mimetype = 'text/tab-separated-values'
                sendReply = True

            if self.path.startswith("/pathway.html"):
                self.path = "/pathway.html"
                mimetype = 'text/html'
                sendReply = True

            if sendReply:
                # Open the static file requested and send it
                f = open(os.curdir + os.sep + self.path)
                self.send_response(200)
                self.send_header('Content-type', mimetype)
                self.end_headers()
                self.wfile.write(f.read().encode())
                f.close()
            return
        except IOError:
            self.send_error(404, 'File Not Found: %s' % self.path)

    #   POST is for submitting data.
    def do_POST(self):
        print("incomming http: ", self.path)
        # Gets the size of data
        content_length = int(self.headers['Content-Length'])
        # Gets the data itself
        post_data = self.rfile.read(content_length).decode("utf-8")
        self.send_response(200)
        # print(urlparse.parse_qs(post_data))
        query = urllib.parse.parse_qs(post_data)

        print(query)

        # set the limit
        if 'limit' in query:
            print(query['limit'])
            limit = int(query['limit'][0])
        else:
            limit = 100
        

        dnidx = upidx = []
        if 'up[]' in query:
            upidx = np.unique(lincs_genes[lincs_genes.pr_gene_symbol.isin(query['up[]'])].index)
        if 'dn[]' in query:
            dnidx = np.unique(lincs_genes[lincs_genes.pr_gene_symbol.isin(query['dn[]'])].index)

        # print(upidx)
        # print(dnidx)

        # print(lincs_signatures.head())
        # print(lincsmeta.head())

        qbindn = np.zeros(length*64, dtype=bool)
        qbinup = np.zeros(length*64, dtype=bool)
        qbindn[dnidx] = True
        qbinup[upidx] = True

        qdn = bool2ints(bool2int, qbindn, length)
        qup = bool2ints(bool2int, qbinup, length)

        qdn = qdn.reshape((-1, length))
        qup = qup.reshape((-1, length))

        print(qup)
        print(qdn)

        bsf.analysis_with_query(lincsup, qup, 'qupup.txt', '')
        bsf.analysis_with_query(lincsdn, qdn, 'qdndn.txt', '')
        bsf.analysis_with_query(lincsup, qdn, 'qupdn.txt', '')
        bsf.analysis_with_query(lincsdn, qup, 'qdnup.txt', '')

        # read a result table (binary file)
        with open("bin_qupup.txt.bin", "rb") as f:
            qupup = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsignatures, 1))
        with open("bin_qdndn.txt.bin", "rb") as f:
            qdndn = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsignatures, 1))
        with open("bin_qupdn.txt.bin", "rb") as f:
            qupdn = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsignatures, 1))
        with open("bin_qdnup.txt.bin", "rb") as f:
            qdnup = np.frombuffer(f.read(), dtype=np.uint32).reshape((nsignatures, 1))
        
        lincs_sigs['upup'] = qupup.flatten()
        lincs_sigs['dndn'] = qdndn.flatten()
        lincs_sigs['updn'] = qupdn.flatten()
        lincs_sigs['dnup'] = qdnup.flatten()
        
        pos_scores = (qupup + qdndn).flatten()
        neg_scores = (qdnup + qupdn).flatten()

        # pos_jar = qupup/np.sum(qbinup) + qdndn/np.sum(qbindn)
        # neg_jar = qdnup/np.sum(qbinup) + qupdn/np.sum(qbindn)

        pos_scores_idx = np.argsort(pos_scores)[-limit:][::-1]
        neg_scores_idx = np.argsort(neg_scores)[-limit:][::-1]

        # pos_jar_idx = np.argsort(pos_jar, axis=0)
        # neg_jar_idx = np.argsort(neg_jar, axis=0)

        print(pos_scores_idx[:5])
        print(pos_scores[pos_scores_idx[:5]])

        print(neg_scores_idx[:5])
        print(neg_scores[neg_scores_idx[:5]])

        pos = lincs_sigs.iloc[pos_scores_idx]
        neg = lincs_sigs.iloc[neg_scores_idx]

        # results
        self.wfile.write(bytes(
            json.dumps({
                'pos': pos.to_json(orient='records'),
                'neg': neg.to_json(orient='records')}), "utf-8"))
        # self.wfile.write(bytes(
        #     json.dumps({
        #         'upup': qupup.flatten().tolist(),
        #         'dndn': qdndn.flatten().tolist(),
        #         'updn': qupdn.flatten().tolist(),
        #         'dnup': qdnup.flatten().tolist()}), "utf-8"))



def run(server_class=HTTPServer, handler_class=BSFServer, port=80):
    server_address = ('', port)
    httpd = server_class(server_address, handler_class)
    print(time.asctime(), "Server Starts - %s:%s" % server_address)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        # f_lincs.close()
        pass
    httpd.server_close()
    print(time.asctime(), "Server Stops - %s:%s" % server_address)


if __name__ == "__main__":
    from sys import argv
    if len(argv) == 2:
        run(port=int(argv[1]))
    else:
        run()

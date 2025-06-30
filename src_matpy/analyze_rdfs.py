import numpy as np
import os
import glob

def ana_rdf(simdir,strpref,colid):

    rdf_file = find_latest_file(glob.glob(simdir + '/' + strpref + '*'))
    if rdf_file == -1:
        return -1, -1
    
    rarr,grarr = extract_gofr(rdf_file,colid)
    nrarr      = compute_nofr(rarr, grarr)
    return rarr,grarr,nrarr

def find_latest_file(rdf_files):
    if not len(rdf_files):
        print('ERROR: No RDF files found in ' + simdir)
        return -1
    rdf_file = max(rdf_files, key=os.path.getmtime)
    return rdf_file
    
        
def extract_gofr(filename,colid):
    data = np.loadtxt(filename, skiprows=1)  # Skips the first header line
    col1 = data[:, 0]  # Column 1
    col2 = data[:, colid]  # Column n
    return col1, col2

def compute_nofr(r, g_r):
    dr = np.gradient(r)  # Compute spacing (dr) between column1 values
    integral_values = 4 * np.pi * r**2 * g_r * dr   # Element-wise integral contribution
    n_r = np.cumsum(integral_values)  # Compute cumulative sum (integral)
    return n_r

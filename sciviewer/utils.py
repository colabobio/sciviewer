import numpy as np
import numpy.linalg as la

def angle_between(v1, v2):
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def load_data(data, expr, umap, gene_names, cell_names,
              embedding_name, pearsonsThreshold, tThreshold,
              pvalueThreshold, maxdisplaygenes_pos, maxdisplaygenes_neg, use_raw):
    data.load(expr, umap, gene_names, cell_names,
              embedding_name, pearsonsThreshold, tThreshold,
              pvalueThreshold, maxdisplaygenes_pos, maxdisplaygenes_neg, use_raw)
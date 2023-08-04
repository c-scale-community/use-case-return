import numpy as np
from multiprocessing import Pool, RawArray
from osgeo import gdal

PROC_OBJS = {}

def read_init(fps, sm):
    """ Helper method for setting the entries of global variable `PROC_OBJS` to be available during multiprocessing. """
    PROC_OBJS['filepaths'] = fps
    PROC_OBJS['shm_map'] = sm

def read_single_file(file_idx):
    """
    Function being responsible to read data from a single GeoTIFF file and check if it is broken or not.
    This function is meant to be executed in parallel on different cores.

    Parameters
    ----------
    file_idx : any
        Index value to access a specific row of the file register. The actual value should come from the Pool's
        mapping function.

    """
    gdal.UseExceptions()
    filepaths = PROC_OBJS['filepaths']
    shm_map = PROC_OBJS['shm_map']
    filepath = filepaths[file_idx]

    try:
        with GeoTiffFile(filepath, mode='r') as gt_file:
            _ = gt_file.read(return_tags=False)
    except:
        shm_data = np.frombuffer(shm_map, dtype=np.uint8)
        shm_data[file_idx] = 1

def remove_broken_files(inventory, n_cores=1):
	"""
	Removes broken GeoTIFF files from internal inventory if they can't be opened.

	Parameters
	----------
	n_cores : int, optional
		Cores used when checking for broken files (defaults to 1).

	"""
	filepaths = inventory['filepath']
	n_files = len(filepaths)
	data_nshm = np.zeros((n_files), dtype=np.uint8)
	c_dtype = np.ctypeslib.as_ctypes_type(data_nshm.dtype)
	shm_rar = RawArray(c_dtype, data_nshm.size)
	shm_data = np.frombuffer(shm_rar, dtype=np.uint8)
	shm_data[:] = data_nshm[:]
	shm_map = shm_rar
	file_idxs = list(range(0, n_files))
	with Pool(n_cores, initializer=read_init, initargs=(filepaths, shm_map)) as p:
		p.map(read_single_file, file_idxs)

	shm_data = np.frombuffer(shm_map, dtype=np.uint8)
	idxs2rem = inventory.index[shm_data.astype(bool)]
	return inventory.drop(index=idxs2rem)
	
	
if __name__ == "__main__":
	# ....
	remove_broken_files(...)
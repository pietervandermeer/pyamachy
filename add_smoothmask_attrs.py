import h5py

f = h5py.File("sdmf_smooth_pyxelmask30.h5","r+")
f.attrs["sdmfVersion"] = "3.2.0"
f.attrs["swVersion"] = "1.0.0"
f.attrs["calibVersion"] = "1.0.0"
f.attrs["dbVersion"] = "1.0.0"
f.close()

